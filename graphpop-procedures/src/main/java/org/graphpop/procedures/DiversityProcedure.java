package org.graphpop.procedures;

import org.graphpop.VectorOps;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Transaction;
import org.neo4j.procedure.Context;
import org.neo4j.procedure.Mode;
import org.neo4j.procedure.Name;
import org.neo4j.procedure.Procedure;

import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

/**
 * Computes diversity statistics for a genomic region and population.
 *
 * <p>FAST PATH: uses only the allele-count arrays stored on Variant nodes,
 * so complexity is O(V x K) independent of sample size.</p>
 *
 * <p>When {@code samples} option is provided, reads packed genotype arrays
 * for on-the-fly recomputation (slower but fully flexible).</p>
 *
 * <p>Usage:
 * <pre>
 * CALL graphpop.diversity('chr22', 16000000, 17000000, 'AFR')
 * YIELD pi, theta_w, tajima_d, fay_wu_h, het_exp, het_obs, fis, n_variants, n_segregating
 * </pre>
 * </p>
 */
public class DiversityProcedure {

    @Context
    public Transaction tx;

    public static class DiversityResult {
        public double pi;
        public double theta_w;
        public double tajima_d;
        public double fay_wu_h;
        public double fay_wu_h_norm;
        public double het_exp;
        public double het_obs;
        public double fis;
        public long n_variants;
        public long n_segregating;
        public long n_polarized;

        public DiversityResult() {}
    }

    @Procedure(name = "graphpop.diversity", mode = Mode.READ)
    @SuppressWarnings("unchecked")
    public Stream<DiversityResult> diversity(
            @Name("chr") String chr,
            @Name("start") long start,
            @Name("end") long end,
            @Name("pop") String pop,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        VariantFilter filter = VariantFilter.fromOptions(options);
        List<String> sampleList = options != null ? (List<String>) options.get("samples") : null;
        boolean useSubset = sampleList != null && !sampleList.isEmpty();

        Map<String, Integer> subsetIndex = null;
        int[] packedIndices = null;
        if (useSubset) {
            subsetIndex = GenotypeLoader.buildSampleIndex(tx, sampleList);
            packedIndices = GenotypeLoader.buildPackedIndices(tx, subsetIndex);
        }

        PopulationContext ctx = null;

        double piSum = 0.0;
        double heSum = 0.0;
        double hoSum = 0.0;
        double thetaHSum = 0.0;
        double piPolarizedSum = 0.0;
        long nVariants = 0;
        long nSegregating = 0;
        long nPolarized = 0;
        int firstAN = -1; // track first observed AN for ploidy-aware n derivation

        var result = tx.execute(
                VariantQuery.build(options),
                VariantQuery.params(options, Map.of("chr", chr, "start", start, "end", end))
        );

        try {
            while (result.hasNext()) {
                Map<String, Object> row = result.next();
                Node variant = (Node) row.get("v");

                int ac, an, hetCount, homAltCount;
                double af;

                if (useSubset) {
                    SampleSubsetComputer.SubsetStats ss =
                            SampleSubsetComputer.compute(variant, subsetIndex, packedIndices);
                    ac = ss.ac; an = ss.an; af = ss.af;
                    hetCount = ss.hetCount; homAltCount = ss.homAltCount;
                } else {
                    if (ctx == null) {
                        ctx = PopulationContext.resolve(tx, variant, pop);
                    }
                    int[] acArr = ArrayUtil.toIntArray(variant.getProperty("ac"));
                    int[] anArr = ArrayUtil.toIntArray(variant.getProperty("an"));
                    double[] afArr = ArrayUtil.toDoubleArray(variant.getProperty("af"));
                    int[] hetArr = ArrayUtil.toIntArray(variant.getProperty("het_count"));
                    int[] homAltArr = ArrayUtil.toIntArray(variant.getProperty("hom_alt_count"));

                    ac = acArr[ctx.index]; an = anArr[ctx.index]; af = afArr[ctx.index];
                    hetCount = hetArr[ctx.index]; homAltCount = homAltArr[ctx.index];
                }

                if (an < 2) continue;

                if (filter.isActive() && !filter.passes(ac, an, af, hetCount, homAltCount, variant)) {
                    continue;
                }

                if (firstAN < 0) firstAN = an;

                nVariants++;
                piSum += VectorOps.piPerSite(af, an);
                heSum += 2.0 * af * (1.0 - af);

                // H_obs: hetCount / nDiploid (approximate for mixed ploidy)
                int nDiploid = an / 2;
                if (nDiploid > 0) {
                    hoSum += (double) hetCount / nDiploid;
                }

                if (ac > 0 && ac < an) {
                    nSegregating++;
                }

                // Fay & Wu's H: requires ancestral allele
                Object ancestralObj = variant.getProperty("ancestral_allele", null);
                if (ancestralObj != null) {
                    String ancestral = (String) ancestralObj;
                    int derivedCount;
                    if ("REF".equals(ancestral)) {
                        derivedCount = ac;
                    } else if ("ALT".equals(ancestral)) {
                        derivedCount = an - ac;
                    } else {
                        continue; // skip unknown polarization within this block
                    }
                    nPolarized++;
                    thetaHSum += FayWuH.thetaHPerSite(derivedCount, an);
                    piPolarizedSum += VectorOps.piPerSite(af, an);
                }
            }
        } finally {
            result.close();
        }

        DiversityResult out = new DiversityResult();
        out.n_variants = nVariants;
        out.n_segregating = nSegregating;
        out.n_polarized = nPolarized;

        if (nVariants == 0) {
            return Stream.of(out);
        }

        // Derive n from actual AN (ploidy-aware: works for haploid/mixed chromosomes)
        int n;
        double a_n, a_n2;
        if (useSubset) {
            // For custom subsets, use first observed AN (accounts for ploidy)
            n = firstAN > 0 ? firstAN : 2 * subsetIndex.size();
            a_n = harmonicNumber(n - 1);
            a_n2 = harmonicNumber2(n - 1);
        } else {
            if (ctx == null) return Stream.of(out);
            // Use first observed AN instead of assuming 2*nSamples (handles non-diploid chr)
            n = firstAN > 0 ? firstAN : 2 * ctx.nSamples;
            a_n = harmonicNumber(n - 1);
            a_n2 = harmonicNumber2(n - 1);
        }

        double L = nVariants;
        out.pi = piSum / L;
        out.het_exp = heSum / L;
        out.het_obs = hoSum / L;
        out.fis = out.het_exp > 0 ? 1.0 - out.het_obs / out.het_exp : 0.0;
        out.theta_w = a_n > 0 ? (nSegregating / a_n) / L : 0.0;

        out.tajima_d = TajimaD.compute(piSum, nSegregating, n, a_n, a_n2);

        // Fay & Wu's H (only if we have polarized sites)
        if (nPolarized > 0) {
            out.fay_wu_h = FayWuH.compute(piPolarizedSum, thetaHSum) / nPolarized;
            out.fay_wu_h_norm = FayWuH.normalizedH(piPolarizedSum, thetaHSum, nPolarized, n, a_n);
        }

        return Stream.of(out);
    }

    private static double harmonicNumber(int m) {
        double h = 0.0;
        for (int i = 1; i <= m; i++) h += 1.0 / i;
        return h;
    }

    private static double harmonicNumber2(int m) {
        double h = 0.0;
        for (int i = 1; i <= m; i++) h += 1.0 / ((double) i * i);
        return h;
    }
}
