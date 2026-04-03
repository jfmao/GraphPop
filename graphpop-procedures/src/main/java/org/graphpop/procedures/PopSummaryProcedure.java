package org.graphpop.procedures;

import org.graphpop.VectorOps;
import org.neo4j.graphdb.Label;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Transaction;
import org.neo4j.procedure.Context;
import org.neo4j.procedure.Description;
import org.neo4j.procedure.Mode;
import org.neo4j.procedure.Name;
import org.neo4j.procedure.Procedure;

import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

/**
 * Whole-chromosome population summary statistics.
 *
 * <p>Computes diversity, neutrality tests, and heterozygosity across an entire
 * chromosome for a population, and writes results as properties on the Population node.</p>
 *
 * <p>Usage:
 * <pre>
 * CALL graphpop.pop_summary('chr22', 'EUR', {variant_type: 'SNP'})
 * YIELD pi, theta_w, tajima_d, fay_wu_h, mean_he, mean_ho, mean_fis,
 *       n_variants, n_segregating, n_polarized
 * </pre>
 * </p>
 */
public class PopSummaryProcedure {

    @Context
    public Transaction tx;

    public static class PopSummaryResult {
        public double pi;
        public double theta_w;
        public double tajima_d;
        public double fay_wu_h;
        public double fay_wu_h_norm;
        public double mean_he;
        public double mean_ho;
        public double mean_fis;
        public long n_variants;
        public long n_segregating;
        public long n_polarized;

        public PopSummaryResult() {}
    }

    @Procedure(name = "graphpop.pop_summary", mode = Mode.WRITE)
    @Description("Compute population-level summary statistics including sample counts, mean heterozygosity, and allele frequency distributions.")
    @SuppressWarnings("unchecked")
    public Stream<PopSummaryResult> popSummary(
            @Name("chr") String chr,
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

        // Query all variants on the chromosome
        var result = tx.execute(
                VariantQuery.buildChromosome(options),
                VariantQuery.params(options, Map.of("chr", chr))
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

                // Fay & Wu's H
                Object ancestralObj = variant.getProperty("ancestral_allele", null);
                if (ancestralObj != null) {
                    String ancestral = (String) ancestralObj;
                    int derivedCount;
                    if ("REF".equals(ancestral)) {
                        derivedCount = ac;
                    } else if ("ALT".equals(ancestral)) {
                        derivedCount = an - ac;
                    } else {
                        continue;
                    }
                    nPolarized++;
                    thetaHSum += FayWuH.thetaHPerSite(derivedCount, an);
                    piPolarizedSum += VectorOps.piPerSite(af, an);
                }
            }
        } finally {
            result.close();
        }

        PopSummaryResult out = new PopSummaryResult();
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
            n = firstAN > 0 ? firstAN : 2 * subsetIndex.size();
            a_n = harmonicNumber(n - 1);
            a_n2 = harmonicNumber2(n - 1);
        } else {
            if (ctx == null) return Stream.of(out);
            n = firstAN > 0 ? firstAN : 2 * ctx.nSamples;
            a_n = harmonicNumber(n - 1);
            a_n2 = harmonicNumber2(n - 1);
        }

        double L = nVariants;
        out.pi = piSum / L;
        out.mean_he = heSum / L;
        out.mean_ho = hoSum / L;
        out.mean_fis = out.mean_he > 0 ? 1.0 - out.mean_ho / out.mean_he : 0.0;
        out.theta_w = a_n > 0 ? (nSegregating / a_n) / L : 0.0;

        out.tajima_d = TajimaD.compute(piSum, nSegregating, n, a_n, a_n2);

        if (nPolarized > 0) {
            out.fay_wu_h = FayWuH.compute(piPolarizedSum, thetaHSum) / nPolarized;
            out.fay_wu_h_norm = FayWuH.normalizedH(piPolarizedSum, thetaHSum, nPolarized, n, a_n);
        }

        // Write summary to Population node
        Node popNode = tx.findNode(Label.label("Population"), "populationId", pop);
        if (popNode != null) {
            popNode.setProperty("summary_pi_" + chr, out.pi);
            popNode.setProperty("summary_theta_w_" + chr, out.theta_w);
            popNode.setProperty("summary_tajima_d_" + chr, out.tajima_d);
            popNode.setProperty("summary_fay_wu_h_" + chr, out.fay_wu_h);
            popNode.setProperty("summary_n_variants_" + chr, out.n_variants);
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
