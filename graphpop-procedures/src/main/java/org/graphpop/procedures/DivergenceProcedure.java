package org.graphpop.procedures;

import org.graphpop.VectorOps;
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
 * Computes population divergence statistics (Fst, Dxy, Da) between two populations.
 *
 * <p>FAST PATH: uses only allele-count arrays on Variant nodes.</p>
 *
 * <p>Usage:
 * <pre>
 * CALL graphpop.divergence('chr22', 16000000, 17000000, 'AFR', 'EUR', {min_af: 0.05})
 * YIELD fst_hudson, dxy, da, n_variants
 * </pre>
 * </p>
 */
public class DivergenceProcedure {

    @Context
    public Transaction tx;

    public static class DivergenceResult {
        public double fst_hudson;
        public double fst_wc;
        public double dxy;
        public double da;
        public double pbs;
        public long n_variants;

        public DivergenceResult() {}
    }

    @Procedure(name = "graphpop.divergence", mode = Mode.READ)
    @Description("Compute Fst (Hudson and Weir-Cockerham), Dxy, and PBS divergence statistics between populations in a genomic region.")
    @SuppressWarnings("unchecked")
    public Stream<DivergenceResult> divergence(
            @Name("chr") String chr,
            @Name("start") long start,
            @Name("end") long end,
            @Name("pop1") String pop1,
            @Name("pop2") String pop2,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        VariantFilter filter = VariantFilter.fromOptions(options);
        List<String> samples1 = options != null ? (List<String>) options.get("samples1") : null;
        List<String> samples2 = options != null ? (List<String>) options.get("samples2") : null;
        String pop3 = options != null ? (String) options.get("pop3") : null;
        boolean useSubset1 = samples1 != null && !samples1.isEmpty();
        boolean useSubset2 = samples2 != null && !samples2.isEmpty();

        Map<String, Integer> subsetIndex1 = useSubset1 ? GenotypeLoader.buildSampleIndex(tx, samples1) : null;
        Map<String, Integer> subsetIndex2 = useSubset2 ? GenotypeLoader.buildSampleIndex(tx, samples2) : null;
        int[] packedIndices1 = useSubset1 ? GenotypeLoader.buildPackedIndices(tx, subsetIndex1) : null;
        int[] packedIndices2 = useSubset2 ? GenotypeLoader.buildPackedIndices(tx, subsetIndex2) : null;

        int popIndex1 = -1;
        int popIndex2 = -1;
        int popIndex3 = -1;

        // Pop1 vs Pop2 accumulators
        double fstNumeratorSum = 0.0;
        double fstDenominatorSum = 0.0;
        double wcA12 = 0.0, wcB12 = 0.0, wcC12 = 0.0;
        double dxySum = 0.0;
        double piWithin1Sum = 0.0;
        double piWithin2Sum = 0.0;
        // PBS accumulators: pop1 vs pop3, pop2 vs pop3
        double wcA13 = 0.0, wcB13 = 0.0, wcC13 = 0.0;
        double wcA23 = 0.0, wcB23 = 0.0, wcC23 = 0.0;
        long nVariants = 0;

        var result = tx.execute(
                VariantQuery.build(options),
                VariantQuery.params(options, Map.of("chr", chr, "start", start, "end", end))
        );

        try {
            while (result.hasNext()) {
                Map<String, Object> row = result.next();
                Node variant = (Node) row.get("v");

                int ac1, an1, ac2, an2, het1, het2, homAlt1;
                double af1, af2;

                if (useSubset1) {
                    SampleSubsetComputer.SubsetStats ss = SampleSubsetComputer.compute(variant, subsetIndex1, packedIndices1);
                    ac1 = ss.ac; an1 = ss.an; af1 = ss.af; het1 = ss.hetCount; homAlt1 = ss.homAltCount;
                } else {
                    if (popIndex1 < 0) {
                        String[] popIds = ArrayUtil.toStringArray(variant.getProperty("pop_ids"));
                        for (int i = 0; i < popIds.length; i++) {
                            if (popIds[i].equals(pop1)) popIndex1 = i;
                            if (popIds[i].equals(pop2)) popIndex2 = i;
                            if (pop3 != null && popIds[i].equals(pop3)) popIndex3 = i;
                        }
                        if (popIndex1 < 0) throw new RuntimeException("Population not found: " + pop1);
                        if (!useSubset2 && popIndex2 < 0) throw new RuntimeException("Population not found: " + pop2);
                        if (pop3 != null && popIndex3 < 0) throw new RuntimeException("Population not found: " + pop3);
                    }
                    int[] acArr = ArrayUtil.toIntArray(variant.getProperty("ac"));
                    int[] anArr = ArrayUtil.toIntArray(variant.getProperty("an"));
                    double[] afArr = ArrayUtil.toDoubleArray(variant.getProperty("af"));
                    int[] hetArr = ArrayUtil.toIntArray(variant.getProperty("het_count"));
                    int[] homAltArr = ArrayUtil.toIntArray(variant.getProperty("hom_alt_count"));
                    ac1 = acArr[popIndex1]; an1 = anArr[popIndex1]; af1 = afArr[popIndex1];
                    het1 = hetArr[popIndex1]; homAlt1 = homAltArr[popIndex1];
                }

                if (useSubset2) {
                    SampleSubsetComputer.SubsetStats ss = SampleSubsetComputer.compute(variant, subsetIndex2, packedIndices2);
                    ac2 = ss.ac; an2 = ss.an; af2 = ss.af; het2 = ss.hetCount;
                } else {
                    int[] acArr = ArrayUtil.toIntArray(variant.getProperty("ac"));
                    int[] anArr = ArrayUtil.toIntArray(variant.getProperty("an"));
                    double[] afArr = ArrayUtil.toDoubleArray(variant.getProperty("af"));
                    int[] hetArr = ArrayUtil.toIntArray(variant.getProperty("het_count"));
                    ac2 = acArr[popIndex2]; an2 = anArr[popIndex2]; af2 = afArr[popIndex2];
                    het2 = hetArr[popIndex2];
                }

                if (an1 < 2 || an2 < 2) continue;

                if (filter.isActive() && !filter.passes(ac1, an1, af1, het1, homAlt1, variant)) continue;

                // Pop3 data (only when PBS requested)
                int ac3 = 0, an3 = 0, het3 = 0;
                if (pop3 != null) {
                    int[] acArr = ArrayUtil.toIntArray(variant.getProperty("ac"));
                    int[] anArr = ArrayUtil.toIntArray(variant.getProperty("an"));
                    int[] hetArr = ArrayUtil.toIntArray(variant.getProperty("het_count"));
                    ac3 = acArr[popIndex3]; an3 = anArr[popIndex3]; het3 = hetArr[popIndex3];
                    if (an3 < 2) continue;
                }

                nVariants++;

                double[] fstComp = VectorOps.hudsonFstComponents(af1, an1, af2, an2);
                fstNumeratorSum += fstComp[0];
                fstDenominatorSum += fstComp[1];

                double[] wcComp = VectorOps.wcFstComponents(ac1, an1, het1, ac2, an2, het2);
                wcA12 += wcComp[0];
                wcB12 += wcComp[1];
                wcC12 += wcComp[2];

                if (pop3 != null) {
                    double[] wc13 = VectorOps.wcFstComponents(ac1, an1, het1, ac3, an3, het3);
                    wcA13 += wc13[0]; wcB13 += wc13[1]; wcC13 += wc13[2];
                    double[] wc23 = VectorOps.wcFstComponents(ac2, an2, het2, ac3, an3, het3);
                    wcA23 += wc23[0]; wcB23 += wc23[1]; wcC23 += wc23[2];
                }

                dxySum += VectorOps.dxyPerSite(af1, af2);

                piWithin1Sum += VectorOps.piPerSite(af1, an1);
                piWithin2Sum += VectorOps.piPerSite(af2, an2);
            }
        } finally {
            result.close();
        }

        DivergenceResult out = new DivergenceResult();
        out.n_variants = nVariants;

        if (nVariants == 0) {
            return Stream.of(out);
        }

        double L = nVariants;
        out.fst_hudson = fstDenominatorSum > 0 ? fstNumeratorSum / fstDenominatorSum : 0.0;
        double wcDenom12 = wcA12 + wcB12 + wcC12;
        out.fst_wc = wcDenom12 > 0 ? wcA12 / wcDenom12 : 0.0;
        out.dxy = dxySum / L;
        out.da = out.dxy - (piWithin1Sum / L + piWithin2Sum / L) / 2.0;

        // PBS: Population Branch Statistic (Yi et al. 2010)
        if (pop3 != null) {
            double fst12 = out.fst_wc;
            double wcDenom13 = wcA13 + wcB13 + wcC13;
            double fst13 = wcDenom13 > 0 ? wcA13 / wcDenom13 : 0.0;
            double wcDenom23 = wcA23 + wcB23 + wcC23;
            double fst23 = wcDenom23 > 0 ? wcA23 / wcDenom23 : 0.0;

            // T = -log(1 - Fst), clamp Fst to avoid log(0)
            double t12 = -Math.log(Math.max(1e-10, 1.0 - fst12));
            double t13 = -Math.log(Math.max(1e-10, 1.0 - fst13));
            double t23 = -Math.log(Math.max(1e-10, 1.0 - fst23));

            // PBS for pop1 (focal population)
            out.pbs = (t12 + t13 - t23) / 2.0;
        }

        return Stream.of(out);
    }
}
