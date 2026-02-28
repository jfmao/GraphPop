package org.graphpop.procedures;

import org.graphpop.VectorOps;
import org.neo4j.graphdb.Label;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Transaction;
import org.neo4j.procedure.Context;
import org.neo4j.procedure.Mode;
import org.neo4j.procedure.Name;
import org.neo4j.procedure.Procedure;

import java.util.Map;
import java.util.stream.Stream;

/**
 * Computes population divergence statistics (Fst, Dxy, Da) between two populations.
 *
 * <p>FAST PATH: uses only allele-count arrays on Variant nodes.</p>
 *
 * <p>Usage:
 * <pre>
 * CALL graphpop.divergence('chr22', 16000000, 17000000, 'AFR', 'EUR')
 * YIELD fst_hudson, dxy, da, n_variants
 * </pre>
 * </p>
 */
public class DivergenceProcedure {

    @Context
    public Transaction tx;

    public static class DivergenceResult {
        public double fst_hudson;
        public double dxy;
        public double da;
        public long n_variants;

        public DivergenceResult() {}
    }

    @Procedure(name = "graphpop.divergence", mode = Mode.READ)
    public Stream<DivergenceResult> divergence(
            @Name("chr") String chr,
            @Name("start") long start,
            @Name("end") long end,
            @Name("pop1") String pop1,
            @Name("pop2") String pop2
    ) {
        int popIndex1 = -1;
        int popIndex2 = -1;

        double fstNumeratorSum = 0.0;
        double fstDenominatorSum = 0.0;
        double dxySum = 0.0;
        double piWithin1Sum = 0.0;
        double piWithin2Sum = 0.0;
        long nVariants = 0;

        var result = tx.execute(
                "MATCH (v:Variant) WHERE v.chr = $chr AND v.pos >= $start AND v.pos <= $end RETURN v",
                Map.of("chr", chr, "start", start, "end", end)
        );

        try {
            while (result.hasNext()) {
                Map<String, Object> row = result.next();
                Node variant = (Node) row.get("v");

                // Resolve population indices on first variant
                if (popIndex1 < 0) {
                    String[] popIds = ArrayUtil.toStringArray(variant.getProperty("pop_ids"));
                    for (int i = 0; i < popIds.length; i++) {
                        if (popIds[i].equals(pop1)) popIndex1 = i;
                        if (popIds[i].equals(pop2)) popIndex2 = i;
                    }
                    if (popIndex1 < 0) throw new RuntimeException("Population not found: " + pop1);
                    if (popIndex2 < 0) throw new RuntimeException("Population not found: " + pop2);
                }

                int[] anArr = ArrayUtil.toIntArray(variant.getProperty("an"));
                double[] afArr = ArrayUtil.toDoubleArray(variant.getProperty("af"));

                int an1 = anArr[popIndex1];
                int an2 = anArr[popIndex2];
                double af1 = afArr[popIndex1];
                double af2 = afArr[popIndex2];

                if (an1 < 2 || an2 < 2) continue;

                nVariants++;

                // Hudson's Fst components
                double[] fstComp = VectorOps.hudsonFstComponents(af1, an1, af2, an2);
                fstNumeratorSum += fstComp[0];
                fstDenominatorSum += fstComp[1];

                // Dxy
                dxySum += VectorOps.dxyPerSite(af1, af2);

                // Within-population pi (for Da)
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
        out.dxy = dxySum / L;
        // Da = Dxy - (pi_within_1 + pi_within_2) / 2
        out.da = out.dxy - (piWithin1Sum / L + piWithin2Sum / L) / 2.0;

        return Stream.of(out);
    }
}
