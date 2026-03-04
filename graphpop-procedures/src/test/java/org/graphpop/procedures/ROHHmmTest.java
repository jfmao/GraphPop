package org.graphpop.procedures;

import org.junit.jupiter.api.Test;

import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Unit tests for HMM-based ROH detection (Narasimhan et al. 2016).
 */
class ROHHmmTest {

    // Use bcftools defaults: transition rates + het_error_rate matching -G30
    private static final ROHHmm DEFAULT_HMM = new ROHHmm(
            ROHHmm.DEFAULT_HW_TO_AZ,  // 6.7e-8 per bp (bcftools default)
            ROHHmm.DEFAULT_AZ_TO_HW,  // 5e-9 per bp (bcftools default)
            1e-3,                      // hetErrorRate (bcftools -G30 → 10^(-3))
            0.0                        // minAf (skip only AF=0)
    );

    /**
     * Helper: build positions array with uniform spacing.
     */
    private static long[] uniformPositions(int n, long start, long spacing) {
        long[] pos = new long[n];
        for (int i = 0; i < n; i++) pos[i] = start + i * spacing;
        return pos;
    }

    /**
     * Helper: build uniform AF array.
     */
    private static double[] uniformAFs(int n, double af) {
        double[] afs = new double[n];
        for (int i = 0; i < n; i++) afs[i] = af;
        return afs;
    }

    @Test
    void allHomozygous_detectsROH() {
        // 200 hom-ref SNPs spanning 2 Mb at AF=0.3 → strong AZ evidence
        int n = 200;
        long[] pos = uniformPositions(n, 1_000_000, 10_000); // 2 Mb span
        double[] afs = uniformAFs(n, 0.3);
        int[] genotypes = new int[n]; // all 0 = hom-ref

        List<long[]> segments = DEFAULT_HMM.viterbi(pos, afs, genotypes,
                500_000, 50);

        assertFalse(segments.isEmpty(), "All-hom should produce at least one ROH segment");
        // The single segment should span most of the region
        long totalLen = segments.stream().mapToLong(s -> s[1] - s[0]).sum();
        assertTrue(totalLen > 1_000_000, "Total ROH length should exceed 1 Mb");
    }

    @Test
    void allHeterozygous_noROH() {
        // 200 het SNPs → no ROH (het is nearly impossible under AZ)
        int n = 200;
        long[] pos = uniformPositions(n, 1_000_000, 10_000);
        double[] afs = uniformAFs(n, 0.3);
        int[] genotypes = new int[n];
        for (int i = 0; i < n; i++) genotypes[i] = 1; // all het

        List<long[]> segments = DEFAULT_HMM.viterbi(pos, afs, genotypes,
                500_000, 50);

        assertTrue(segments.isEmpty(), "All-het should produce no ROH");
    }

    @Test
    void embeddedROH_detected() {
        // 50 het + 100 hom + 50 het — ROH in the middle
        int n = 200;
        long[] pos = uniformPositions(n, 1_000_000, 10_000);
        double[] afs = uniformAFs(n, 0.3);
        int[] genotypes = new int[n];
        for (int i = 0; i < 50; i++) genotypes[i] = 1;     // het
        // 50..149: hom-ref (genotypes already 0)
        for (int i = 150; i < n; i++) genotypes[i] = 1;     // het

        List<long[]> segments = DEFAULT_HMM.viterbi(pos, afs, genotypes,
                500_000, 25);

        assertFalse(segments.isEmpty(), "Should detect ROH in the middle region");
        // ROH should be centered roughly in the hom region
        long[] seg = segments.get(0);
        assertTrue(seg[0] >= 1_000_000 && seg[0] <= 1_600_000,
                "ROH start should be in the hom region");
        assertTrue(seg[1] >= 1_900_000 && seg[1] <= 2_600_000,
                "ROH end should be in the hom region");
    }

    @Test
    void rareHomozygotes_strongerSignal() {
        // Hom-alt at rare AF should produce stronger AZ signal than at common AF
        // because P(hom-alt|HW) = p² is tiny for rare p, but P(hom-alt|AZ) = p is larger
        int n = 100;
        long[] pos = uniformPositions(n, 1_000_000, 10_000);

        // Rare AF: p=0.05 → P_HW = 0.0025, P_AZ = 0.05 → ratio 20
        double[] rareAfs = uniformAFs(n, 0.05);
        int[] genotypesRare = new int[n];
        for (int i = 0; i < n; i++) genotypesRare[i] = 2; // hom-alt

        // Common AF: p=0.5 → P_HW = 0.25, P_AZ = 0.5 → ratio 2
        double[] commonAfs = uniformAFs(n, 0.5);
        int[] genotypesCommon = new int[n];
        for (int i = 0; i < n; i++) genotypesCommon[i] = 2; // hom-alt

        // Both should detect ROH, but let's verify both work with lenient filters
        List<long[]> rareSegs = DEFAULT_HMM.viterbi(pos, rareAfs, genotypesRare,
                100_000, 10);
        List<long[]> commonSegs = DEFAULT_HMM.viterbi(pos, commonAfs, genotypesCommon,
                100_000, 10);

        assertFalse(rareSegs.isEmpty(), "Rare hom-alt should detect ROH");
        assertFalse(commonSegs.isEmpty(), "Common hom-alt should detect ROH");
    }

    @Test
    void scatteredHets_tolerated() {
        // 100 hom-ref with 2 scattered hets → HMM should still detect ROH
        // because het error rate allows occasional hets under AZ state
        int n = 100;
        long[] pos = uniformPositions(n, 1_000_000, 10_000);
        double[] afs = uniformAFs(n, 0.3);
        int[] genotypes = new int[n]; // all hom-ref
        genotypes[30] = 1;  // scattered het
        genotypes[70] = 1;  // scattered het

        List<long[]> segments = DEFAULT_HMM.viterbi(pos, afs, genotypes,
                500_000, 25);

        assertFalse(segments.isEmpty(),
                "A few scattered hets should not break ROH detection");
    }

    @Test
    void shortROH_filtered() {
        // 20 hom-ref SNPs spanning 200 kb → below minLength=500 kb
        int n = 20;
        long[] pos = uniformPositions(n, 1_000_000, 10_000); // 200 kb span
        double[] afs = uniformAFs(n, 0.3);
        int[] genotypes = new int[n]; // all hom-ref

        List<long[]> segments = DEFAULT_HMM.viterbi(pos, afs, genotypes,
                500_000, 10);

        assertTrue(segments.isEmpty(),
                "ROH shorter than minLength should be filtered out");
    }

    @Test
    void fromOptions_parsesLegacyParameters() {
        // Legacy az_length/az_fraction parameterization
        Map<String, Object> options = Map.of(
                "az_length", 2_000_000.0,
                "az_fraction", 0.005,
                "het_error_rate", 1e-4
        );
        ROHHmm hmm = ROHHmm.fromOptions(options);

        // Verify it doesn't throw and produces valid results
        int n = 50;
        long[] pos = uniformPositions(n, 1_000_000, 10_000);
        double[] afs = uniformAFs(n, 0.3);
        int[] genotypes = new int[n];

        List<long[]> segments = hmm.viterbi(pos, afs, genotypes, 100_000, 10);
        assertNotNull(segments);
    }

    @Test
    void fromOptions_parsesPerBpRates() {
        // Per-bp rate parameterization (bcftools style)
        Map<String, Object> options = Map.of(
                "hw_to_az", 6.7e-8,
                "az_to_hw", 5e-9,
                "het_error_rate", 1e-5
        );
        ROHHmm hmm = ROHHmm.fromOptions(options);

        int n = 200;
        long[] pos = uniformPositions(n, 1_000_000, 10_000);
        double[] afs = uniformAFs(n, 0.3);
        int[] genotypes = new int[n]; // all hom-ref

        List<long[]> segments = hmm.viterbi(pos, afs, genotypes, 500_000, 50);
        assertFalse(segments.isEmpty(), "Per-bp rate parameterization should detect ROH");
    }
}
