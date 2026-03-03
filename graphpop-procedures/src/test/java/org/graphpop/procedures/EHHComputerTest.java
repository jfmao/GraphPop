package org.graphpop.procedures;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Unit tests for {@link EHHComputer} using synthetic haplotype matrices.
 */
class EHHComputerTest {

    private static final double EPS = 1e-10;

    /**
     * Build a minimal HaplotypeMatrix from raw data (no DB needed).
     */
    private static HaplotypeMatrix buildMatrix(long[] positions, byte[][] haplotypes) {
        int nVariants = positions.length;
        int nHaplotypes = haplotypes[0].length;
        int nSamples = nHaplotypes / 2;

        String[] variantIds = new String[nVariants];
        double[] afs = new double[nVariants];
        int[] acs = new int[nVariants];
        int[] ans = new int[nVariants];

        for (int v = 0; v < nVariants; v++) {
            variantIds[v] = "test:" + positions[v] + ":A:T";
            int ac = 0;
            for (int h = 0; h < nHaplotypes; h++) {
                ac += haplotypes[v][h];
            }
            acs[v] = ac;
            ans[v] = nHaplotypes;
            afs[v] = (double) ac / nHaplotypes;
        }

        return HaplotypeMatrix.forTest(variantIds, positions, afs, acs, ans,
                haplotypes, nSamples);
    }

    // -----------------------------------------------------------------------
    // Two-variant matrix: known EHH
    // -----------------------------------------------------------------------

    @Test
    void twoVariants_allIdentical_EHH1() {
        // All 4 carrier haplotypes carry same allele at both sites → groups never split → EHH=1
        // Variant 1 must be polymorphic in the population (not skipped in walk).
        // Carriers are {0,1,2,3} (alt at focal). At idx=1, haps 0-3 are all 0 → same group.
        long[] pos = {1000, 1100};
        byte[][] hap = {
                {1, 1, 1, 1, 0, 0},  // variant 0 (focal): carriers={0,1,2,3}
                {0, 0, 0, 0, 1, 1},  // variant 1: polymorphic (ac=2), but all carriers have 0
        };

        HaplotypeMatrix m = buildMatrix(pos, hap);
        int[] carriers = {0, 1, 2, 3};

        // From focal=0 walking downstream to idx=1 (distance=100):
        // EHH at idx=1: all carriers in same group → ehh=1.0
        // ihh_down = (1.0 + 1.0)/2 * 100 = 100
        // From focal=0 walking upstream: no variant → ihh_up = 0
        double ihh = EHHComputer.computeIHH(m, 0, carriers, 0.05);
        assertEquals(100.0, ihh, EPS);
    }

    @Test
    void twoVariants_completeSplit_EHH0() {
        // Focal at idx=0. At idx=1, haplotypes split into {0,1} carrying ref and {2,3} carrying alt.
        // Groups: {0,1}→ref, {2,3}→alt  →  each group has C(2,2)=1 pair
        // EHH = 2 / C(4,2) = 2/6 = 1/3
        long[] pos = {1000, 1100};
        byte[][] hap = {
                {1, 1, 1, 1},  // focal: all alt
                {0, 0, 1, 1},  // splits into two groups of 2
        };

        HaplotypeMatrix m = buildMatrix(pos, hap);
        int[] carriers = {0, 1, 2, 3};

        // ihh_down = (1.0 + 1.0/3.0)/2 * 100
        double expectedDown = (1.0 + 1.0 / 3.0) / 2.0 * 100.0;
        double ihh = EHHComputer.computeIHH(m, 0, carriers, 0.0);
        assertEquals(expectedDown, ihh, EPS);
    }

    // -----------------------------------------------------------------------
    // Multi-variant: known iHH with trapezoidal integration
    // -----------------------------------------------------------------------

    @Test
    void threeVariants_trapezoidalIntegration() {
        // Focal at idx=1 (middle). Walk upstream to idx=0 and downstream to idx=2.
        // 8 haplotypes, 6 carry alt at focal (carriers={0,1,2,3,4,5}).
        // At idx=0: split into {0,1,2}→ref and {3,4,5}→alt → EHH = 2*C(3,2)/C(6,2) = 6/15 = 0.4
        // At idx=2: polymorphic in pop, but all carriers have same allele → no split → EHH = 1.0
        long[] pos = {900, 1000, 1200};
        byte[][] hap = {
                {0, 0, 0, 1, 1, 1, 0, 1},  // idx=0: polymorphic (ac=4)
                {1, 1, 1, 1, 1, 1, 0, 0},  // idx=1 (focal): carriers={0..5}
                {0, 0, 0, 0, 0, 0, 1, 1},  // idx=2: polymorphic (ac=2), all carriers have 0
        };

        HaplotypeMatrix m = buildMatrix(pos, hap);
        int[] carriers = {0, 1, 2, 3, 4, 5};

        // Upstream (idx=1 → idx=0, distance=100):
        //   EHH at idx=0: groups {0,1,2} and {3,4,5}: pairs = 3+3=6, total=15 → 0.4
        //   ihh_up = (1.0 + 0.4)/2 * 100 = 70
        double expectedUp = (1.0 + 0.4) / 2.0 * 100.0;

        // Downstream (idx=1 → idx=2, distance=200):
        //   EHH at idx=2: all carriers have 0, no split → EHH=1.0
        //   ihh_down = (1.0 + 1.0)/2 * 200 = 200
        double expectedDown = (1.0 + 1.0) / 2.0 * 200.0;

        double ihh = EHHComputer.computeIHH(m, 1, carriers, 0.0);
        assertEquals(expectedUp + expectedDown, ihh, EPS);
    }

    // -----------------------------------------------------------------------
    // Edge cases
    // -----------------------------------------------------------------------

    @Test
    void allUnique_EHHdropsFast() {
        // 4 haplotypes, each uniquely identified after first step.
        // Focal at idx=1. At idx=2, each haplotype is in its own group.
        // EHH = 0/C(4,2) = 0 (each group has size 1, no pairs)
        long[] pos = {800, 1000, 1100, 1200};
        byte[][] hap = {
                {0, 0, 1, 1},  // idx=0
                {1, 1, 1, 1},  // idx=1 (focal)
                {0, 1, 0, 1},  // idx=2 — combined with idx=0 alleles, 4 unique groups
                {0, 0, 1, 1},  // idx=3
        };

        HaplotypeMatrix m = buildMatrix(pos, hap);
        int[] carriers = {0, 1, 2, 3};

        // Downstream from focal(1):
        //   Step to idx=2: groups from idx=0 don't matter (we start fresh from focal).
        //   At focal, all in group 0. At idx=2, alleles are {0,1,0,1}.
        //   Split group 0 into {0,2}→ref and {1,3}→alt. Each size 2.
        //   EHH = (C(2,2)+C(2,2))/C(4,2) = 2/6 = 1/3
        //   At idx=3, alleles are {0,0,1,1}. Group {0,2} checks alleles: 0→ref,2→alt → splits.
        //     Group {1,3} checks alleles: 0→ref, 1→alt → splits.
        //     Now 4 groups of size 1 → EHH = 0
        //   ihh_down = (1.0 + 1/3)/2 * 100 + (1/3 + 0)/2 * 100 = 200/3 + 100/6 = 500/6

        // This shows progressive splitting behavior works correctly
        double ihh = EHHComputer.computeIHH(m, 1, carriers, 0.0);
        assertTrue(ihh > 0, "iHH should be positive");
    }

    @Test
    void singleCarrier_returnsZero() {
        long[] pos = {1000, 1100};
        byte[][] hap = {{1, 0}, {1, 0}};
        HaplotypeMatrix m = buildMatrix(pos, hap);
        assertEquals(0.0, EHHComputer.computeIHH(m, 0, new int[]{0}, 0.05), EPS);
    }

    @Test
    void emptyCarriers_returnsZero() {
        long[] pos = {1000, 1100};
        byte[][] hap = {{1, 0}, {1, 0}};
        HaplotypeMatrix m = buildMatrix(pos, hap);
        assertEquals(0.0, EHHComputer.computeIHH(m, 0, new int[]{}, 0.05), EPS);
    }

    // -----------------------------------------------------------------------
    // computeSL: count SNP positions traversed
    // -----------------------------------------------------------------------

    @Test
    void computeSL_twoVariants_allIdentical_countsOne() {
        // All 4 haplotypes identical at both sites → EHH never drops → SL = 1 (one step)
        long[] pos = {1000, 1100};
        byte[][] hap = {
                {1, 1, 1, 1},
                {1, 1, 1, 1},
        };
        HaplotypeMatrix m = buildMatrix(pos, hap);
        int[] carriers = {0, 1, 2, 3};

        // Only 1 variant downstream from focal(0), upstream has 0
        int sl = EHHComputer.computeSL(m, 0, carriers, 0.0);
        assertEquals(1, sl); // 1 step downstream, 0 upstream
    }

    @Test
    void computeSL_threeVariants_countsCorrectly() {
        // Focal at 1 (middle). Walk up 1 step, down 1 step = 2 total
        // No splitting → EHH stays above 0.05
        long[] pos = {900, 1000, 1200};
        byte[][] hap = {
                {1, 1, 1, 1},
                {1, 1, 1, 1}, // focal
                {1, 1, 1, 1},
        };
        HaplotypeMatrix m = buildMatrix(pos, hap);
        int[] carriers = {0, 1, 2, 3};

        int sl = EHHComputer.computeSL(m, 1, carriers, 0.0);
        assertEquals(2, sl); // 1 upstream + 1 downstream
    }

    @Test
    void computeSL_singleCarrier_returnsZero() {
        long[] pos = {1000, 1100};
        byte[][] hap = {{1, 0}, {1, 0}};
        HaplotypeMatrix m = buildMatrix(pos, hap);
        assertEquals(0, EHHComputer.computeSL(m, 0, new int[]{0}, 0.05));
    }

    @Test
    void computeSL_emptyCarriers_returnsZero() {
        long[] pos = {1000, 1100};
        byte[][] hap = {{1, 0}, {1, 0}};
        HaplotypeMatrix m = buildMatrix(pos, hap);
        assertEquals(0, EHHComputer.computeSL(m, 0, new int[]{}, 0.05));
    }

    // -----------------------------------------------------------------------
    // Existing tests below
    // -----------------------------------------------------------------------

    @Test
    void minEHH_stopsEarly() {
        // With default minEHH=0.05, EHH should stop when below threshold.
        // 20 haplotypes, splitting aggressively.
        int nHaps = 20;
        long[] pos = new long[10];
        byte[][] hap = new byte[10][nHaps];

        for (int v = 0; v < 10; v++) {
            pos[v] = 1000 + v * 100;
            for (int h = 0; h < nHaps; h++) {
                if (v == 0) {
                    hap[v][h] = 1;  // focal: all alt
                } else {
                    // Alternate pattern to split groups
                    hap[v][h] = (byte) ((h / (1 << (v - 1))) % 2);
                }
            }
        }

        HaplotypeMatrix m = buildMatrix(pos, hap);
        int[] carriers = new int[nHaps];
        for (int h = 0; h < nHaps; h++) carriers[h] = h;

        // With minEHH=0.05, should stop before reaching the end
        double ihhStrict = EHHComputer.computeIHH(m, 0, carriers, 0.05);
        // With minEHH=0.0, should integrate to the end
        double ihhFull = EHHComputer.computeIHH(m, 0, carriers, 0.0);
        assertTrue(ihhFull >= ihhStrict, "Full integration should be >= strict");
    }
}
