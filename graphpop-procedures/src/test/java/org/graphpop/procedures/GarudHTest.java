package org.graphpop.procedures;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Unit tests for Garud et al. (2015) haplotype homozygosity statistics.
 */
class GarudHTest {

    private static final double EPS = 1e-10;

    /** Pack unpacked byte[][] (1 byte/haplotype) to bit-packed format. */
    private static byte[][] pack(byte[][] unpacked) {
        byte[][] packed = new byte[unpacked.length][];
        for (int v = 0; v < unpacked.length; v++) {
            packed[v] = HaplotypeMatrix.packRow(unpacked[v]);
        }
        return packed;
    }

    @Test
    void allIdenticalHaplotypes() {
        // 10 haplotypes, 5 variants, all ref → 1 unique haplotype
        int nHap = 10;
        int nVar = 5;
        byte[][] hap = new byte[nVar][nHap]; // all zeros

        GarudH.Result r = GarudH.compute(pack(hap), 0, nVar - 1, nHap);

        assertNotNull(r);
        assertEquals(1.0, r.h1, EPS, "All identical → H1 = 1.0");
        assertEquals(1.0, r.h12, EPS, "All identical → H12 = 1.0");
        assertEquals(0.0, r.h2_h1, EPS, "All identical → H2/H1 = 0");
        assertEquals(1, r.n_haplotypes);
        assertEquals(0.0, r.hap_diversity, EPS, "All identical → Hd = 0");
    }

    @Test
    void allUniqueHaplotypes() {
        // 8 haplotypes, 3 variants, each haplotype different
        int nHap = 8;
        int nVar = 3;
        byte[][] hap = new byte[nVar][nHap];
        // Assign unique 3-bit patterns: 000, 001, 010, 011, 100, 101, 110, 111
        for (int h = 0; h < nHap; h++) {
            for (int v = 0; v < nVar; v++) {
                hap[v][h] = (byte) ((h >> (nVar - 1 - v)) & 1);
            }
        }

        GarudH.Result r = GarudH.compute(pack(hap), 0, nVar - 1, nHap);

        assertNotNull(r);
        assertEquals(8, r.n_haplotypes);
        // H1 = 8 × (1/8)² = 8/64 = 1/8
        assertEquals(1.0 / 8.0, r.h1, EPS);
        // H12 = (1/8 + 1/8)² + 6 × (1/8)² = 4/64 + 6/64 = 10/64 = 5/32
        assertEquals(5.0 / 32.0, r.h12, EPS);
        // H2/H1 = (1/8 - 1/64)/(1/8) = (7/64)/(1/8) = 7/8
        assertEquals(7.0 / 8.0, r.h2_h1, EPS);
        // Hd = (8/7) × (1 - 1/8) = (8/7) × (7/8) = 1.0
        assertEquals(1.0, r.hap_diversity, EPS);
    }

    @Test
    void hardSweepPattern() {
        // 1 dominant haplotype (8 copies) + 2 rare (1 copy each) = 10 haplotypes
        // Hard sweep: one haplotype dominates
        int nHap = 10;
        int nVar = 4;
        byte[][] hap = new byte[nVar][nHap];
        // First 8 haplotypes: all 0000 (dominant)
        // Haplotype 8: 1000
        hap[0][8] = 1;
        // Haplotype 9: 0100
        hap[1][9] = 1;

        GarudH.Result r = GarudH.compute(pack(hap), 0, nVar - 1, nHap);

        assertNotNull(r);
        assertEquals(3, r.n_haplotypes);
        // f1=0.8, f2=0.1, f3=0.1
        // H1 = 0.64 + 0.01 + 0.01 = 0.66
        assertEquals(0.66, r.h1, EPS);
        // H12 = (0.8+0.1)² + 0.1² = 0.81 + 0.01 = 0.82
        assertEquals(0.82, r.h12, EPS);
        // H2/H1 = (0.66 - 0.64)/0.66 = 0.02/0.66 ≈ 0.0303
        assertEquals(0.02 / 0.66, r.h2_h1, EPS);
    }

    @Test
    void softSweepPattern() {
        // 2 common haplotypes (4 copies each) + 2 rare (1 each) = 10 haplotypes
        int nHap = 10;
        int nVar = 3;
        byte[][] hap = new byte[nVar][nHap];
        // First 4: pattern 000
        // Next 4: pattern 111
        for (int h = 4; h < 8; h++) {
            for (int v = 0; v < nVar; v++) hap[v][h] = 1;
        }
        // Haplotype 8: 100
        hap[0][8] = 1;
        // Haplotype 9: 010
        hap[1][9] = 1;

        GarudH.Result r = GarudH.compute(pack(hap), 0, nVar - 1, nHap);

        assertNotNull(r);
        assertEquals(4, r.n_haplotypes);
        // f1=0.4, f2=0.4, f3=0.1, f4=0.1
        // H1 = 0.16 + 0.16 + 0.01 + 0.01 = 0.34
        assertEquals(0.34, r.h1, EPS);
        // H12 = (0.4+0.4)² + 0.01 + 0.01 = 0.64 + 0.02 = 0.66
        assertEquals(0.66, r.h12, EPS);
        // H2/H1 = (0.34 - 0.16)/0.34 = 0.18/0.34 ≈ 0.529
        double expectedH2H1 = 0.18 / 0.34;
        assertEquals(expectedH2H1, r.h2_h1, EPS);
        // Soft sweep: H2/H1 much higher than hard sweep
        assertTrue(r.h2_h1 > 0.4, "Soft sweep should have elevated H2/H1");
    }

    @Test
    void subrangeComputation() {
        // Verify that compute works with a sub-range of variants
        int nHap = 6;
        int nVar = 6;
        byte[][] hap = new byte[nVar][nHap];
        // Make all unique within vars 2-4 (sub-range)
        for (int h = 0; h < nHap; h++) {
            hap[2][h] = (byte) ((h >> 2) & 1);
            hap[3][h] = (byte) ((h >> 1) & 1);
            hap[4][h] = (byte) (h & 1);
        }

        // Full range vs sub-range should give different results
        GarudH.Result rFull = GarudH.compute(pack(hap), 0, nVar - 1, nHap);
        GarudH.Result rSub = GarudH.compute(pack(hap), 2, 4, nHap);

        assertNotNull(rFull);
        assertNotNull(rSub);
        // Sub-range has 6 unique haplotypes (3 bits, 6 < 8 possible)
        assertEquals(6, rSub.n_haplotypes);
        // Each haplotype frequency = 1/6
        assertEquals(6.0 * (1.0 / 36.0), rSub.h1, EPS);
    }
}
