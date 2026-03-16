package org.graphpop.procedures;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Unit tests for ROH detection logic.
 *
 * <p>Tests the core sliding-window homozygosity scanning algorithm
 * using synthetic haplotype matrices. ROHProcedure itself requires
 * a Neo4j transaction, so we test the underlying logic indirectly.</p>
 */
class ROHScannerTest {

    /**
     * Build a HaplotypeMatrix for a single sample from two haplotype arrays.
     */
    private static HaplotypeMatrix buildSingleSample(long[] positions, int[] hap0, int[] hap1) {
        int nVariants = positions.length;
        int nSamples = 1;
        byte[][] haplotypes = new byte[nVariants][2];
        String[] variantIds = new String[nVariants];
        int[] acs = new int[nVariants];
        int[] ans = new int[nVariants];
        double[] afs = new double[nVariants];

        for (int v = 0; v < nVariants; v++) {
            haplotypes[v][0] = (byte) hap0[v];
            haplotypes[v][1] = (byte) hap1[v];
            variantIds[v] = "test:" + positions[v] + ":A:T";
            acs[v] = hap0[v] + hap1[v];
            ans[v] = 2;
            afs[v] = acs[v] / 2.0;
        }

        return HaplotypeMatrix.forTest(variantIds, positions, afs, acs, ans, haplotypes, nSamples);
    }

    /**
     * Simulate the ROH sliding window scan for a single sample.
     */
    private static int[] scanROH(HaplotypeMatrix matrix, int sampleIdx, int windowSnps, int maxHet) {
        int nVariants = matrix.nVariants;
        boolean[] isHet = new boolean[nVariants];
        for (int v = 0; v < nVariants; v++) {
            isHet[v] = matrix.hap(v, 2 * sampleIdx) != matrix.hap(v, 2 * sampleIdx + 1);
        }

        boolean[] homWindow = new boolean[nVariants];
        int hetCount = 0;
        for (int v = 0; v < windowSnps && v < nVariants; v++) {
            if (isHet[v]) hetCount++;
        }
        if (hetCount <= maxHet) {
            for (int v = 0; v < windowSnps && v < nVariants; v++) homWindow[v] = true;
        }
        for (int v = 1; v <= nVariants - windowSnps; v++) {
            if (isHet[v - 1]) hetCount--;
            if (isHet[v + windowSnps - 1]) hetCount++;
            if (hetCount <= maxHet) {
                for (int w = v; w < v + windowSnps; w++) homWindow[w] = true;
            }
        }

        // Count contiguous segments
        int nSegs = 0;
        boolean inSeg = false;
        for (boolean b : homWindow) {
            if (b && !inSeg) { nSegs++; inSeg = true; }
            else if (!b) inSeg = false;
        }

        int trueCount = 0;
        for (boolean b : homWindow) if (b) trueCount++;
        return new int[]{nSegs, trueCount};
    }

    @Test
    void allHomozygous_oneSegment() {
        int n = 100;
        long[] pos = new long[n];
        int[] hap0 = new int[n];
        int[] hap1 = new int[n];
        for (int i = 0; i < n; i++) {
            pos[i] = 1000000 + i * 1000;
            hap0[i] = 0; hap1[i] = 0;  // all homozygous ref
        }

        HaplotypeMatrix m = buildSingleSample(pos, hap0, hap1);
        int[] result = scanROH(m, 0, 10, 0);
        assertEquals(1, result[0], "All homozygous should give one segment");
    }

    @Test
    void allHeterozygous_noSegment() {
        int n = 100;
        long[] pos = new long[n];
        int[] hap0 = new int[n];
        int[] hap1 = new int[n];
        for (int i = 0; i < n; i++) {
            pos[i] = 1000000 + i * 1000;
            hap0[i] = 0; hap1[i] = 1;  // all heterozygous
        }

        HaplotypeMatrix m = buildSingleSample(pos, hap0, hap1);
        int[] result = scanROH(m, 0, 10, 0);
        assertEquals(0, result[0], "All heterozygous should give no segments");
    }

    @Test
    void mixedPattern_twoSegments() {
        // Pattern: 40 hom + 20 het + 40 hom
        int n = 100;
        long[] pos = new long[n];
        int[] hap0 = new int[n];
        int[] hap1 = new int[n];
        for (int i = 0; i < n; i++) {
            pos[i] = 1000000 + i * 1000;
            if (i >= 40 && i < 60) {
                hap0[i] = 0; hap1[i] = 1;  // het
            } else {
                hap0[i] = 0; hap1[i] = 0;  // hom
            }
        }

        HaplotypeMatrix m = buildSingleSample(pos, hap0, hap1);
        int[] result = scanROH(m, 0, 10, 0);
        assertEquals(2, result[0], "Should detect two homozygous segments");
    }

    @Test
    void maxHet_allowsSomeHeterozygosity() {
        // Mostly homozygous with a few scattered hets
        int n = 100;
        long[] pos = new long[n];
        int[] hap0 = new int[n];
        int[] hap1 = new int[n];
        for (int i = 0; i < n; i++) {
            pos[i] = 1000000 + i * 1000;
            if (i == 25 || i == 75) {
                hap0[i] = 0; hap1[i] = 1;  // sparse hets
            } else {
                hap0[i] = 0; hap1[i] = 0;
            }
        }

        HaplotypeMatrix m = buildSingleSample(pos, hap0, hap1);

        // With maxHet=0, the het sites should break the segments
        int[] strict = scanROH(m, 0, 10, 0);

        // With maxHet=1, windows containing 1 het are still considered homozygous
        int[] lenient = scanROH(m, 0, 10, 1);

        assertTrue(lenient[1] >= strict[1],
                "Allowing 1 het per window should give same or more homozygous positions");
    }

    @Test
    void shortRegion_lessThanWindowSize() {
        // Only 5 variants but window requires 10
        // The initial window check runs 0..min(10,5)=5 and marks them if
        // het count is low enough. All hom → they get marked. This is correct
        // behavior — ROHProcedure checks nVariants < windowSnps before scanning,
        // but the raw scan algorithm handles short regions by scanning what's available.
        long[] pos = {1000, 2000, 3000, 4000, 5000};
        int[] hap0 = {0, 0, 0, 0, 0};
        int[] hap1 = {0, 0, 0, 0, 0};

        HaplotypeMatrix m = buildSingleSample(pos, hap0, hap1);
        int[] result = scanROH(m, 0, 10, 0);
        // With 5 hom variants and window_snps=10, initial window still scans
        // and marks all as homozygous, creating 1 segment
        assertEquals(1, result[0], "Short all-hom region still detected as one segment");
    }
}
