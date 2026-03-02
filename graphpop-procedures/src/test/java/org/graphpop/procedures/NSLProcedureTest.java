package org.graphpop.procedures;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Unit tests for the pairwise SSL nSL computation in {@link NSLProcedure}.
 *
 * <p>Tests verify that scanDirection() produces correct mean SSL values
 * matching the Ferrer-Admetlla et al. (2014) algorithm.</p>
 */
class NSLProcedureTest {

    private static final double EPS = 1e-10;

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
    // Forward scan: basic SSL accumulation
    // -----------------------------------------------------------------------

    @Test
    void twoIdenticalVariants_forwardScan() {
        // 4 haplotypes, 2 focal variants. All haplotypes identical across both sites.
        // At variant 0: pairs (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)
        //   All carry allele 0 at site 0 → ssl becomes 1 for all pairs
        //   meanSL0 = 6*1/6 = 1.0, meanSL1 = 0 (no allele-1 pairs)
        // At variant 1: ssl becomes 2 for all pairs, meanSL0 = 6*2/6 = 2.0
        long[] pos = {1000, 1100};
        byte[][] hap = {
                {0, 0, 0, 0},
                {0, 0, 0, 0},
        };
        HaplotypeMatrix m = buildMatrix(pos, hap);
        int[] focalIndices = {0, 1};
        double[] sl0 = new double[2];
        double[] sl1 = new double[2];

        NSLProcedure.scanDirection(m, focalIndices, 4, sl0, sl1, +1);

        assertEquals(1.0, sl0[0], EPS); // first variant: all 6 pairs, ssl=1
        assertEquals(0.0, sl1[0], EPS); // no allele-1 pairs
        assertEquals(2.0, sl0[1], EPS); // second variant: all 6 pairs, ssl=2
        assertEquals(0.0, sl1[1], EPS);
    }

    @Test
    void mismatchResetsSsl() {
        // 4 haplotypes, 3 focal variants.
        // Site 0: all carry 0 → ssl=1 for all pairs
        // Site 1: haps split (0,1)→allele0, (2,3)→allele1
        //   pairs (0,1): same allele 0 at both sites → ssl = 2, contributes to sl0
        //   pairs (2,3): same allele 1 at both sites → ssl was 1 at site0 (but they had allele0),
        //     now they have allele1 → mismatch from previous (same allele, but ssl state carries)
        //   Actually: at site0, all have allele 0, ssl=1.
        //   At site1: hap0=0, hap1=0, hap2=1, hap3=1
        //   Pair (0,1): both 0 → ssl=1+1=2, sl0 += 2
        //   Pair (0,2): 0 vs 1 → mismatch → ssl=0
        //   Pair (0,3): 0 vs 1 → mismatch → ssl=0
        //   Pair (1,2): 0 vs 1 → mismatch → ssl=0
        //   Pair (1,3): 0 vs 1 → mismatch → ssl=0
        //   Pair (2,3): both 1 → ssl=1+1=2, sl1 += 2
        //   meanSL0 = 2/1 = 2.0, meanSL1 = 2/1 = 2.0

        long[] pos = {1000, 1100, 1200};
        byte[][] hap = {
                {0, 0, 0, 0},
                {0, 0, 1, 1},
                {0, 0, 1, 1},
        };
        HaplotypeMatrix m = buildMatrix(pos, hap);
        int[] focalIndices = {0, 1, 2};
        double[] sl0 = new double[3];
        double[] sl1 = new double[3];

        NSLProcedure.scanDirection(m, focalIndices, 4, sl0, sl1, +1);

        // Site 0: all pairs have allele 0, ssl=1 each, 6 pairs → mean=1.0
        assertEquals(1.0, sl0[0], EPS);
        assertEquals(0.0, sl1[0], EPS);

        // Site 1: pair(0,1) ssl=2, others reset. pair(2,3) ssl=2
        assertEquals(2.0, sl0[1], EPS); // 1 pair (0,1) with ssl=2 → mean=2/1=2
        assertEquals(2.0, sl1[1], EPS); // 1 pair (2,3) with ssl=2 → mean=2/1=2

        // Site 2: same alleles as site 1.
        // pair(0,1): still both 0 → ssl=3, sl0 = 3/1 = 3
        // pair(2,3): still both 1 → ssl=3, sl1 = 3/1 = 3
        // Mismatched pairs: (0,2),(0,3),(1,2),(1,3) still mismatch → ssl stays 0
        assertEquals(3.0, sl0[2], EPS);
        assertEquals(3.0, sl1[2], EPS);
    }

    @Test
    void backwardScan_reverseOrder() {
        // Same data as mismatchResetsSsl but scanning backward.
        // Backward visits sites in order: 2, 1, 0
        long[] pos = {1000, 1100, 1200};
        byte[][] hap = {
                {0, 0, 0, 0},
                {0, 0, 1, 1},
                {0, 0, 1, 1},
        };
        HaplotypeMatrix m = buildMatrix(pos, hap);
        int[] focalIndices = {0, 1, 2};
        double[] sl0 = new double[3];
        double[] sl1 = new double[3];

        NSLProcedure.scanDirection(m, focalIndices, 4, sl0, sl1, -1);

        // Backward: first visit site 2, then 1, then 0
        // Site 2 (first visited): ssl=1 for matching pairs
        //   hap = {0, 0, 1, 1}
        //   pair(0,1): both 0, ssl=1 → sl0 sum=1, n00=1
        //   pair(0,2): 0 vs 1 → mismatch
        //   pair(2,3): both 1, ssl=1 → sl1 sum=1, n11=1
        assertEquals(1.0, sl0[2], EPS);
        assertEquals(1.0, sl1[2], EPS);

        // Site 1: same alleles as site 2 → ssl accumulates
        //   pair(0,1): both 0, ssl=2 → sl0=2/1
        //   pair(2,3): both 1, ssl=2 → sl1=2/1
        assertEquals(2.0, sl0[1], EPS);
        assertEquals(2.0, sl1[1], EPS);

        // Site 0: all allele 0 → pair(0,1): ssl=3, pair(2,3) had allele 1 now allele 0: mismatch, ssl=0
        //   pair(0,1): both 0, ssl=3 → contributes 3
        //   pair(0,2): 0 vs 0 → same! ssl was 0 (mismatch at site1) → ssl=1
        //   pair(0,3): same → ssl=1
        //   pair(1,2): same → ssl=1
        //   pair(1,3): same → ssl=1
        //   pair(2,3): 0 vs 0 → they had allele 1 at sites 1,2 but now allele 0: ssl was 2, but
        //     wait — at sites 2,1 pair(2,3) both had allele 1, so ssl=2. Now at site 0, both have allele 0 → same → ssl=3
        //   So all pairs now match: (0,1)→ssl=3, (0,2)→ssl=1, (0,3)→ssl=1, (1,2)→ssl=1, (1,3)→ssl=1, (2,3)→ssl=3
        //   All are allele 0: n00=6, sum=3+1+1+1+1+3=10, mean=10/6
        assertEquals(10.0 / 6.0, sl0[0], EPS);
        assertEquals(0.0, sl1[0], EPS); // no allele-1 at site 0
    }

    @Test
    void nslCombinedForwardBackward() {
        // 6 haplotypes, 5 focal variants, simple pattern
        // This tests the full nSL computation for a middle variant
        long[] pos = {1000, 1100, 1200, 1300, 1400};
        byte[][] hap = {
                {0, 0, 0, 1, 1, 1},  // site 0: 3 ref, 3 alt
                {0, 0, 0, 1, 1, 1},  // site 1: same
                {0, 0, 0, 1, 1, 1},  // site 2: same (focal for test)
                {0, 0, 0, 1, 1, 1},  // site 3: same
                {0, 0, 0, 1, 1, 1},  // site 4: same
        };
        HaplotypeMatrix m = buildMatrix(pos, hap);
        int[] focalIndices = {0, 1, 2, 3, 4};
        int nHaps = 6;

        double[] sl0Fwd = new double[5];
        double[] sl1Fwd = new double[5];
        double[] sl0Rev = new double[5];
        double[] sl1Rev = new double[5];

        NSLProcedure.scanDirection(m, focalIndices, nHaps, sl0Fwd, sl1Fwd, +1);
        NSLProcedure.scanDirection(m, focalIndices, nHaps, sl0Rev, sl1Rev, -1);

        // For the middle variant (index 2):
        // Forward: sites 0, 1, 2 → ssl accumulates to 3 for matching pairs
        // All allele-0 pairs: (0,1),(0,2),(1,2) = 3 pairs, each ssl=3 → mean=3
        // All allele-1 pairs: (3,4),(3,5),(4,5) = 3 pairs, each ssl=3 → mean=3
        assertEquals(3.0, sl0Fwd[2], EPS);
        assertEquals(3.0, sl1Fwd[2], EPS);

        // Backward visits: 4, 3, 2 → ssl accumulates to 3 for matching pairs at site 2
        assertEquals(3.0, sl0Rev[2], EPS);
        assertEquals(3.0, sl1Rev[2], EPS);

        // Combined: SL0 = 3+3 = 6, SL1 = 3+3 = 6 → nSL = log(6/6) = 0
        double sl0 = sl0Fwd[2] + sl0Rev[2];
        double sl1 = sl1Fwd[2] + sl1Rev[2];
        double nsl = Math.log(sl1 / sl0);
        assertEquals(0.0, nsl, EPS, "Symmetric haplotypes should give nSL=0");
    }

    @Test
    void nslAsymmetric_positiveScore() {
        // Allele 1 carriers are more homogeneous → SL1 > SL0 → nSL > 0
        // 4 haplotypes: hap0,hap1 carry allele 0; hap2,hap3 carry allele 1
        // At flanking sites, allele-0 carriers diverge but allele-1 carriers stay identical
        long[] pos = {1000, 1100, 1200};
        byte[][] hap = {
                {0, 1, 1, 1},  // site 0: hap0 differs from hap1 among allele-0 carriers
                {0, 0, 1, 1},  // site 1: focal — hap0,1 carry 0; hap2,3 carry 1
                {0, 1, 1, 1},  // site 2: hap0 differs from hap1 again
        };
        HaplotypeMatrix m = buildMatrix(pos, hap);
        int[] focalIndices = {0, 1, 2};
        int nHaps = 4;

        double[] sl0Fwd = new double[3];
        double[] sl1Fwd = new double[3];
        double[] sl0Rev = new double[3];
        double[] sl1Rev = new double[3];

        NSLProcedure.scanDirection(m, focalIndices, nHaps, sl0Fwd, sl1Fwd, +1);
        NSLProcedure.scanDirection(m, focalIndices, nHaps, sl0Rev, sl1Rev, -1);

        // Forward scan visits: site 0, site 1, site 2
        // At site 1 (focal, fi=1):
        //   pair(0,1): alleles 0,0 → same → ssl was: at site0 pair(0,1) had alleles 0,1 → mismatch → ssl=0
        //     so ssl=0+1=1 → contributes 1 to sl0
        //   pair(2,3): alleles 1,1 → same → ssl was: at site0 pair(2,3) had alleles 1,1 → ssl=1
        //     so ssl=1+1=2 → contributes 2 to sl1
        //   pair(0,2): 0 vs 1 → mismatch
        //   pair(0,3): 0 vs 1 → mismatch
        //   pair(1,2): 0 vs 1 → mismatch
        //   pair(1,3): 0 vs 1 → mismatch
        // sl0Fwd[1] = 1/1 = 1.0, sl1Fwd[1] = 2/1 = 2.0
        assertEquals(1.0, sl0Fwd[1], EPS);
        assertEquals(2.0, sl1Fwd[1], EPS);

        // Backward scan visits: site 2, site 1, site 0
        // At site 1 (fi=1):
        //   Backward from site 2: pair(0,1) at site 2 had alleles 0,1 → mismatch → ssl=0
        //     at site 1: both 0 → ssl=1
        //   pair(2,3) at site 2 had alleles 1,1 → ssl=1; at site 1: both 1 → ssl=2
        // sl0Rev[1] = 1/1 = 1.0, sl1Rev[1] = 2/1 = 2.0
        assertEquals(1.0, sl0Rev[1], EPS);
        assertEquals(2.0, sl1Rev[1], EPS);

        // nSL for focal site 1: log((2+2)/(1+1)) = log(4/2) = log(2) > 0
        double sl0 = sl0Fwd[1] + sl0Rev[1];
        double sl1 = sl1Fwd[1] + sl1Rev[1];
        double nsl = Math.log(sl1 / sl0);
        assertEquals(Math.log(2.0), nsl, EPS, "Allele 1 more homogeneous → positive nSL");
    }

    @Test
    void nslAsymmetric_negativeScore() {
        // Allele 0 carriers more homogeneous → SL0 > SL1 → nSL < 0
        long[] pos = {1000, 1100, 1200};
        byte[][] hap = {
                {0, 0, 0, 1},  // allele-1 carriers (hap2,hap3) diverge at flanking sites
                {0, 0, 1, 1},  // focal: hap0,1 carry 0; hap2,3 carry 1
                {0, 0, 0, 1},  // allele-1 carriers diverge again
        };
        HaplotypeMatrix m = buildMatrix(pos, hap);
        int[] focalIndices = {0, 1, 2};
        int nHaps = 4;

        double[] sl0Fwd = new double[3];
        double[] sl1Fwd = new double[3];
        double[] sl0Rev = new double[3];
        double[] sl1Rev = new double[3];

        NSLProcedure.scanDirection(m, focalIndices, nHaps, sl0Fwd, sl1Fwd, +1);
        NSLProcedure.scanDirection(m, focalIndices, nHaps, sl0Rev, sl1Rev, -1);

        // At focal (site 1), forward:
        //   pair(0,1): both allele 0 everywhere → ssl accumulates from site 0 → ssl=2
        //   pair(2,3): allele 1,1 at site 1 but at site 0 they had 0,1 → mismatch → ssl=0 then ssl=1
        // sl0Fwd[1] = 2.0, sl1Fwd[1] = 1.0
        // Similarly backward
        double sl0 = sl0Fwd[1] + sl0Rev[1];
        double sl1 = sl1Fwd[1] + sl1Rev[1];
        double nsl = Math.log(sl1 / sl0);
        assertTrue(nsl < 0, "Allele 0 more homogeneous → negative nSL, got " + nsl);
        assertEquals(Math.log(0.5), nsl, EPS);
    }

    @Test
    void singlePairPerAllele_minimalCase() {
        // Minimal case: 4 haplotypes, 2 per allele, single pair each
        // 2 focal variants
        long[] pos = {1000, 1100};
        byte[][] hap = {
                {0, 0, 1, 1},
                {0, 0, 1, 1},
        };
        HaplotypeMatrix m = buildMatrix(pos, hap);
        int[] focalIndices = {0, 1};
        int nHaps = 4;

        double[] sl0Fwd = new double[2];
        double[] sl1Fwd = new double[2];

        NSLProcedure.scanDirection(m, focalIndices, nHaps, sl0Fwd, sl1Fwd, +1);

        // Site 0: pair(0,1) both 0 → ssl=1, pair(2,3) both 1 → ssl=1
        assertEquals(1.0, sl0Fwd[0], EPS);
        assertEquals(1.0, sl1Fwd[0], EPS);

        // Site 1: pair(0,1) ssl=2, pair(2,3) ssl=2
        assertEquals(2.0, sl0Fwd[1], EPS);
        assertEquals(2.0, sl1Fwd[1], EPS);
    }
}
