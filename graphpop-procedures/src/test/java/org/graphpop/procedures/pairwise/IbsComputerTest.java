package org.graphpop.procedures.pairwise;

import org.graphpop.procedures.PackedGenotypeReader;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Hand-computed unit tests for {@link IbsComputer} on the same
 * 4-sample × 8-variant fixture as {@link KingRobustComputerTest}.
 *
 * <pre>
 *      V1 V2 V3 V4 V5 V6 V7 V8
 * S1:  1  0  2  1  0  1  2  1
 * S2:  1  1  1  1  1  0  2  0
 * S3:  2  2  0  0  2  0  1  2
 * S4:  0  0  0  0  0  0  0  0
 * </pre>
 *
 * <p>Hand-computed pair counters for each pair (ibs0, ibs1, ibs2,
 * n_used) and IBS = (2*ibs2 + ibs1) / (2 * n_used):</p>
 * <ul>
 *   <li>(S1, S2): diffs={0,1,1,0,1,1,0,1} -&gt; ibs2=3, ibs1=5, ibs0=0
 *       -&gt; IBS = (6+5)/16 = 11/16 = 0.6875</li>
 *   <li>(S1, S3): diffs={1,2,2,1,2,1,1,1} -&gt; ibs2=0, ibs1=5, ibs0=3
 *       -&gt; IBS = (0+5)/16 = 5/16 = 0.3125</li>
 *   <li>(S1, S4): diffs={1,0,2,1,0,1,2,1} -&gt; ibs2=2, ibs1=4, ibs0=2
 *       -&gt; IBS = (4+4)/16 = 8/16 = 0.5</li>
 *   <li>(S2, S3): diffs={1,1,1,1,1,0,1,2} -&gt; ibs2=1, ibs1=6, ibs0=1
 *       -&gt; IBS = (2+6)/16 = 8/16 = 0.5</li>
 *   <li>(S2, S4): diffs={1,1,1,1,1,0,2,0} -&gt; ibs2=2, ibs1=5, ibs0=1
 *       -&gt; IBS = (4+5)/16 = 9/16 = 0.5625</li>
 *   <li>(S3, S4): diffs={2,2,0,0,2,0,1,2} -&gt; ibs2=3, ibs1=1, ibs0=4
 *       -&gt; IBS = (6+1)/16 = 7/16 = 0.4375</li>
 * </ul>
 */
class IbsComputerTest {

    private static final double EPS = 1e-12;

    private static byte[][] fixture() {
        int[][] gt = {
            { 1, 0, 2, 1, 0, 1, 2, 1 },
            { 1, 1, 1, 1, 1, 0, 2, 0 },
            { 2, 2, 0, 0, 2, 0, 1, 2 },
            { 0, 0, 0, 0, 0, 0, 0, 0 },
        };
        int nS = gt.length;
        int nV = gt[0].length;
        byte[][] out = new byte[nV][PackedGenotypeReader.gtPackedLength(nS)];
        for (int v = 0; v < nV; v++)
            for (int s = 0; s < nS; s++)
                PackedGenotypeReader.setGenotype(out[v], s, gt[s][v]);
        return out;
    }

    @Test
    void s1s2_isHighIBS() {
        IbsComputer.PairAccumulator p = new IbsComputer.PairAccumulator();
        for (byte[] gt : fixture()) {
            p.accumulate(PackedGenotypeReader.genotype(gt, 0),
                          PackedGenotypeReader.genotype(gt, 1));
        }
        assertEquals(3, p.ibs2);
        assertEquals(5, p.ibs1);
        assertEquals(0, p.ibs0);
        assertEquals(8, p.nUsed);
        assertEquals(11.0 / 16.0, p.ibs(), EPS);
    }

    @Test
    void s1s3_isLowIBS() {
        IbsComputer.PairAccumulator p = new IbsComputer.PairAccumulator();
        for (byte[] gt : fixture()) {
            p.accumulate(PackedGenotypeReader.genotype(gt, 0),
                          PackedGenotypeReader.genotype(gt, 2));
        }
        assertEquals(0, p.ibs2);
        assertEquals(5, p.ibs1);
        assertEquals(3, p.ibs0);
        assertEquals(5.0 / 16.0, p.ibs(), EPS);
    }

    @Test
    void s2s3_handCheckedCounters() {
        IbsComputer.PairAccumulator p = new IbsComputer.PairAccumulator();
        for (byte[] gt : fixture()) {
            p.accumulate(PackedGenotypeReader.genotype(gt, 1),
                          PackedGenotypeReader.genotype(gt, 2));
        }
        assertEquals(1, p.ibs2);
        assertEquals(6, p.ibs1);
        assertEquals(1, p.ibs0);
        assertEquals(8.0 / 16.0, p.ibs(), EPS);
    }

    @Test
    void selfPairIsOne_eachSample() {
        for (int s = 0; s < 4; s++) {
            IbsComputer.PairAccumulator p = new IbsComputer.PairAccumulator();
            for (byte[] gt : fixture()) {
                int g = PackedGenotypeReader.genotype(gt, s);
                p.accumulate(g, g);
            }
            assertEquals(1.0, p.ibs(), EPS, "self-pair S" + (s + 1));
        }
    }

    @Test
    void missingGenotypeIgnored() {
        IbsComputer.PairAccumulator p = new IbsComputer.PairAccumulator();
        byte[] v1 = new byte[1];
        PackedGenotypeReader.setGenotype(v1, 0, PackedGenotypeReader.GT_HET);
        PackedGenotypeReader.setGenotype(v1, 1, PackedGenotypeReader.GT_HET);
        byte[] v2 = new byte[1];
        PackedGenotypeReader.setGenotype(v2, 0, PackedGenotypeReader.GT_MISSING);
        PackedGenotypeReader.setGenotype(v2, 1, PackedGenotypeReader.GT_HOM_ALT);

        p.accumulate(PackedGenotypeReader.genotype(v1, 0),
                      PackedGenotypeReader.genotype(v1, 1));
        p.accumulate(PackedGenotypeReader.genotype(v2, 0),
                      PackedGenotypeReader.genotype(v2, 1));
        assertEquals(1, p.nUsed, "missing variant must be excluded");
        assertEquals(1.0, p.ibs(), EPS);
    }

    @Test
    void zeroVariants_isNaN() {
        IbsComputer.PairAccumulator p = new IbsComputer.PairAccumulator();
        assertTrue(Double.isNaN(p.ibs()));
    }

    @Test
    void computeAllPairs_returnsExpectedMatrix() {
        byte[][] gt = fixture();
        IbsComputer.Result r = IbsComputer.computeAllPairs(gt, 4,
            new int[]{0, 1, 2, 3});

        // Diagonal == 1 for every sample.
        for (int i = 0; i < 4; i++) assertEquals(1.0, r.ibs(i, i), EPS);

        // Hand-computed off-diagonals.
        assertEquals(11.0 / 16.0, r.ibs(0, 1), EPS);
        assertEquals(5.0 / 16.0, r.ibs(0, 2), EPS);
        assertEquals(8.0 / 16.0, r.ibs(0, 3), EPS);
        assertEquals(8.0 / 16.0, r.ibs(1, 2), EPS);
        assertEquals(9.0 / 16.0, r.ibs(1, 3), EPS);
        assertEquals(7.0 / 16.0, r.ibs(2, 3), EPS);

        // Symmetry.
        assertEquals(r.ibs(0, 1), r.ibs(1, 0), EPS);
        assertEquals(r.ibs(2, 3), r.ibs(3, 2), EPS);

        // Counters present.
        assertEquals(8, r.nUsed(0, 1));
        assertEquals(0, r.ibs0(0, 1));
        assertEquals(3, r.ibs2(0, 1));
    }
}
