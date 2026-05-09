package org.graphpop.procedures.pairwise;

import org.graphpop.procedures.PackedGenotypeReader;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Hand-computed unit tests for KING-robust kinship.
 *
 * <p>Manichaikul et al. 2010 eq. 9:
 * {@code phi_ij = (HetHet - 2*IBS0) / (2 * min(N_Aa(i), N_Aa(j)))}.</p>
 *
 * <p>Fixture: 4 samples x 8 variants. Genotype matrix (rows=samples, cols=variants):
 * <pre>
 *      V1 V2 V3 V4 V5 V6 V7 V8
 * S1:  1  0  2  1  0  1  2  1     N_Aa(S1) = 4
 * S2:  1  1  1  1  1  0  2  0     N_Aa(S2) = 5
 * S3:  2  2  0  0  2  0  1  2     N_Aa(S3) = 1
 * S4:  0  0  0  0  0  0  0  0     N_Aa(S4) = 0
 * </pre>
 * Expected (hand-computed):
 * <ul>
 *   <li>(S1,S2): HetHet=2 (V1,V4), IBS0=0 -&gt; phi = (2-0)/(2*4) = 0.25</li>
 *   <li>(S1,S3): HetHet=0, IBS0=3 (V2,V3,V5) -&gt; phi = (0-6)/(2*1) = -3.0</li>
 *   <li>(S2,S3): HetHet=0, IBS0=1 (V8) -&gt; phi = -2/(2*1) = -1.0</li>
 *   <li>(S1,S4),(S2,S4),(S3,S4): min N_Aa = 0 -&gt; NaN (undefined)</li>
 * </ul>
 */
class KingRobustComputerTest {

    private static final double EPS = 1e-12;

    /** Build the 4-sample x 8-variant fixture as one packed byte[] per variant. */
    private static byte[][] fixtureGtPacked() {
        int[][] gt = {
            // V1, V2, V3, V4, V5, V6, V7, V8
            { 1, 0, 2, 1, 0, 1, 2, 1 },  // S1
            { 1, 1, 1, 1, 1, 0, 2, 0 },  // S2
            { 2, 2, 0, 0, 2, 0, 1, 2 },  // S3
            { 0, 0, 0, 0, 0, 0, 0, 0 },  // S4
        };
        int nSamples = gt.length;
        int nVariants = gt[0].length;
        byte[][] out = new byte[nVariants][PackedGenotypeReader.gtPackedLength(nSamples)];
        for (int v = 0; v < nVariants; v++) {
            for (int s = 0; s < nSamples; s++) {
                PackedGenotypeReader.setGenotype(out[v], s, gt[s][v]);
            }
        }
        return out;
    }

    @Test
    void pairStats_S1S2_parentChild() {
        KingRobustComputer.PairAccumulator pair = new KingRobustComputer.PairAccumulator();
        for (byte[] gtPacked : fixtureGtPacked()) {
            int gtI = PackedGenotypeReader.genotype(gtPacked, 0);
            int gtJ = PackedGenotypeReader.genotype(gtPacked, 1);
            pair.accumulate(gtI, gtJ);
        }
        assertEquals(2, pair.hetHet);
        assertEquals(0, pair.ibs0);
        assertEquals(4, pair.nAaI);
        assertEquals(5, pair.nAaJ);
        assertEquals(8, pair.nUsed);
        assertEquals(0.25, pair.phi(), EPS);
    }

    @Test
    void pairStats_S1S3_unrelatedDistant() {
        KingRobustComputer.PairAccumulator pair = new KingRobustComputer.PairAccumulator();
        for (byte[] gtPacked : fixtureGtPacked()) {
            int gtI = PackedGenotypeReader.genotype(gtPacked, 0);
            int gtJ = PackedGenotypeReader.genotype(gtPacked, 2);
            pair.accumulate(gtI, gtJ);
        }
        assertEquals(0, pair.hetHet);
        assertEquals(3, pair.ibs0);
        assertEquals(4, pair.nAaI);
        assertEquals(1, pair.nAaJ);
        assertEquals(-3.0, pair.phi(), EPS);
    }

    @Test
    void pairStats_S2S3() {
        KingRobustComputer.PairAccumulator pair = new KingRobustComputer.PairAccumulator();
        for (byte[] gtPacked : fixtureGtPacked()) {
            int gtI = PackedGenotypeReader.genotype(gtPacked, 1);
            int gtJ = PackedGenotypeReader.genotype(gtPacked, 2);
            pair.accumulate(gtI, gtJ);
        }
        assertEquals(0, pair.hetHet);
        assertEquals(1, pair.ibs0);
        assertEquals(5, pair.nAaI);
        assertEquals(1, pair.nAaJ);
        assertEquals(-1.0, pair.phi(), EPS);
    }

    @Test
    void pairStats_anyPairWithS4_undefined() {
        // S4 has zero heterozygous sites => min(N_Aa) = 0 => phi is NaN
        for (int j : new int[]{0, 1, 2}) {
            KingRobustComputer.PairAccumulator pair = new KingRobustComputer.PairAccumulator();
            for (byte[] gtPacked : fixtureGtPacked()) {
                int gtI = PackedGenotypeReader.genotype(gtPacked, j);
                int gtJ = PackedGenotypeReader.genotype(gtPacked, 3);
                pair.accumulate(gtI, gtJ);
            }
            assertEquals(0, pair.nAaJ, "S4 must have zero het sites");
            assertTrue(Double.isNaN(pair.phi()),
                "phi must be NaN when min(N_Aa) is zero (pair S" + (j + 1) + ", S4)");
        }
    }

    @Test
    void pairStats_selfPair_isHalf() {
        // KING formula applied to a sample against itself:
        //   HetHet = N_Aa, IBS0 = 0  ->  phi = N_Aa / (2*N_Aa) = 0.5
        KingRobustComputer.PairAccumulator pair = new KingRobustComputer.PairAccumulator();
        for (byte[] gtPacked : fixtureGtPacked()) {
            int gt = PackedGenotypeReader.genotype(gtPacked, 0);
            pair.accumulate(gt, gt);
        }
        assertEquals(4, pair.hetHet);
        assertEquals(0, pair.ibs0);
        assertEquals(0.5, pair.phi(), EPS);
    }

    @Test
    void pairStats_missingGenotypesIgnored() {
        // Variants where either sample is missing should not affect counters.
        KingRobustComputer.PairAccumulator pair = new KingRobustComputer.PairAccumulator();
        // Two variants: first both het, second one missing.
        byte[] v1 = new byte[1];
        PackedGenotypeReader.setGenotype(v1, 0, PackedGenotypeReader.GT_HET);
        PackedGenotypeReader.setGenotype(v1, 1, PackedGenotypeReader.GT_HET);
        byte[] v2 = new byte[1];
        PackedGenotypeReader.setGenotype(v2, 0, PackedGenotypeReader.GT_MISSING);
        PackedGenotypeReader.setGenotype(v2, 1, PackedGenotypeReader.GT_HOM_ALT);

        for (byte[] gtPacked : new byte[][]{v1, v2}) {
            int gtI = PackedGenotypeReader.genotype(gtPacked, 0);
            int gtJ = PackedGenotypeReader.genotype(gtPacked, 1);
            pair.accumulate(gtI, gtJ);
        }
        assertEquals(1, pair.hetHet);
        assertEquals(0, pair.ibs0);
        assertEquals(1, pair.nAaI);
        assertEquals(1, pair.nAaJ);
        assertEquals(1, pair.nUsed, "only V1 is fully called");
        assertEquals(0.5, pair.phi(), EPS);
    }

    // ---- Higher-level API: computeAllPairs over a sample matrix ----

    @Test
    void computeAllPairs_returnsExpectedMatrix() {
        byte[][] gtPackedPerVariant = fixtureGtPacked();
        int nSamples = 4;
        int[] packedIndices = {0, 1, 2, 3};

        KingRobustComputer.Result result =
            KingRobustComputer.computeAllPairs(gtPackedPerVariant, nSamples, packedIndices);

        // Diagonal: 0.5 for samples with at least one het site; NaN for S4 (no het sites)
        assertEquals(0.5, result.phi(0, 0), EPS);
        assertEquals(0.5, result.phi(1, 1), EPS);
        assertEquals(0.5, result.phi(2, 2), EPS);
        assertTrue(Double.isNaN(result.phi(3, 3)), "S4 has no het sites; diagonal must be NaN");
        // Off-diagonal hand-computed values
        assertEquals(0.25, result.phi(0, 1), EPS);
        assertEquals(0.25, result.phi(1, 0), EPS, "phi must be symmetric");
        assertEquals(-3.0, result.phi(0, 2), EPS);
        assertEquals(-1.0, result.phi(1, 2), EPS);
        // Pairs involving S4 are NaN
        assertTrue(Double.isNaN(result.phi(0, 3)));
        assertTrue(Double.isNaN(result.phi(1, 3)));
        assertTrue(Double.isNaN(result.phi(2, 3)));

        assertEquals(8, result.nUsed(0, 1));
        assertEquals(8, result.nUsed(0, 2));
    }
}
