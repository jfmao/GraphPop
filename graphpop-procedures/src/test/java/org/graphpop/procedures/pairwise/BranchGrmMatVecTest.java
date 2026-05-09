package org.graphpop.procedures.pairwise;

import org.junit.jupiter.api.Test;

import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Pure-Java unit tests for {@link BranchGrmMatVec}. Validates the
 * Algorithm V identity {@code G v == M_full · v} on hand-built ARGs
 * by comparing against {@link BranchGrmComputer} (which builds the
 * full matrix).
 */
class BranchGrmMatVecTest {

    private static final double EPS = 1e-9;

    private static ARG balancedFour() {
        // Same ARG as BranchGrmComputerTest#balancedBinaryTree.
        int nNodes = 7;
        double[] time = {0, 0, 0, 0, 0.5, 0.5, 1.0};
        boolean[] isSample = {true, true, true, true, false, false, false};
        long[] flags = {1, 1, 1, 1, 0, 0, 0};
        int[] tskitId = {0, 1, 2, 3, 4, 5, 6};
        int[] edgeParent = {4, 4, 5, 5, 6, 6};
        int[] edgeChild = {0, 1, 2, 3, 4, 5};
        long[] edgeStart = {0, 0, 0, 0, 0, 0};
        long[] edgeEnd = {100, 100, 100, 100, 100, 100};
        int[] sampleNodes = {0, 1, 2, 3};
        return new ARG("test", nNodes, 6, time, isSample, flags, tskitId,
                edgeParent, edgeChild, edgeStart, edgeEnd, sampleNodes, 100L);
    }

    private static double[] matMulVec(double[][] m, double[] v) {
        int n = v.length;
        double[] out = new double[n];
        for (int i = 0; i < n; i++) {
            double s = 0.0;
            for (int j = 0; j < n; j++) s += m[i][j] * v[j];
            out[i] = s;
        }
        return out;
    }

    @Test
    void matVec_matchesFullMatrixOnBalancedTree() {
        ARG arg = balancedFour();
        BranchGrmComputer.Result full = BranchGrmComputer.compute(arg, 0, 100);
        double[] v = {0.5, -0.3, 0.1, 0.7};
        double[] expected = matMulVec(full.matrix, v);
        double[] actual = BranchGrmMatVec.apply(arg, 0, 100, v);
        assertArrayEquals(expected, actual, EPS);
    }

    @Test
    void matVec_zeroVectorReturnsZero() {
        ARG arg = balancedFour();
        double[] gv = BranchGrmMatVec.apply(arg, 0, 100, new double[]{0, 0, 0, 0});
        for (double x : gv) assertEquals(0.0, x, 1e-15);
    }

    @Test
    void matVec_constantVectorReturnsZero() {
        // For any doubly-centred matrix B, B·1 = 0, so a constant
        // vector becomes its mean (zero after centring) and Gv = 0.
        ARG arg = balancedFour();
        double[] gv = BranchGrmMatVec.apply(arg, 0, 100, new double[]{1, 1, 1, 1});
        for (double x : gv) assertEquals(0.0, x, 1e-12);
    }

    @Test
    void matVec_rejectsWrongLength() {
        ARG arg = balancedFour();
        assertThrows(IllegalArgumentException.class,
            () -> BranchGrmMatVec.apply(arg, 0, 100, new double[]{1, 2, 3}));
    }

    @Test
    void matVec_composesWithWeightFn() {
        ARG arg = balancedFour();
        BranchWeightFn half = (a, b, c, d, e, f) -> 0.5;

        BranchGrmComputer.Result full =
                BranchGrmComputer.compute(arg, 0, 100, half);
        double[] v = {1.0, 2.0, 3.0, 4.0};
        double[] expected = matMulVec(full.matrix, v);
        double[] actual = BranchGrmMatVec.apply(arg, 0, 100, half, v);
        assertArrayEquals(expected, actual, EPS);
    }

    @Test
    void matVec_block_matchesPerVectorCalls() {
        ARG arg = balancedFour();
        Random rng = new Random(42);
        double[][] vectors = new double[3][];
        for (int i = 0; i < 3; i++) {
            vectors[i] = new double[]{rng.nextGaussian(), rng.nextGaussian(),
                    rng.nextGaussian(), rng.nextGaussian()};
        }
        double[][] block = BranchGrmMatVec.applyBlock(
                arg, 0, 100, BranchWeightFn.UNIT, vectors);
        for (int i = 0; i < 3; i++) {
            double[] one = BranchGrmMatVec.apply(
                    arg, 0, 100, BranchWeightFn.UNIT, vectors[i]);
            assertArrayEquals(one, block[i], 1e-15,
                "block[" + i + "] should match a single-vector call");
        }
    }
}
