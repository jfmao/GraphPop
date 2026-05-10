package org.graphpop.procedures.pairwise;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Unit tests for {@link BranchGrmLanczos}: identity tests against
 * Commons Math {@link EigenDecomposition} on the full matrix
 * produced by {@link BranchGrmComputer}.
 */
class BranchGrmLanczosTest {

    /** Balanced 4-leaf tree used elsewhere in the kinship suite. */
    private static ARG balancedFour() {
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

    /** Sorted descending eigenvalues from a full symmetric matrix. */
    private static double[] sortedEvalsDesc(double[][] m) {
        RealMatrix mr = new Array2DRowRealMatrix(m, false);
        double[] e = new EigenDecomposition(mr).getRealEigenvalues();
        java.util.Arrays.sort(e);
        // Reverse for descending order.
        double[] desc = new double[e.length];
        for (int i = 0; i < e.length; i++) desc[i] = e[e.length - 1 - i];
        return desc;
    }

    @Test
    void topEigenvalues_matchFullMatrix() {
        ARG arg = balancedFour();
        BranchGrmComputer.Result full = BranchGrmComputer.compute(arg, 0, 100);
        double[] expected = sortedEvalsDesc(full.matrix);

        // The balanced 4-sample tree has degenerate eigenvalues (within-
        // block contrast pairs). Lanczos collapses them into a single
        // direction, so top-2 distinct eigenvalues is what's testable
        // here. The integration test uses a 20-sample msprime fixture
        // with all-distinct eigenvalues.
        BranchGrmLanczos.Result r = BranchGrmLanczos.runTopK(
                arg, 0, 100, 2, 3, /*seed=*/42L);
        assertEquals(2, r.eigenvalues.length);
        // Lanczos returns the largest eigenvalue plus one of the
        // degenerate pair (which has the same value). Both must match
        // the corresponding sorted-desc reference entries.
        for (int i = 0; i < 2; i++) {
            assertEquals(expected[i], r.eigenvalues[i], 1e-9,
                "eigenvalue " + i);
        }
    }

    @Test
    void topEigenvector_matchesFullMatrix_upToSign() {
        ARG arg = balancedFour();
        BranchGrmComputer.Result full = BranchGrmComputer.compute(arg, 0, 100);
        RealMatrix mr = new Array2DRowRealMatrix(full.matrix, false);
        EigenDecomposition ed = new EigenDecomposition(mr);

        // Find index of the largest eigenvalue.
        double[] e = ed.getRealEigenvalues();
        int argMax = 0;
        for (int i = 1; i < e.length; i++) if (e[i] > e[argMax]) argMax = i;
        double[] expected = ed.getEigenvector(argMax).toArray();

        BranchGrmLanczos.Result r = BranchGrmLanczos.runTopK(
                arg, 0, 100, 1, 3, 42L);
        // |v_lanczos · v_full| should be close to 1.
        double dot = 0;
        for (int i = 0; i < expected.length; i++)
            dot += r.eigenvectors[0][i] * expected[i];
        assertTrue(Math.abs(dot) > 0.999,
            "|<v_lanczos, v_full>| = " + dot);
    }

    @Test
    void allReturnedPCsHaveZeroSum() {
        ARG arg = balancedFour();
        BranchGrmLanczos.Result r = BranchGrmLanczos.runTopK(
                arg, 0, 100, 2, 3, 42L);
        for (int j = 0; j < r.eigenvectors.length; j++) {
            double sum = 0;
            for (double x : r.eigenvectors[j]) sum += x;
            assertEquals(0.0, sum, 1e-9, "PC " + j + " sum non-zero");
        }
    }

    @Test
    void returnedPCsAreOrthogonal() {
        ARG arg = balancedFour();
        BranchGrmLanczos.Result r = BranchGrmLanczos.runTopK(
                arg, 0, 100, 2, 3, 42L);
        for (int i = 0; i < r.eigenvectors.length; i++) {
            for (int j = i + 1; j < r.eigenvectors.length; j++) {
                double dot = 0;
                for (int k = 0; k < r.eigenvectors[i].length; k++)
                    dot += r.eigenvectors[i][k] * r.eigenvectors[j][k];
                assertEquals(0.0, dot, 1e-8,
                    "PC " + i + " · PC " + j);
            }
        }
    }

    @Test
    void rejectsKEqualToOrAboveN() {
        ARG arg = balancedFour();
        assertThrows(IllegalArgumentException.class,
            () -> BranchGrmLanczos.runTopK(arg, 0, 100, 4, 4, 42L));
    }

    @Test
    void deterministicSeed() {
        ARG arg = balancedFour();
        BranchGrmLanczos.Result r1 = BranchGrmLanczos.runTopK(arg, 0, 100, 2, 4, 42L);
        BranchGrmLanczos.Result r2 = BranchGrmLanczos.runTopK(arg, 0, 100, 2, 4, 42L);
        assertArrayEquals(r1.eigenvalues, r2.eigenvalues, 1e-15);
        for (int j = 0; j < r1.eigenvectors.length; j++)
            assertArrayEquals(r1.eigenvectors[j], r2.eigenvectors[j], 1e-15);
    }
}
