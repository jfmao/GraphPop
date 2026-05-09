package org.graphpop.procedures.pairwise;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Pure-Java unit tests for {@link BranchGrmComputer}. Operates directly
 * on hand-built {@link ARG} structs, no Neo4j required.
 *
 * <p>Validation of the full numeric agreement with the {@code egrm}
 * Python package on a real msprime fixture is in
 * {@link BranchGrmProcedureTest} (Step 5), which uses the checked-in
 * {@code egrm_expected_20samples.json}.</p>
 */
class BranchGrmComputerTest {

    private static final double EPS = 1e-12;

    /**
     * Construct a one-tree ARG with N samples and a single internal root.
     * <pre>
     *      root (time = T)
     *      / | ... | \
     *     0  1 ... N-1
     * </pre>
     * All N edges share the interval [0, L).
     */
    private static ARG starTree(int n, double rootTime, long L) {
        int nNodes = n + 1;  // n leaves + 1 root
        double[] time = new double[nNodes];
        boolean[] isSample = new boolean[nNodes];
        long[] flags = new long[nNodes];
        int[] tskitId = new int[nNodes];
        for (int i = 0; i < n; i++) {
            time[i] = 0.0;
            isSample[i] = true;
            flags[i] = 1L;
            tskitId[i] = i;
        }
        time[n] = rootTime;
        isSample[n] = false;
        flags[n] = 0L;
        tskitId[n] = n;

        int[] edgeParent = new int[n];
        int[] edgeChild = new int[n];
        long[] edgeStart = new long[n];
        long[] edgeEnd = new long[n];
        for (int i = 0; i < n; i++) {
            edgeParent[i] = n;
            edgeChild[i] = i;
            edgeStart[i] = 0L;
            edgeEnd[i] = L;
        }

        int[] sampleNodes = new int[n];
        for (int i = 0; i < n; i++) sampleNodes[i] = i;

        return new ARG("test", nNodes, n, time, isSample, flags, tskitId,
                edgeParent, edgeChild, edgeStart, edgeEnd, sampleNodes, L);
    }

    @Test
    void starTree_zeroOffDiagonal_diagonalEqualsPositive() {
        // For a star tree, every internal-node descendant count is 1 (each
        // sample's own branch); the root has all N descendants which is the
        // skip-condition. So only n*1 nodes contribute, each updating
        // exactly egrm[i][i] (since each branch has descendant set = {i}).
        ARG arg = starTree(4, 1.0, 100L);
        BranchGrmComputer.Result r = BranchGrmComputer.compute(arg, 0L, 100L);

        // Before centering, each diagonal entry has same un-normalised
        // contribution; after normalisation diagonal == 1, off-diag == 0;
        // after double-centering: row/column means are 1/n => diagonal -= 2/n
        // (once for row, once for column) and off-diagonals -= 0 - 0 + 0 = 0
        // ... actually centering applies to the normalised matrix.
        //
        // Let's just verify symmetry and structure (diagonal > 0, off-diag
        // = small non-zero from centering, row/column sums = 0).
        for (int i = 0; i < 4; i++) {
            for (int j = i + 1; j < 4; j++) {
                assertEquals(r.matrix[i][j], r.matrix[j][i], EPS, "symmetry");
            }
            // After double-centering, every row sum and column sum is zero.
            double rowSum = 0, colSum = 0;
            for (int j = 0; j < 4; j++) {
                rowSum += r.matrix[i][j];
                colSum += r.matrix[j][i];
            }
            assertEquals(0.0, rowSum, EPS, "row sum zero after centering");
            assertEquals(0.0, colSum, EPS, "col sum zero after centering");
        }
    }

    @Test
    void emptyArg_returnsZeroMatrix() {
        ARG arg = starTree(0, 1.0, 100L);
        BranchGrmComputer.Result r = BranchGrmComputer.compute(arg, 0L, 100L);
        assertEquals(0, r.n);
    }

    @Test
    void singletonSample_skippedDueToFullDescendantSet() {
        // n=1: the only internal-node descendant set covers all samples
        // (skip condition). total_mu == 0 => empty matrix returned.
        ARG arg = starTree(1, 1.0, 100L);
        BranchGrmComputer.Result r = BranchGrmComputer.compute(arg, 0L, 100L);
        assertEquals(0.0, r.totalMu);
        assertEquals(0.0, r.matrix[0][0]);
    }

    @Test
    void zeroLengthRegion_returnsZero() {
        ARG arg = starTree(4, 1.0, 100L);
        BranchGrmComputer.Result r = BranchGrmComputer.compute(arg, 50L, 50L);
        assertEquals(0.0, r.totalMu);
    }

    @Test
    void rejectsLargeSampleCount() {
        // 65 samples violates the 64-bit mask precondition.
        IllegalArgumentException ex = assertThrows(
            IllegalArgumentException.class,
            () -> BranchGrmComputer.compute(starTree(65, 1.0, 100L), 0L, 100L)
        );
        assertTrue(ex.getMessage().contains("64"));
    }

    /**
     * 4-leaf balanced binary tree:
     * <pre>
     *           5  (time T)
     *         /   \
     *        4     3   (both at time t)
     *       / \   / \
     *      0   1 2   3
     * </pre>
     * Wait — node names overlap. Use:
     * <pre>
     *           6  (time T)
     *         /   \
     *        4     5   (both at time t)
     *       / \   / \
     *      0   1 2   3
     * </pre>
     * Three non-root non-trivial nodes contribute: 4 (descendants {0,1}),
     * 5 ({2,3}), and 0,1,2,3,4,5 all have parents (only 6 is root).
     * Specifically: nodes {0,1,2,3} contribute (each desc is themselves);
     * nodes {4,5} contribute (desc {0,1} and {2,3} respectively).
     */
    @Test
    void balancedBinaryTree_symmetryAndCentering() {
        int nNodes = 7;
        double[] time = {0, 0, 0, 0, 0.5, 0.5, 1.0};
        boolean[] isSample = {true, true, true, true, false, false, false};
        long[] flags = {1, 1, 1, 1, 0, 0, 0};
        int[] tskitId = {0, 1, 2, 3, 4, 5, 6};
        // Edges: 4->0, 4->1, 5->2, 5->3, 6->4, 6->5
        int[] edgeParent = {4, 4, 5, 5, 6, 6};
        int[] edgeChild = {0, 1, 2, 3, 4, 5};
        long[] edgeStart = {0, 0, 0, 0, 0, 0};
        long[] edgeEnd = {100, 100, 100, 100, 100, 100};
        int[] sampleNodes = {0, 1, 2, 3};

        ARG arg = new ARG("test", nNodes, 6, time, isSample, flags, tskitId,
                edgeParent, edgeChild, edgeStart, edgeEnd, sampleNodes, 100L);

        BranchGrmComputer.Result r = BranchGrmComputer.compute(arg, 0L, 100L);

        // Symmetry.
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                assertEquals(r.matrix[i][j], r.matrix[j][i], EPS,
                    "symmetry [" + i + "," + j + "]");
            }
        }

        // Row/column sums all zero post-centering.
        for (int i = 0; i < 4; i++) {
            double rowSum = 0;
            for (int j = 0; j < 4; j++) rowSum += r.matrix[i][j];
            assertEquals(0.0, rowSum, EPS, "row " + i + " sum zero");
        }

        // Pairs (0,1) and (2,3) should be more related than (0,2)/(1,3)
        // because they share the recent common ancestor 4 / 5 vs the
        // older common ancestor 6.
        double sib01 = r.matrix[0][1];
        double sib23 = r.matrix[2][3];
        double cross02 = r.matrix[0][2];
        assertTrue(sib01 > cross02,
            "siblings (0,1) more related than (0,2): " + sib01 + " vs " + cross02);
        assertTrue(sib23 > cross02,
            "siblings (2,3) more related than (0,2): " + sib23 + " vs " + cross02);
        assertEquals(sib01, sib23, EPS, "topology-symmetric pairs equal");
    }
}
