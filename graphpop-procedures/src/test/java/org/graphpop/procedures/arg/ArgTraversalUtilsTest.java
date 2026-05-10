package org.graphpop.procedures.arg;

import org.graphpop.procedures.pairwise.ARG;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Unit tests on a hand-built 4-leaf balanced tree:
 *
 * <pre>
 *           7  (t=2.0)
 *          / \
 *         5   6  (t=1.0)
 *        /|   |\
 *       0 1   2 3  (t=0)
 * </pre>
 *
 * Edges all span [0, 100]. Packed indices == tskit ids.
 */
class ArgTraversalUtilsTest {

    private static ARG fourLeafTree() {
        int n = 8;
        double[] time   = {0, 0, 0, 0, 0, 1.0, 1.0, 2.0};
        boolean[] isSam = {true, true, true, true, false, false, false, false};
        long[] flags    = {1, 1, 1, 1, 0, 0, 0, 0};
        int[] tskitIds  = {0, 1, 2, 3, 4, 5, 6, 7};
        // Note: index 4 unused (no node 4 in the tree); keep nNodes=8
        // for symmetry with tskit's id space, but ignore it.
        // 6 edges: (5,0), (5,1), (6,2), (6,3), (7,5), (7,6).
        int[] eParent = {5, 5, 6, 6, 7, 7};
        int[] eChild  = {0, 1, 2, 3, 5, 6};
        long[] eStart = {0, 0, 0, 0, 0, 0};
        long[] eEnd   = {100, 100, 100, 100, 100, 100};
        int[] sampleNodes = {0, 1, 2, 3};
        return new ARG("test", n, eParent.length,
                time, isSam, flags, tskitIds,
                eParent, eChild, eStart, eEnd,
                sampleNodes, 100L);
    }

    @Test
    void marginalTree_inside_span_returns_full_tree() {
        ARG arg = fourLeafTree();
        int[] parent = ArgTraversalUtils.marginalTree(arg, 50);
        assertEquals(5, parent[0]);
        assertEquals(5, parent[1]);
        assertEquals(6, parent[2]);
        assertEquals(6, parent[3]);
        assertEquals(7, parent[5]);
        assertEquals(7, parent[6]);
        assertEquals(-1, parent[7]);  // root
    }

    @Test
    void marginalTree_outside_span_all_disconnected() {
        ARG arg = fourLeafTree();
        int[] parent = ArgTraversalUtils.marginalTree(arg, 200);
        for (int p : parent) assertEquals(-1, p);
    }

    @Test
    void marginalTree_at_left_endpoint_inclusive() {
        ARG arg = fourLeafTree();
        int[] parent = ArgTraversalUtils.marginalTree(arg, 0);
        assertEquals(5, parent[0]);  // edges have start <= position
    }

    @Test
    void marginalTree_at_right_endpoint_exclusive() {
        ARG arg = fourLeafTree();
        int[] parent = ArgTraversalUtils.marginalTree(arg, 100);
        // end=100 means position 100 is OUTSIDE the interval [0,100).
        assertEquals(-1, parent[0]);
    }

    @Test
    void mrca_within_left_subtree() {
        ARG arg = fourLeafTree();
        int[] parent = ArgTraversalUtils.marginalTree(arg, 50);
        assertEquals(5, ArgTraversalUtils.mrca(0, 1, parent));
    }

    @Test
    void mrca_within_right_subtree() {
        ARG arg = fourLeafTree();
        int[] parent = ArgTraversalUtils.marginalTree(arg, 50);
        assertEquals(6, ArgTraversalUtils.mrca(2, 3, parent));
    }

    @Test
    void mrca_across_subtrees_is_root() {
        ARG arg = fourLeafTree();
        int[] parent = ArgTraversalUtils.marginalTree(arg, 50);
        assertEquals(7, ArgTraversalUtils.mrca(0, 2, parent));
        assertEquals(7, ArgTraversalUtils.mrca(1, 3, parent));
    }

    @Test
    void mrca_self_returns_self() {
        ARG arg = fourLeafTree();
        int[] parent = ArgTraversalUtils.marginalTree(arg, 50);
        assertEquals(0, ArgTraversalUtils.mrca(0, 0, parent));
    }

    @Test
    void mrca_when_no_tree_returns_minus_one() {
        ARG arg = fourLeafTree();
        int[] parent = ArgTraversalUtils.marginalTree(arg, 200);
        assertEquals(-1, ArgTraversalUtils.mrca(0, 1, parent));
    }

    @Test
    void descendantCounts_full_focal_set() {
        ARG arg = fourLeafTree();
        int[] parent = ArgTraversalUtils.marginalTree(arg, 50);
        int[] counts = ArgTraversalUtils.descendantCounts(
                arg, parent, new int[]{0, 1, 2, 3});
        assertEquals(1, counts[0]);
        assertEquals(1, counts[1]);
        assertEquals(1, counts[2]);
        assertEquals(1, counts[3]);
        assertEquals(2, counts[5]);
        assertEquals(2, counts[6]);
        assertEquals(4, counts[7]);
    }

    @Test
    void descendantCounts_partial_focal_set() {
        ARG arg = fourLeafTree();
        int[] parent = ArgTraversalUtils.marginalTree(arg, 50);
        int[] counts = ArgTraversalUtils.descendantCounts(
                arg, parent, new int[]{0, 2});
        assertEquals(1, counts[0]);
        assertEquals(0, counts[1]);
        assertEquals(1, counts[2]);
        assertEquals(0, counts[3]);
        assertEquals(1, counts[5]);
        assertEquals(1, counts[6]);
        assertEquals(2, counts[7]);
    }

    @Test
    void descendantCounts_singleton_focal() {
        ARG arg = fourLeafTree();
        int[] parent = ArgTraversalUtils.marginalTree(arg, 50);
        int[] counts = ArgTraversalUtils.descendantCounts(
                arg, parent, new int[]{0});
        assertEquals(1, counts[0]);
        assertEquals(1, counts[5]);
        assertEquals(1, counts[7]);
        assertEquals(0, counts[1]);
        assertEquals(0, counts[2]);
        assertEquals(0, counts[6]);
    }

    @Test
    void breakpoints_returns_zero_and_sequence_length() {
        ARG arg = fourLeafTree();
        long[] bp = ArgTraversalUtils.breakpoints(arg);
        assertEquals(0L, bp[0]);
        assertEquals(100L, bp[bp.length - 1]);
    }
}
