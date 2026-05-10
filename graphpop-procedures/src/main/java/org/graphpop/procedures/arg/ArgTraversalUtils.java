package org.graphpop.procedures.arg;

import org.graphpop.procedures.pairwise.ARG;

import java.util.Arrays;

/**
 * Marginal-tree primitives shared by the ARG-statistics procedures (M6).
 *
 * <p>All helpers are pure functions over a packed {@link ARG}. The
 * {@code parent[]} arrays are indexed by packed node-id (see
 * {@link ARG}). They reuse the same conventions as
 * {@code BranchIbdComputer} and {@code BranchGrmComputer}.</p>
 */
public final class ArgTraversalUtils {

    private ArgTraversalUtils() {}

    /**
     * Build the marginal-tree parent[] at a given genomic position.
     * {@code parent[i]} is the packed index of {@code i}'s parent in the
     * marginal tree containing {@code position}, or {@code -1} if
     * {@code i} is a root (or the position falls outside any edge's
     * interval).
     */
    public static int[] marginalTree(ARG arg, long position) {
        int[] parent = new int[arg.nNodes];
        Arrays.fill(parent, -1);
        for (int e = 0; e < arg.nEdges; e++) {
            if (arg.edgeStart[e] <= position && arg.edgeEnd[e] > position) {
                parent[arg.edgeChild[e]] = arg.edgeParent[e];
            }
        }
        return parent;
    }

    /**
     * Count focal-sample descendants at every packed node in the marginal
     * tree defined by {@code parent[]}. {@code counts[i]} is the number
     * of focal nodes whose ancestor chain passes through {@code i}
     * (each focal node itself is counted at its own position).
     *
     * <p>Focal nodes outside the marginal tree (i.e., with
     * {@code parent[focal] == -1} and {@code focal} is not actually a
     * leaf of this tree) are still walked from themselves only; this
     * mirrors how IBD/branch-GRM iterate samples.</p>
     */
    public static int[] descendantCounts(ARG arg, int[] parent, int[] focalNodes) {
        int[] counts = new int[arg.nNodes];
        for (int leaf : focalNodes) {
            int cur = leaf;
            while (cur != -1) {
                counts[cur]++;
                cur = parent[cur];
            }
        }
        return counts;
    }

    /**
     * MRCA packed index of two leaves in a marginal tree. Returns
     * {@code -1} if no common ancestor exists in the tree.
     */
    public static int mrca(int leafA, int leafB, int[] parent) {
        boolean[] seen = new boolean[parent.length];
        int cur = leafA;
        while (cur != -1) { seen[cur] = true; cur = parent[cur]; }
        int result = -1;
        cur = leafB;
        while (cur != -1) {
            if (seen[cur]) { result = cur; break; }
            cur = parent[cur];
        }
        return result;
    }

    /**
     * Convenience: sorted distinct edge-interval breakpoints over the
     * full sequence span. The marginal tree is constant on every
     * sub-interval {@code [breakpoints[k], breakpoints[k+1])}.
     */
    public static long[] breakpoints(ARG arg) {
        long hi = (arg.sequenceLength == Long.MAX_VALUE) ? Long.MAX_VALUE
                : arg.sequenceLength;
        return arg.breakpointsClipped(0L, hi);
    }
}
