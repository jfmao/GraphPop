package org.graphpop.procedures.pairwise;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * ARG-derived IBD-segment caller.
 *
 * <p>Sweep-line over marginal-tree intervals. For each pair of
 * sample-flagged TreeNodes, tracks the current MRCA across breakpoints
 * and emits a segment whenever the MRCA changes (or whenever the
 * MRCA's time exceeds {@code max_tmrca}).</p>
 *
 * <p>Algorithmic equivalence with
 * {@code tskit.TreeSequence.ibd_segments(within=samples,
 *                                          max_time=max_tmrca,
 *                                          min_span=min_length_bp,
 *                                          store_pairs=True,
 *                                          store_segments=True)}
 * is the headline correctness gate (see {@code BranchIbdProcedureTest}).</p>
 */
public final class BranchIbdComputer {

    private BranchIbdComputer() {}

    public static final double NO_MAX_TMRCA = Double.POSITIVE_INFINITY;
    public static final long NO_MIN_LENGTH = 0L;

    /**
     * Compute IBD segments for every pair of sample-flagged
     * TreeNodes in the ARG.
     *
     * @param arg            loaded ARG
     * @param regionStart    bp inclusive
     * @param regionEnd      bp exclusive ({@link Long#MAX_VALUE} for full)
     * @param maxTmrca       drop segments whose MRCA time exceeds this
     *                       (in generations); {@link #NO_MAX_TMRCA} for none
     * @param minLengthBp    drop segments shorter than this
     */
    public static List<IbdSegment> compute(ARG arg,
                                            long regionStart, long regionEnd,
                                            double maxTmrca, long minLengthBp) {
        final int n = arg.nSamples();
        if (n < 2) return List.of();
        final long lo = Math.max(regionStart, 0L);
        final long hi = (regionEnd == Long.MAX_VALUE) ? arg.sequenceLength
                : Math.min(regionEnd, arg.sequenceLength);
        if (hi <= lo) return List.of();

        final long[] breakpoints = arg.breakpointsClipped(lo, hi);
        if (breakpoints.length < 2) return List.of();

        final int nNodes = arg.nNodes;
        final int nIntervals = breakpoints.length - 1;

        // For each interval k, build the marginal-tree parent[] once and
        // store; needed to scan O(N^2) pairs per interval.
        int[][] parentByInterval = new int[nIntervals][nNodes];
        for (int k = 0; k < nIntervals; k++) {
            long b0 = breakpoints[k];
            int[] parent = parentByInterval[k];
            Arrays.fill(parent, -1);
            for (int e = 0; e < arg.nEdges; e++) {
                if (arg.edgeStart[e] <= b0 && arg.edgeEnd[e] > b0) {
                    parent[arg.edgeChild[e]] = arg.edgeParent[e];
                }
            }
        }

        // Sample tskit indices (packed indexing already maps to [0, n)).
        // For canonical pair iteration we use the packed sample index,
        // converted to/from tskit node id at emit time.
        final int[] sampleNodes = arg.sampleNodes;

        // Per-pair "open" segment state. Flatten by (a, b).
        // openMrca[idx] = -1   ⇒ no open segment.
        // openPath[idx] = the path-signature for the currently-open segment
        //   (i_chain ++ reverse(j_chain[1:])). A segment ends when this
        //   signature changes -- matching tskit.TreeSequence.ibd_segments,
        //   which splits on path identity (not just MRCA identity).
        final int[] openMrca = new int[n * n];
        final long[] openStart = new long[n * n];
        final int[][] openPath = new int[n * n][];
        Arrays.fill(openMrca, -1);

        List<IbdSegment> out = new ArrayList<>();
        boolean[] seen = new boolean[nNodes];

        for (int k = 0; k < nIntervals; k++) {
            long b0 = breakpoints[k];
            long b1 = breakpoints[k + 1];
            int[] parent = parentByInterval[k];

            for (int b = 1; b < n; b++) {
                int sB = sampleNodes[b];
                for (int a = 0; a < b; a++) {
                    int sA = sampleNodes[a];

                    int mrca = findMRCA(sA, sB, parent, seen);
                    double mrcaTime = (mrca >= 0)
                            ? arg.time[mrca]
                            : Double.POSITIVE_INFINITY;
                    boolean acceptable = mrca >= 0 && mrcaTime <= maxTmrca;
                    int[] path = acceptable
                            ? buildPath(sA, sB, mrca, parent)
                            : null;

                    int idx = a * n + b;
                    int prevMrca = openMrca[idx];
                    int[] prevPath = openPath[idx];

                    boolean pathChanged = (prevPath == null) != (path == null)
                            || (path != null && !Arrays.equals(path, prevPath));

                    if (pathChanged) {
                        if (prevMrca != -1) {
                            long segStart = openStart[idx];
                            long segEnd = b0;
                            if (segEnd - segStart >= minLengthBp) {
                                out.add(new IbdSegment(a, b, segStart, segEnd,
                                        arg.tskitNodeId[prevMrca],
                                        arg.time[prevMrca]));
                            }
                        }
                        if (acceptable) {
                            openMrca[idx] = mrca;
                            openStart[idx] = b0;
                            openPath[idx] = path;
                        } else {
                            openMrca[idx] = -1;
                            openPath[idx] = null;
                        }
                    }
                    // else: same path continues -- nothing to do.
                }
            }
        }

        // Close any segments still open at hi.
        long lastEnd = breakpoints[breakpoints.length - 1];
        for (int b = 1; b < n; b++) {
            for (int a = 0; a < b; a++) {
                int idx = a * n + b;
                int mrca = openMrca[idx];
                if (mrca != -1) {
                    long segStart = openStart[idx];
                    long segEnd = lastEnd;
                    if (segEnd - segStart >= minLengthBp) {
                        out.add(new IbdSegment(a, b, segStart, segEnd,
                                arg.tskitNodeId[mrca],
                                arg.time[mrca]));
                    }
                }
            }
        }
        return out;
    }

    /**
     * Find the MRCA of two leaf nodes in a single marginal tree
     * defined by {@code parent[]}. {@code seen} is reused scratch
     * (cleared internally). Returns -1 if no common ancestor.
     */
    private static int findMRCA(int a, int b, int[] parent, boolean[] seen) {
        int cur = a;
        while (cur != -1) {
            seen[cur] = true;
            cur = parent[cur];
        }
        int mrca = -1;
        cur = b;
        while (cur != -1) {
            if (seen[cur]) { mrca = cur; break; }
            cur = parent[cur];
        }
        cur = a;
        while (cur != -1) {
            seen[cur] = false;
            cur = parent[cur];
        }
        return mrca;
    }

    /**
     * Build the path identifier from leaf {@code a} up to {@code mrca}
     * and back down to leaf {@code b}. The resulting array uniquely
     * identifies the topology between the two leaves under the
     * current marginal tree. Two intervals share an IBD segment iff
     * they produce the same path.
     */
    private static int[] buildPath(int a, int b, int mrca, int[] parent) {
        // Walk from a up to mrca.
        int lenA = 0;
        int cur = a;
        while (cur != mrca) { lenA++; cur = parent[cur]; }
        int lenB = 0;
        cur = b;
        while (cur != mrca) { lenB++; cur = parent[cur]; }
        // Path: a -> ... -> mrca -> ... -> b
        int[] path = new int[lenA + 1 + lenB];
        cur = a;
        for (int i = 0; i < lenA; i++) {
            path[i] = cur;
            cur = parent[cur];
        }
        path[lenA] = mrca;
        // Walk b's chain into a buffer, then reverse-fill the tail.
        int[] bChain = new int[lenB];
        cur = b;
        for (int i = 0; i < lenB; i++) {
            bChain[i] = cur;
            cur = parent[cur];
        }
        for (int i = 0; i < lenB; i++) {
            path[lenA + 1 + i] = bChain[lenB - 1 - i];
        }
        return path;
    }
}
