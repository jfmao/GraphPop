package org.graphpop.procedures.pairwise;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Packed in-memory ancestral recombination graph for a single
 * {@code :ARGRun}. Loaded by {@link ARGTraversal#load} and consumed by
 * {@link BranchGrmComputer} and the conditional / posterior variants.
 *
 * <p>Indexing convention: every node has a packed index in
 * {@code [0, nNodes)} that mirrors insertion order from Cypher; edges
 * reference parent/child by packed index (not tskit node-id).</p>
 */
public final class ARG {

    public final int nNodes;
    public final int nEdges;

    /** Time-of-node (generations). Length {@link #nNodes}. */
    public final double[] time;

    /** True when the node is a sample-flagged leaf. Length {@link #nNodes}. */
    public final boolean[] isSample;

    /** tskit node flags. Length {@link #nNodes}. */
    public final long[] flags;

    /**
     * tskit local node-id stored on the {@code :TreeNode}. Length
     * {@link #nNodes}.
     */
    public final int[] tskitNodeId;

    // ---- edges (parallel arrays of length nEdges) ----

    /** Packed index of edge parent. */
    public final int[] edgeParent;

    /** Packed index of edge child. */
    public final int[] edgeChild;

    /** Edge interval start (bp). */
    public final long[] edgeStart;

    /** Edge interval end (bp). */
    public final long[] edgeEnd;

    /**
     * Packed indices of sample-flagged nodes, in tskit's natural order
     * (haplotype 0, haplotype 1, ... for diploid samples).
     */
    public final int[] sampleNodes;

    /**
     * Run-id this ARG was loaded from.
     */
    public final String runId;

    /**
     * Sequence length covered by this ARG (bp). May be {@link Long#MAX_VALUE}
     * if loaded over an open interval.
     */
    public final long sequenceLength;

    public ARG(String runId, int nNodes, int nEdges,
               double[] time, boolean[] isSample, long[] flags, int[] tskitNodeId,
               int[] edgeParent, int[] edgeChild,
               long[] edgeStart, long[] edgeEnd,
               int[] sampleNodes, long sequenceLength) {
        this.runId = runId;
        this.nNodes = nNodes;
        this.nEdges = nEdges;
        this.time = time;
        this.isSample = isSample;
        this.flags = flags;
        this.tskitNodeId = tskitNodeId;
        this.edgeParent = edgeParent;
        this.edgeChild = edgeChild;
        this.edgeStart = edgeStart;
        this.edgeEnd = edgeEnd;
        this.sampleNodes = sampleNodes;
        this.sequenceLength = sequenceLength;
    }

    /** Number of sample-flagged nodes (haplotypes). */
    public int nSamples() { return sampleNodes.length; }

    /**
     * Build a map from tskit node-id to packed index. O(nNodes) space; called
     * lazily by callers that need the inverse mapping.
     */
    public Map<Integer, Integer> tskitIdToIndexMap() {
        Map<Integer, Integer> m = new HashMap<>(nNodes * 2);
        for (int i = 0; i < nNodes; i++) m.put(tskitNodeId[i], i);
        return m;
    }

    /**
     * Sorted unique breakpoints (edge starts and ends) clipped to {@code [lo, hi]}.
     * Always includes {@code lo} as the first element. Used as the sweep-line
     * pivot for marginal-tree iteration.
     */
    public long[] breakpointsClipped(long lo, long hi) {
        long[] raw = new long[nEdges * 2];
        int n = 0;
        for (int e = 0; e < nEdges; e++) {
            long s = edgeStart[e], t = edgeEnd[e];
            if (s >= lo && s <= hi) raw[n++] = s;
            if (t >= lo && t <= hi) raw[n++] = t;
        }
        long[] trimmed = Arrays.copyOf(raw, n);
        Arrays.sort(trimmed);
        // Deduplicate.
        int w = 0;
        for (int i = 0; i < n; i++) {
            if (w == 0 || trimmed[i] != trimmed[w - 1]) {
                trimmed[w++] = trimmed[i];
            }
        }
        // Trim to actual length.
        long[] dedup = Arrays.copyOf(trimmed, w);
        // Ensure lo prepended.
        if (w == 0 || dedup[0] != lo) {
            long[] withLo = new long[w + 1];
            withLo[0] = lo;
            System.arraycopy(dedup, 0, withLo, 1, w);
            dedup = withLo;
            w++;
        }
        // Ensure hi appended.
        if (dedup[w - 1] != hi) {
            long[] withHi = Arrays.copyOf(dedup, w + 1);
            withHi[w] = hi;
            dedup = withHi;
        }
        return dedup;
    }
}
