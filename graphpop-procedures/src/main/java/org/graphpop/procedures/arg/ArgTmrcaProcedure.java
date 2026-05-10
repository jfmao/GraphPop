package org.graphpop.procedures.arg;


import org.graphpop.procedures.pairwise.ARG;
import org.graphpop.procedures.pairwise.ARGTraversal;
import org.neo4j.graphdb.Result;
import org.neo4j.graphdb.Transaction;
import org.neo4j.procedure.Context;
import org.neo4j.procedure.Description;
import org.neo4j.procedure.Mode;
import org.neo4j.procedure.Name;
import org.neo4j.procedure.Procedure;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

/**
 * TMRCA between two samples in an ARG (M6 part A).
 *
 * <pre>
 * CALL graphpop.arg.tmrca($runId, $sampleA, $sampleB,
 *                          {position: 1_000_000})
 *   YIELD sample_a, sample_b, position, tmrca, mean_tmrca,
 *         mrca_node_id, runId
 * </pre>
 *
 * <p>Two modes:</p>
 *
 * <ul>
 *   <li><b>Single-position</b> — {@code options.position} set: returns
 *       one row per haplotype-pair with the LCA time at that position
 *       in the marginal tree.</li>
 *   <li><b>Window-mean</b> — {@code options.window_start} /
 *       {@code window_end} set (or both omitted for whole-genome
 *       mean): returns one row per pair with the span-weighted mean
 *       LCA time across overlapping marginal trees.</li>
 * </ul>
 *
 * <p>If a sample has multiple haplotypes (i.e. multiple
 * {@code :REPRESENTS} edges), every haplotype pair is enumerated.</p>
 */
public class ArgTmrcaProcedure {

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.arg.tmrca", mode = Mode.READ)
    @Description("TMRCA between two samples in an ARG. "
            + "options.position for single-position; "
            + "options.window_start/window_end for window-mean; "
            + "no options for whole-genome span-weighted mean.")
    public Stream<TmrcaResult> tmrca(
            @Name("runId") String runId,
            @Name("sampleA") String sampleA,
            @Name("sampleB") String sampleB,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        if (options == null) options = new HashMap<>();
        Long position = getNullableLong(options, "position");
        Long windowStart = getNullableLong(options, "window_start");
        Long windowEnd = getNullableLong(options, "window_end");

        ARG arg = ARGTraversal.load(tx, runId, 0L, Long.MAX_VALUE);
        if (arg.nNodes == 0) return Stream.empty();

        // Resolve sample-id -> packed leaf indices.
        int[] leavesA = leafIndicesForSample(tx, runId, sampleA, arg);
        int[] leavesB = leafIndicesForSample(tx, runId, sampleB, arg);
        if (leavesA.length == 0 || leavesB.length == 0) return Stream.empty();

        List<TmrcaResult> rows = new ArrayList<>();

        if (position != null) {
            int[] parent = ArgTraversalUtils.marginalTree(arg, position);
            for (int la : leavesA) {
                for (int lb : leavesB) {
                    int m = ArgTraversalUtils.mrca(la, lb, parent);
                    if (m < 0) continue;
                    rows.add(new TmrcaResult(
                            sampleA, sampleB, position,
                            arg.time[m], Double.NaN,
                            arg.tskitNodeId[m], runId));
                }
            }
            return rows.stream();
        }

        long lo = (windowStart != null) ? windowStart : 0L;
        long hi = (windowEnd != null) ? windowEnd : arg.sequenceLength;
        if (hi <= lo) return Stream.empty();

        long[] bp = arg.breakpointsClipped(lo, hi);
        for (int la : leavesA) {
            for (int lb : leavesB) {
                double weightedSum = 0.0;
                long span = 0L;
                long lastMrcaTskit = -1L;
                for (int k = 0; k < bp.length - 1; k++) {
                    long b0 = bp[k], b1 = bp[k + 1];
                    if (b1 <= b0) continue;
                    int[] parent = ArgTraversalUtils.marginalTree(arg, b0);
                    int m = ArgTraversalUtils.mrca(la, lb, parent);
                    if (m < 0) continue;
                    double t = arg.time[m];
                    long s = b1 - b0;
                    weightedSum += t * s;
                    span += s;
                    lastMrcaTskit = arg.tskitNodeId[m];
                }
                if (span == 0) continue;
                rows.add(new TmrcaResult(
                        sampleA, sampleB, -1L,
                        Double.NaN, weightedSum / span,
                        lastMrcaTskit, runId));
            }
        }
        return rows.stream();
    }

    private static int[] leafIndicesForSample(Transaction tx, String runId,
                                                String sampleId, ARG arg) {
        // Walk REPRESENTS to find every TreeNode tied to this sample
        // for this run, then map tskit node id -> packed index.
        List<Integer> leaves = new ArrayList<>();
        Result r = tx.execute(
                "MATCH (n:TreeNode {runId: $runId})-[:REPRESENTS]->"
              + "(:Sample {sampleId: $sid}) "
              + "RETURN n.nodeId AS nodeId",
                Map.of("runId", runId, "sid", sampleId));
        try {
            Map<Integer, Integer> idMap = arg.tskitIdToIndexMap();
            while (r.hasNext()) {
                int tskitId = ((Number) r.next().get("nodeId")).intValue();
                Integer packed = idMap.get(tskitId);
                if (packed != null) leaves.add(packed);
            }
        } finally {
            r.close();
        }
        int[] out = new int[leaves.size()];
        for (int i = 0; i < out.length; i++) out[i] = leaves.get(i);
        return out;
    }

    private static Long getNullableLong(Map<String, Object> opts, String key) {
        Object v = opts.get(key);
        return (v instanceof Number) ? ((Number) v).longValue() : null;
    }
}
