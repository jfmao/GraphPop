package org.graphpop.procedures.selection;

import org.graphpop.procedures.arg.ArgTraversalUtils;
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
 * Per-window branch-length outlier scan (M8 part B).
 *
 * <p>Speidel et al. 2019: a sweep in a focal sample set pulls all
 * lineages tight at the sweep position, dropping the total branch
 * length on the marginal trees overlapping the sweep window. This
 * procedure slides a window over the genome, computes total branch
 * length restricted to the focal set per window, and emits a
 * z-score against the genome-wide null. Strong negative z = sweep
 * candidate.</p>
 *
 * <pre>
 * CALL graphpop.selection.branch_outlier_scan($runId, $sampleIds,
 *         {window_size: 10000, step: 5000})
 *   YIELD start, end, total_branch_length, mean_total, sd_total,
 *         z_score, n_samples, runId
 * </pre>
 */
public class BranchOutlierScanProcedure {

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.selection.branch_outlier_scan", mode = Mode.READ)
    @Description("Per-window branch-length z-score for sweep detection. "
            + "Total branch length restricted to focal samples; null = "
            + "genome-wide (mean, sd) of per-window totals.")
    public Stream<BranchOutlierScanResult> branchOutlierScan(
            @Name("runId") String runId,
            @Name("sampleIds") List<String> sampleIds,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        if (options == null) options = new HashMap<>();
        long windowSize = getLong(options, "window_size", 10_000L);
        long step = getLong(options, "step", windowSize);
        if (windowSize <= 0 || step <= 0) {
            throw new IllegalArgumentException(
                    "window_size and step must be positive");
        }

        ARG arg = ARGTraversal.load(tx, runId, 0L, Long.MAX_VALUE);
        if (arg.nNodes == 0 || sampleIds == null || sampleIds.isEmpty()) {
            return Stream.empty();
        }
        int[] focal = resolveFocal(tx, runId, sampleIds, arg);
        if (focal.length < 2) return Stream.empty();

        long[][] windows = SelectionUtils.slidingWindows(
                arg.sequenceLength, windowSize, step);
        if (windows.length == 0) return Stream.empty();

        // Pass 1: per-window total branch length (focal-restricted).
        double[] totals = new double[windows.length];
        for (int w = 0; w < windows.length; w++) {
            totals[w] = totalBranchLengthOver(arg, focal,
                    windows[w][0], windows[w][1]);
        }

        // Genome-wide null: mean / sd of per-window totals.
        SelectionUtils.Welford null_ = new SelectionUtils.Welford();
        for (double t : totals) null_.add(t);

        List<BranchOutlierScanResult> rows = new ArrayList<>(windows.length);
        for (int w = 0; w < windows.length; w++) {
            double z = SelectionUtils.zScore(totals[w], null_);
            rows.add(new BranchOutlierScanResult(
                    windows[w][0], windows[w][1],
                    totals[w], null_.mean(), null_.sd(), z,
                    focal.length, runId));
        }
        return rows.stream();
    }

    /**
     * Sum, over marginal trees overlapping [lo, hi), of branch lengths
     * times the per-tree span fraction within the window. Restricted
     * to branches with at least one focal-sample descendant; matches
     * tskit's {@code TreeSequence.tree_sequence_total_branch_length}
     * filtered to the focal set.
     */
    private static double totalBranchLengthOver(ARG arg, int[] focal,
                                                  long lo, long hi) {
        long[] bp = arg.breakpointsClipped(lo, hi);
        if (bp.length < 2) return 0.0;
        double accum = 0.0;
        for (int k = 0; k < bp.length - 1; k++) {
            long b0 = bp[k];
            long b1 = bp[k + 1];
            if (b1 <= b0) continue;
            int[] parent = ArgTraversalUtils.marginalTree(arg, b0);
            int[] desc = ArgTraversalUtils.descendantCounts(arg, parent, focal);
            double treeSum = 0.0;
            for (int e = 0; e < arg.nEdges; e++) {
                if (arg.edgeStart[e] > b0 || arg.edgeEnd[e] <= b0) continue;
                int p = arg.edgeParent[e];
                int c = arg.edgeChild[e];
                if (desc[c] == 0) continue;
                double tb = arg.time[p] - arg.time[c];
                treeSum += tb;
            }
            long span = b1 - b0;
            accum += treeSum * span;
        }
        return accum / (hi - lo);
    }

    private static int[] resolveFocal(Transaction tx, String runId,
                                       List<String> sampleIds, ARG arg) {
        Result r = tx.execute(
                "MATCH (n:TreeNode {runId: $runId})-[:REPRESENTS]->"
              + "(s:Sample) WHERE s.sampleId IN $sids "
              + "RETURN n.nodeId AS nodeId",
                Map.of("runId", runId, "sids", sampleIds));
        List<Integer> packed = new ArrayList<>();
        try {
            Map<Integer, Integer> idMap = arg.tskitIdToIndexMap();
            while (r.hasNext()) {
                int tskitId = ((Number) r.next().get("nodeId")).intValue();
                Integer p = idMap.get(tskitId);
                if (p != null) packed.add(p);
            }
        } finally {
            r.close();
        }
        int[] out = new int[packed.size()];
        for (int i = 0; i < out.length; i++) out[i] = packed.get(i);
        return out;
    }

    private static long getLong(Map<String, Object> opts, String key, long def) {
        Object v = opts.get(key);
        return (v instanceof Number) ? ((Number) v).longValue() : def;
    }
}
