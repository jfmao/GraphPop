package org.graphpop.procedures.recombination;

import org.graphpop.procedures.arg.ArgTraversalUtils;
import org.graphpop.procedures.pairwise.ARG;
import org.graphpop.procedures.pairwise.ARGTraversal;
import org.graphpop.procedures.selection.SelectionUtils;
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
 * Per-window ARG breakpoint density + Hudson-scaled recombination
 * rate (M11 part A).
 *
 * <pre>
 * CALL graphpop.recombination.arg_breakpoints($runId,
 *         {window_size: 10000, step: 5000})
 *   YIELD start, end, n_breakpoints, n_marginal_trees,
 *         total_branch_length, rho_per_bp, runId
 * </pre>
 *
 * <p>Algorithm per window {@code [a, b)}:</p>
 *
 * <ul>
 *   <li><b>n_breakpoints</b>: count of distinct {@code :PARENT_OF}
 *       boundary positions strictly inside {@code (a, b)} (window
 *       endpoints themselves are not breakpoints).</li>
 *   <li><b>total_branch_length</b>: span-weighted sum of branch
 *       lengths over all marginal trees overlapping the window
 *       (mirrors {@code ArgBranchDiversityProcedure}'s inner loop).</li>
 *   <li><b>rho_per_bp</b>: {@code n_breakpoints / total_branch_length}
 *       — Hudson 1983 per-bp recombination probability scaled by
 *       the lineage-time at risk.</li>
 * </ul>
 *
 * <p>Reuses {@link ArgTraversalUtils} (M6) and
 * {@link SelectionUtils#slidingWindows} (M8); no new primitives.</p>
 */
public class ArgBreakpointsProcedure {

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.recombination.arg_breakpoints",
               mode = Mode.READ)
    @Description("Per-window ARG breakpoint density + Hudson-scaled "
            + "ρ-per-bp. rho = breakpoints / total_branch_length.")
    public Stream<ArgBreakpointsResult> argBreakpoints(
            @Name("runId") String runId,
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
        if (arg.nNodes == 0) return Stream.empty();

        long[][] windows = SelectionUtils.slidingWindows(
                arg.sequenceLength, windowSize, step);
        if (windows.length == 0) return Stream.empty();

        List<ArgBreakpointsResult> rows = new ArrayList<>(windows.length);
        for (long[] win : windows) {
            long lo = win[0];
            long hi = win[1];
            // breakpoints strictly inside (lo, hi):  edges whose boundary
            // (start or end) falls in (lo, hi) and represents a transition
            // between marginal trees.
            long[] bp = arg.breakpointsClipped(lo, hi);
            long nBreakpoints = 0L;
            for (long b : bp) {
                if (b > lo && b < hi) nBreakpoints++;
            }
            long nMarginalTrees = Math.max(1L, nBreakpoints + 1L);

            double total = 0.0;
            for (int k = 0; k < bp.length - 1; k++) {
                long b0 = bp[k];
                long b1 = bp[k + 1];
                if (b1 <= b0) continue;
                int[] parent = ArgTraversalUtils.marginalTree(arg, b0);
                double treeSum = 0.0;
                for (int e = 0; e < arg.nEdges; e++) {
                    if (arg.edgeStart[e] > b0 || arg.edgeEnd[e] <= b0) {
                        continue;
                    }
                    int p = arg.edgeParent[e];
                    int c = arg.edgeChild[e];
                    treeSum += arg.time[p] - arg.time[c];
                }
                total += treeSum * (b1 - b0);
            }
            // total is in units of (generations × bp); divide by window
            // length to get per-bp lineage time. ρ_per_bp =
            // n_breakpoints / total_lineage_time_in_window.
            double rhoPerBp = (total > 0)
                    ? (double) nBreakpoints / total
                    : 0.0;

            rows.add(new ArgBreakpointsResult(
                    lo, hi, nBreakpoints, nMarginalTrees,
                    total / (hi - lo), rhoPerBp, runId));
        }
        return rows.stream();
    }

    private static long getLong(Map<String, Object> opts, String key, long def) {
        Object v = opts.get(key);
        return (v instanceof Number) ? ((Number) v).longValue() : def;
    }
}
