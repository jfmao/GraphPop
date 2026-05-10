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
 * Branch-mode diversity from an ARG (M6 part B).
 *
 * <pre>
 * CALL graphpop.arg.branch_diversity($runId, $sampleIds,
 *                                     {mode: 'pi', windows: null})
 *   YIELD start, end, branch_pi, n_samples, mode, runId
 * </pre>
 *
 * <p>For {@code mode = 'pi'} (the v1 default) reproduces
 * {@code tskit.TreeSequence.diversity(samples, mode='branch',
 * windows=...)} bit-for-bit:</p>
 *
 * <pre>
 * pi_branch(window) = (1 / span(window))
 *                   × Σ_{tree ∩ window} span(tree ∩ window)
 *                   × Σ_{branch} t_b × 2 × k_b × (n - k_b) / (n × (n - 1))
 * </pre>
 *
 * <p>where {@code n} = focal-sample count, {@code k_b} = focal-sample
 * descendants of branch {@code b}'s child, {@code t_b} =
 * {@code parent.time - child.time}.</p>
 */
public class ArgBranchDiversityProcedure {

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.arg.branch_diversity", mode = Mode.READ)
    @Description("Branch-mode diversity from an ARG. mode='pi' (v1 default); "
            + "options.windows for per-window output.")
    @SuppressWarnings("unchecked")
    public Stream<BranchDiversityResult> branchDiversity(
            @Name("runId") String runId,
            @Name("sampleIds") List<String> sampleIds,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        if (options == null) options = new HashMap<>();
        String mode = (String) options.getOrDefault("mode", "pi");
        Object windowsObj = options.get("windows");

        if (!"pi".equals(mode)) {
            throw new IllegalArgumentException(
                    "Only mode='pi' is supported in v1; got: " + mode);
        }

        ARG arg = ARGTraversal.load(tx, runId, 0L, Long.MAX_VALUE);
        if (arg.nNodes == 0 || sampleIds == null || sampleIds.isEmpty()) {
            return Stream.empty();
        }

        int[] focal = resolveFocal(tx, runId, sampleIds, arg);
        if (focal.length < 2) return Stream.empty();
        final int n = focal.length;
        final double pairs = (double) n * (n - 1);

        long[] windowEdges;
        if (windowsObj instanceof List) {
            List<?> ws = (List<?>) windowsObj;
            windowEdges = new long[ws.size()];
            for (int i = 0; i < ws.size(); i++) {
                windowEdges[i] = ((Number) ws.get(i)).longValue();
            }
        } else {
            windowEdges = new long[]{0L, arg.sequenceLength};
        }

        List<BranchDiversityResult> rows = new ArrayList<>(windowEdges.length - 1);
        for (int w = 0; w < windowEdges.length - 1; w++) {
            long lo = windowEdges[w];
            long hi = windowEdges[w + 1];
            if (hi <= lo) {
                rows.add(new BranchDiversityResult(lo, hi, 0.0, n, mode, runId));
                continue;
            }
            double pi = piBranchOver(arg, focal, lo, hi, pairs);
            rows.add(new BranchDiversityResult(lo, hi, pi, n, mode, runId));
        }
        return rows.stream();
    }

    private static double piBranchOver(ARG arg, int[] focal, long lo, long hi,
                                        double pairs) {
        long[] bp = arg.breakpointsClipped(lo, hi);
        if (bp.length < 2) return 0.0;
        double accum = 0.0;
        for (int k = 0; k < bp.length - 1; k++) {
            long b0 = bp[k];
            long b1 = bp[k + 1];
            if (b1 <= b0) continue;
            int[] parent = ArgTraversalUtils.marginalTree(arg, b0);
            int[] desc = ArgTraversalUtils.descendantCounts(arg, parent, focal);
            double treeContribution = 0.0;
            for (int e = 0; e < arg.nEdges; e++) {
                if (arg.edgeStart[e] > b0 || arg.edgeEnd[e] <= b0) continue;
                int p = arg.edgeParent[e];
                int c = arg.edgeChild[e];
                int kb = desc[c];
                int nk = focal.length - kb;
                if (kb == 0 || nk == 0) continue;
                double tb = arg.time[p] - arg.time[c];
                treeContribution += tb * 2.0 * kb * nk / pairs;
            }
            long span = b1 - b0;
            accum += treeContribution * span;
        }
        // Normalise by window length to match tskit's per-bp diversity.
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
}
