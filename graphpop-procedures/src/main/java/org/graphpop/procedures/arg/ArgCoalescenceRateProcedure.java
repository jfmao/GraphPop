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
 * Per-bin coalescence rate from an ARG (M6 part C).
 *
 * <p>Speidel/tsdate-style estimator:</p>
 *
 * <pre>
 * T(bin) = ∫ k_t × (k_t - 1) / 2 dt  summed over marginal trees,
 *           weighted by tree-span / sequence_length
 * C(bin) = focal-sample-lineage coalescent events with parent_time
 *           ∈ bin, weighted by tree-span / sequence_length
 * rate(bin) = C(bin) / T(bin)
 * </pre>
 *
 * <p>k_t is the count of branches alive at time t with at least one
 * focal-sample descendant. Aggregation lives in
 * {@link ArgCoalescenceRateComputer} and is shared with M7's
 * demographic-inference procedure.</p>
 *
 * <pre>
 * CALL graphpop.arg.coalescence_rate($runId, $sampleIds,
 *                                     {time_bins: [0,100,1000,10000,100000]})
 *   YIELD time_lo, time_hi, n_coalescent_events, lineage_pair_time,
 *         rate, runId
 * </pre>
 */
public class ArgCoalescenceRateProcedure {

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.arg.coalescence_rate", mode = Mode.READ)
    @Description("Per-bin coalescence rate from an ARG. Mirrors the "
            + "Speidel/tsdate estimator: rate(bin) = events / "
            + "lineage_pair_time, both span-weighted across marginal trees.")
    @SuppressWarnings("unchecked")
    public Stream<CoalescenceRateResult> coalescenceRate(
            @Name("runId") String runId,
            @Name("sampleIds") List<String> sampleIds,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        if (options == null) options = new HashMap<>();
        Object binsObj = options.get("time_bins");
        if (!(binsObj instanceof List)) {
            throw new IllegalArgumentException(
                    "options.time_bins is required and must be a list of "
                  + "monotonically increasing numbers");
        }
        List<?> binsList = (List<?>) binsObj;
        double[] bins = new double[binsList.size()];
        for (int i = 0; i < bins.length; i++) {
            bins[i] = ((Number) binsList.get(i)).doubleValue();
        }
        ArgCoalescenceRateComputer.validateBins(bins);

        ARG arg = ARGTraversal.load(tx, runId, 0L, Long.MAX_VALUE);
        if (arg.nNodes == 0 || sampleIds == null || sampleIds.isEmpty()) {
            return Stream.empty();
        }
        int[] focal = resolveFocal(tx, runId, sampleIds, arg);
        if (focal.length < 2) return Stream.empty();

        List<ArgCoalescenceRateComputer.BinResult> bins2 =
                ArgCoalescenceRateComputer.compute(arg, focal, bins);
        List<CoalescenceRateResult> rows = new ArrayList<>(bins2.size());
        for (ArgCoalescenceRateComputer.BinResult b : bins2) {
            rows.add(new CoalescenceRateResult(
                    b.timeLo, b.timeHi, b.events, b.pairTime, b.rate, runId));
        }
        return rows.stream();
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
