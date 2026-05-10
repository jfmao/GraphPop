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
import java.util.TreeSet;
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
 * focal-sample descendant. The reference Python implementation lives
 * at {@code build_egrm_fixture.py:coalescence_rate_reference} and is
 * mirrored here bit-for-bit.</p>
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
        if (bins.length < 2) return Stream.empty();
        for (int i = 1; i < bins.length; i++) {
            if (bins[i] <= bins[i - 1]) {
                throw new IllegalArgumentException(
                        "time_bins must be strictly increasing");
            }
        }

        ARG arg = ARGTraversal.load(tx, runId, 0L, Long.MAX_VALUE);
        if (arg.nNodes == 0 || sampleIds == null || sampleIds.isEmpty()) {
            return Stream.empty();
        }
        int[] focal = resolveFocal(tx, runId, sampleIds, arg);
        if (focal.length < 2) return Stream.empty();

        int nBins = bins.length - 1;
        double[] events = new double[nBins];
        double[] pairTime = new double[nBins];
        double seqLen = arg.sequenceLength;

        long[] bp = arg.breakpointsClipped(0L, arg.sequenceLength);
        for (int k = 0; k < bp.length - 1; k++) {
            long b0 = bp[k];
            long b1 = bp[k + 1];
            if (b1 <= b0) continue;
            double s = (b1 - b0) / seqLen;

            int[] parent = ArgTraversalUtils.marginalTree(arg, b0);
            int[] desc = ArgTraversalUtils.descendantCounts(arg, parent, focal);

            // Build children lists for this marginal tree.
            Map<Integer, List<Integer>> children = new HashMap<>();
            for (int i = 0; i < arg.nNodes; i++) {
                if (parent[i] != -1) {
                    children.computeIfAbsent(parent[i], x -> new ArrayList<>())
                            .add(i);
                }
            }

            // Coalescence events: at each internal node, the focal-pair
            // events here = C(total, 2) - Σ C(c_i, 2).
            for (int u = 0; u < arg.nNodes; u++) {
                List<Integer> cs = children.get(u);
                if (cs == null || cs.size() < 2) continue;
                int total = 0;
                int sumChooseTwo = 0;
                for (int c : cs) {
                    int ck = desc[c];
                    total += ck;
                    sumChooseTwo += ck * (ck - 1) / 2;
                }
                if (total < 2) continue;
                int cross = (total * (total - 1) / 2) - sumChooseTwo;
                if (cross == 0) continue;
                int bi = binIndex(bins, arg.time[u]);
                if (bi >= 0) events[bi] += cross * s;
            }

            // Lineage-pair time: piecewise-constant k_t over node-time grid.
            TreeSet<Double> timeSet = new TreeSet<>();
            for (int i = 0; i < arg.nNodes; i++) {
                if (parent[i] != -1 || children.containsKey(i)) {
                    timeSet.add(arg.time[i]);
                }
            }
            Double[] timeArr = timeSet.toArray(new Double[0]);
            for (int t = 0; t < timeArr.length - 1; t++) {
                double tLo = timeArr[t];
                double tHi = timeArr[t + 1];
                int kT = 0;
                for (int u = 0; u < arg.nNodes; u++) {
                    int p = parent[u];
                    if (p == -1) continue;
                    double ct = arg.time[u];
                    double pt = arg.time[p];
                    if (ct <= tLo && pt >= tHi && desc[u] > 0) {
                        kT++;
                    }
                }
                double pairs = kT * (kT - 1) / 2.0;
                if (pairs == 0) continue;
                for (int bi = 0; bi < nBins; bi++) {
                    double lo = Math.max(bins[bi], tLo);
                    double hi = Math.min(bins[bi + 1], tHi);
                    if (hi > lo) pairTime[bi] += pairs * (hi - lo) * s;
                }
            }
        }

        List<CoalescenceRateResult> rows = new ArrayList<>(nBins);
        for (int bi = 0; bi < nBins; bi++) {
            double rate = (pairTime[bi] > 0) ? events[bi] / pairTime[bi] : 0.0;
            rows.add(new CoalescenceRateResult(
                    bins[bi], bins[bi + 1],
                    events[bi], pairTime[bi], rate, runId));
        }
        return rows.stream();
    }

    private static int binIndex(double[] bins, double t) {
        for (int i = 0; i < bins.length - 1; i++) {
            if (bins[i] <= t && t < bins[i + 1]) return i;
        }
        return -1;
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
