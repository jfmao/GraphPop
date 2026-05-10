package org.graphpop.procedures.demography;

import org.graphpop.procedures.arg.ArgCoalescenceRateComputer;
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
 * Demographic-inference procedure (M7).
 *
 * <p>Closed-form per-bin Ne(t) inversion of M6's coalescence-rate
 * estimator:</p>
 *
 * <pre>
 * Ne(bin) = 1 / (ploidy × rate(bin))
 * SE(Ne)  = (1 / (ploidy · rate²)) · √(rate / lineage_pair_time)
 *           [delta method on Poisson event count]
 * </pre>
 *
 * <p>When {@code rate(bin) == 0} the bin yields no events; the row
 * carries {@code ne = +Infinity}, {@code ne_se = NaN}, and
 * {@code flag = "no_events"}.</p>
 *
 * <pre>
 * CALL graphpop.demography.ne_trajectory($runId, $sampleIds,
 *     {time_bins: [0, 0.25, 0.5, 1, 2, 1e9], ploidy: 2})
 *   YIELD time_lo, time_hi, n_coalescent_events, lineage_pair_time,
 *         rate, ne, ne_se, flag, runId
 * </pre>
 *
 * <p>For convenience, {@code options.population} resolves
 * {@code sampleIds} from {@code :Sample.population}; the explicit
 * {@code sampleIds} list overrides it when both are passed.</p>
 */
public class NeTrajectoryProcedure {

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.demography.ne_trajectory", mode = Mode.READ)
    @Description("Closed-form Ne(t) trajectory by inverting "
            + "graphpop.arg.coalescence_rate per bin. Ne = 1/(ploidy*rate); "
            + "SE via delta method.")
    @SuppressWarnings("unchecked")
    public Stream<NeTrajectoryResult> neTrajectory(
            @Name("runId") String runId,
            @Name(value = "sampleIds", defaultValue = "[]") List<String> sampleIds,
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

        long ploidy = 2L;
        Object pObj = options.get("ploidy");
        if (pObj instanceof Number) {
            ploidy = ((Number) pObj).longValue();
            if (ploidy < 1L) {
                throw new IllegalArgumentException("ploidy must be >= 1");
            }
        }
        Object popObj = options.get("population");

        ARG arg = ARGTraversal.load(tx, runId, 0L, Long.MAX_VALUE);
        if (arg.nNodes == 0) return Stream.empty();

        // Resolve focal sample-id set: explicit > population > nothing.
        List<String> focalIds = (sampleIds != null && !sampleIds.isEmpty())
                ? sampleIds
                : resolvePopulationSamples(tx, runId, popObj);
        if (focalIds == null || focalIds.isEmpty()) return Stream.empty();

        int[] focal = resolveFocal(tx, runId, focalIds, arg);
        if (focal.length < 2) return Stream.empty();

        List<ArgCoalescenceRateComputer.BinResult> binResults =
                ArgCoalescenceRateComputer.compute(arg, focal, bins);

        List<NeTrajectoryResult> rows = new ArrayList<>(binResults.size());
        double pl = ploidy;
        for (ArgCoalescenceRateComputer.BinResult b : binResults) {
            double ne;
            double neSe;
            String flag;
            if (b.rate <= 0.0 || b.pairTime <= 0.0) {
                ne = Double.POSITIVE_INFINITY;
                neSe = Double.NaN;
                flag = "no_events";
            } else {
                ne = 1.0 / (pl * b.rate);
                // delta method:  SE(Ne) = (1 / (pl * rate^2)) * sqrt(rate / T)
                neSe = (1.0 / (pl * b.rate * b.rate))
                       * Math.sqrt(b.rate / b.pairTime);
                flag = "ok";
            }
            rows.add(new NeTrajectoryResult(
                    b.timeLo, b.timeHi,
                    b.events, b.pairTime, b.rate,
                    ne, neSe, flag, runId));
        }
        return rows.stream();
    }

    private static List<String> resolvePopulationSamples(Transaction tx,
                                                          String runId,
                                                          Object popObj) {
        if (!(popObj instanceof String)) return List.of();
        String pop = (String) popObj;
        // Pull samples in this run that have :Sample.population = pop.
        Result r = tx.execute(
                "MATCH (:TreeNode {runId: $runId})-[:REPRESENTS]->"
              + "(s:Sample {population: $pop}) "
              + "RETURN DISTINCT s.sampleId AS sid",
                Map.of("runId", runId, "pop", pop));
        List<String> out = new ArrayList<>();
        try {
            while (r.hasNext()) {
                out.add((String) r.next().get("sid"));
            }
        } finally {
            r.close();
        }
        return out;
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
