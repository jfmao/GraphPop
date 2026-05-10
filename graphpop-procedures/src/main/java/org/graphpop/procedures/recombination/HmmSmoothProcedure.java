package org.graphpop.procedures.recombination;

import org.graphpop.procedures.selection.SelectionUtils;
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
 * Pyrho-style HMM smoothing of per-window ρ (M13.B).
 *
 * <p>For each sliding window, compute the M11 moment estimator
 * (mean r² over in-window pairs → Hudson 1985 bisection); then
 * run a discrete-state HMM over the resulting log10(ρ̂)
 * sequence to share strength across adjacent windows. Output:
 * raw and smoothed ρ per window plus the Viterbi state.</p>
 *
 * <pre>
 * CALL graphpop.recombination.hmm_smooth($sampleIds,
 *         {window_size: 10000, step: 5000, n_states: 20,
 *          state_log_lo: -10, state_log_hi: -4,
 *          emission_sd: 0.5, switch_rate: 0.1})
 *   YIELD start, end, rho_per_bp, rho_smoothed, hmm_state,
 *         n_variant_pairs, runId
 * </pre>
 */
public class HmmSmoothProcedure {

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.recombination.hmm_smooth", mode = Mode.READ)
    @Description("Pyrho-style HMM smoothing of per-window ρ. "
            + "Discrete-state HMM over log10(ρ̂) with Gaussian "
            + "emissions and banded-walk transition prior.")
    public Stream<HmmSmoothResult> hmmSmooth(
            @Name("sampleIds") List<String> sampleIds,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        if (options == null) options = new HashMap<>();
        if (sampleIds == null || sampleIds.isEmpty()) return Stream.empty();
        long windowSize = getLong(options, "window_size", 10_000L);
        long step = getLong(options, "step", windowSize);
        double minMaf = getDouble(options, "min_maf", 0.05);
        long maxPairDistance = getLong(options, "max_pair_distance", 5_000L);
        int nStates = (int) getLong(options, "n_states", 20L);
        double stateLogLo = getDouble(options, "state_log_lo", -10.0);
        double stateLogHi = getDouble(options, "state_log_hi", -4.0);
        double emissionSd = getDouble(options, "emission_sd", 0.5);
        double switchRate = getDouble(options, "switch_rate", 0.1);

        if (windowSize <= 0 || step <= 0 || maxPairDistance <= 0) {
            throw new IllegalArgumentException(
                    "window_size, step, max_pair_distance must be positive");
        }

        long seqEnd = maxVariantPosition(tx);
        if (seqEnd <= 0) return Stream.empty();
        long[][] windows = SelectionUtils.slidingWindows(
                seqEnd + 1, windowSize, step);
        if (windows.length == 0) return Stream.empty();
        String runId = peekRunId(tx);

        // Pass 1: per-window raw ρ estimates via the moment estimator.
        double[] log10Rho = new double[windows.length];
        double[] rawRho = new double[windows.length];
        long[] nPairs = new long[windows.length];
        for (int w = 0; w < windows.length; w++) {
            long lo = windows[w][0];
            long hi = windows[w][1];
            List<LdPairLoader.Pair> pairs = LdPairLoader.loadPairs(
                    tx, sampleIds, lo, hi, minMaf, maxPairDistance);
            nPairs[w] = pairs.size();
            if (pairs.isEmpty()) {
                rawRho[w] = 0.0;
                log10Rho[w] = Double.NaN;
                continue;
            }
            double sumR2 = 0.0;
            double sumD = 0.0;
            for (LdPairLoader.Pair p : pairs) {
                sumR2 += p.r2;
                sumD += p.distance;
            }
            double meanR2 = sumR2 / pairs.size();
            double meanD = sumD / pairs.size();
            double rho = HudsonRecombination.solveRhoFromMeanR2(meanR2, meanD);
            if (Double.isNaN(rho) || rho <= 0.0) {
                rawRho[w] = 0.0;
                log10Rho[w] = Double.NaN;
            } else {
                rawRho[w] = rho;
                log10Rho[w] = Math.log10(rho);
            }
        }

        // Pass 2: HMM smoothing on log10 scale.
        HmmSmoothComputer.Result smooth = HmmSmoothComputer.smooth(
                log10Rho, nStates, stateLogLo, stateLogHi,
                emissionSd, switchRate);

        List<HmmSmoothResult> rows = new ArrayList<>(windows.length);
        for (int w = 0; w < windows.length; w++) {
            double smoothed = Math.pow(10.0, smooth.smoothedLog10Rho[w]);
            rows.add(new HmmSmoothResult(
                    windows[w][0], windows[w][1],
                    rawRho[w], smoothed,
                    (long) smooth.viterbiStates[w],
                    nPairs[w], runId));
        }
        return rows.stream();
    }

    private static long maxVariantPosition(Transaction tx) {
        Result r = tx.execute(
                "MATCH (v:Variant) RETURN max(v.position) AS p");
        try {
            if (r.hasNext()) {
                Object p = r.next().get("p");
                if (p instanceof Number) return ((Number) p).longValue();
            }
        } finally {
            r.close();
        }
        return 0L;
    }

    private static String peekRunId(Transaction tx) {
        Result r = tx.execute(
                "MATCH ()-[m:MUTATED_ON]->() RETURN m.runId AS rid LIMIT 1");
        try {
            if (r.hasNext()) {
                Object rid = r.next().get("rid");
                if (rid != null) return (String) rid;
            }
        } finally {
            r.close();
        }
        return "";
    }

    private static long getLong(Map<String, Object> opts, String key, long def) {
        Object v = opts.get(key);
        return (v instanceof Number) ? ((Number) v).longValue() : def;
    }

    private static double getDouble(Map<String, Object> opts, String key,
                                    double def) {
        Object v = opts.get(key);
        return (v instanceof Number) ? ((Number) v).doubleValue() : def;
    }
}
