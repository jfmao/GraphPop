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
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

/**
 * Stratified LD-decay recombination-map (M13.C).
 *
 * <pre>
 * CALL graphpop.recombination.stratified_ld_decay($sampleIds,
 *         {stratify_by: 'population',
 *          window_size: 10000, step: 5000,
 *          max_pair_distance: 5000, min_maf: 0.05})
 *   YIELD start, end, stratum, n_variant_pairs, mean_r2,
 *         mean_pair_distance, rho_per_bp, n_samples,
 *         stratify_by, runId
 * </pre>
 *
 * <p>{@code stratify_by} ∈ {@code "population"} or {@code "sex"}.
 * Splits the focal sample set by the named {@code :Sample}
 * property and runs the M11 Hudson-Kaplan moment estimator
 * independently per stratum. Useful for sex-averaged versus
 * sex-specific ρ-maps and for per-population ρ in admixed cohorts.</p>
 *
 * <p>The {@code stratify_by = 'ancestry_block'} mode that walks
 * M4.B {@code :HAS_ANCESTRY} edges is deferred to a v2.</p>
 */
public class StratifiedLdDecayProcedure {

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.recombination.stratified_ld_decay",
               mode = Mode.READ)
    @Description("Per-(window, stratum) Hudson-Kaplan moment estimator "
            + "on :CARRIES r² over a focal sample set split by "
            + ":Sample.population or :Sample.sex.")
    public Stream<StratifiedLdDecayResult> stratifiedLdDecay(
            @Name("sampleIds") List<String> sampleIds,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        if (options == null) options = new HashMap<>();
        if (sampleIds == null || sampleIds.isEmpty()) return Stream.empty();
        String stratifyBy = ((String) options.getOrDefault(
                "stratify_by", "population")).toLowerCase();
        if (!"population".equals(stratifyBy) && !"sex".equals(stratifyBy)) {
            throw new IllegalArgumentException(
                    "v1 supports stratify_by ∈ {population, sex}; got: "
                  + stratifyBy);
        }
        long windowSize = getLong(options, "window_size", 10_000L);
        long step = getLong(options, "step", windowSize);
        double minMaf = getDouble(options, "min_maf", 0.05);
        long maxPairDistance = getLong(options, "max_pair_distance", 5_000L);

        if (windowSize <= 0 || step <= 0 || maxPairDistance <= 0) {
            throw new IllegalArgumentException(
                    "window_size, step, max_pair_distance must be positive");
        }

        Map<String, List<String>> stratumToSamples =
                resolveStrata(tx, sampleIds, stratifyBy);
        if (stratumToSamples.isEmpty()) return Stream.empty();

        long seqEnd = maxVariantPosition(tx);
        if (seqEnd <= 0) return Stream.empty();
        long[][] windows = SelectionUtils.slidingWindows(
                seqEnd + 1, windowSize, step);
        if (windows.length == 0) return Stream.empty();

        String runId = peekRunId(tx);

        List<StratifiedLdDecayResult> rows = new ArrayList<>();
        for (Map.Entry<String, List<String>> e : stratumToSamples.entrySet()) {
            String stratum = e.getKey();
            List<String> stratumSids = e.getValue();
            if (stratumSids.size() < 2) continue;  // need ≥ 2 to estimate r²
            for (long[] win : windows) {
                long lo = win[0];
                long hi = win[1];
                List<LdPairLoader.Pair> pairs = LdPairLoader.loadPairs(
                        tx, stratumSids, lo, hi, minMaf, maxPairDistance);
                double sumR2 = 0.0;
                double sumD = 0.0;
                for (LdPairLoader.Pair p : pairs) {
                    sumR2 += p.r2;
                    sumD += p.distance;
                }
                double meanR2 = pairs.isEmpty()
                        ? Double.NaN : sumR2 / pairs.size();
                double meanD = pairs.isEmpty()
                        ? Double.NaN : sumD / pairs.size();
                double rho = 0.0;
                if (!pairs.isEmpty() && meanD > 0.0) {
                    rho = HudsonRecombination.solveRhoFromMeanR2(meanR2, meanD);
                    if (Double.isNaN(rho)) rho = 0.0;
                }
                rows.add(new StratifiedLdDecayResult(
                        lo, hi, stratum, pairs.size(),
                        meanR2, meanD, rho, stratumSids.size(),
                        stratifyBy, runId));
            }
        }
        return rows.stream();
    }

    private static Map<String, List<String>> resolveStrata(
            Transaction tx, List<String> sampleIds, String stratifyBy) {
        String prop = "population".equals(stratifyBy) ? "population" : "sex";
        Result r = tx.execute(
                "MATCH (s:Sample) WHERE s.sampleId IN $sids "
              + "RETURN s.sampleId AS sid, s." + prop + " AS stratum",
                Map.of("sids", sampleIds));
        // Linked to preserve a stable output order.
        Map<String, List<String>> out = new LinkedHashMap<>();
        try {
            while (r.hasNext()) {
                Map<String, Object> row = r.next();
                String sid = (String) row.get("sid");
                Object stratumObj = row.get("stratum");
                if (stratumObj == null) continue;
                String stratum = String.valueOf(stratumObj);
                out.computeIfAbsent(stratum, k -> new ArrayList<>()).add(sid);
            }
        } finally {
            r.close();
        }
        return out;
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
