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
 * Recombination-hotspot detection with Benjamini-Hochberg FDR (M13.D).
 *
 * <pre>
 * CALL graphpop.recombination.hotspots($sampleIds,
 *         {window_size: 10000, step: 5000, fdr_q: 0.05,
 *          method: 'ld_decay'})
 *   YIELD start, end, rho_per_bp, z_score, p_value, adj_p_value,
 *         is_hotspot, n_variant_pairs, method, runId
 * </pre>
 *
 * <p>Pipeline (v1):</p>
 *
 * <ol>
 *   <li>Compute per-window ρ_per_bp via M11's
 *       {@code LdDecayProcedure} logic (the only method exposed in
 *       v1 here; {@code arg_breakpoints} reserved for v2).</li>
 *   <li>Genome-wide null = {@code (μ, σ)} of {@code log10(ρ_per_bp)}
 *       across all non-zero-ρ windows (via Welford).</li>
 *   <li>Per-window z-score → one-sided p-value (right tail =
 *       hotspot) via the standard normal CDF complement.</li>
 *   <li>Benjamini-Hochberg FDR adjustment at {@code fdr_q}.</li>
 *   <li>Flag windows with {@code adj_p_value < fdr_q} as hotspots.</li>
 * </ol>
 *
 * <p>Output is ordered by genomic position (not adjusted p-value).</p>
 */
public class HotspotsProcedure {

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.recombination.hotspots", mode = Mode.READ)
    @Description("Recombination-hotspot detection with Benjamini-Hochberg "
            + "FDR over per-window ρ (LD-decay moment estimator).")
    public Stream<HotspotResult> hotspots(
            @Name("sampleIds") List<String> sampleIds,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        if (options == null) options = new HashMap<>();
        if (sampleIds == null || sampleIds.isEmpty()) return Stream.empty();
        String method = ((String) options.getOrDefault("method", "ld_decay"))
                .toLowerCase();
        if (!"ld_decay".equals(method)) {
            throw new IllegalArgumentException(
                    "v1 only supports method='ld_decay'; got: " + method);
        }
        long windowSize = getLong(options, "window_size", 10_000L);
        long step = getLong(options, "step", windowSize);
        double minMaf = getDouble(options, "min_maf", 0.05);
        long maxPairDistance = getLong(options, "max_pair_distance", 5_000L);
        double fdrQ = getDouble(options, "fdr_q", 0.05);
        if (windowSize <= 0 || step <= 0 || maxPairDistance <= 0
                || fdrQ <= 0.0 || fdrQ >= 1.0) {
            throw new IllegalArgumentException(
                    "window_size / step / max_pair_distance must be positive "
                  + "and fdr_q ∈ (0, 1)");
        }

        long seqEnd = maxVariantPosition(tx);
        if (seqEnd <= 0) return Stream.empty();
        long[][] windows = SelectionUtils.slidingWindows(
                seqEnd + 1, windowSize, step);
        if (windows.length == 0) return Stream.empty();
        String runId = peekRunId(tx);

        // Pass 1: per-window ρ.
        double[] rho = new double[windows.length];
        long[] nPairs = new long[windows.length];
        for (int w = 0; w < windows.length; w++) {
            List<LdPairLoader.Pair> pairs = LdPairLoader.loadPairs(
                    tx, sampleIds, windows[w][0], windows[w][1],
                    minMaf, maxPairDistance);
            nPairs[w] = pairs.size();
            if (pairs.isEmpty()) {
                rho[w] = 0.0;
                continue;
            }
            double sumR2 = 0.0;
            double sumD = 0.0;
            for (LdPairLoader.Pair p : pairs) {
                sumR2 += p.r2;
                sumD += p.distance;
            }
            double mr2 = sumR2 / pairs.size();
            double md = sumD / pairs.size();
            double rhoHat = (md > 0.0)
                    ? HudsonRecombination.solveRhoFromMeanR2(mr2, md)
                    : 0.0;
            rho[w] = Double.isNaN(rhoHat) ? 0.0 : rhoHat;
        }

        // Pass 2: genome-wide null over log10(ρ) for ρ > 0.
        SelectionUtils.Welford nullStats = new SelectionUtils.Welford();
        for (double r : rho) {
            if (r > 0.0) nullStats.add(Math.log10(r));
        }
        double mean = nullStats.mean();
        double sd = nullStats.sd();

        // Pass 3: z-score + one-sided p-values.
        double[] pValues = new double[windows.length];
        double[] zScores = new double[windows.length];
        for (int w = 0; w < windows.length; w++) {
            if (rho[w] <= 0.0 || Double.isNaN(sd) || sd <= 0.0) {
                zScores[w] = 0.0;
                pValues[w] = 1.0;
            } else {
                double z = (Math.log10(rho[w]) - mean) / sd;
                zScores[w] = z;
                pValues[w] = oneSidedUpperTailNormal(z);
            }
        }

        // Pass 4: Benjamini-Hochberg adjusted p-values.
        double[] adjP = benjaminiHochberg(pValues);

        // Pass 5: emit rows in genomic order.
        List<HotspotResult> rows = new ArrayList<>(windows.length);
        for (int w = 0; w < windows.length; w++) {
            boolean isHotspot = adjP[w] < fdrQ;
            rows.add(new HotspotResult(
                    windows[w][0], windows[w][1],
                    rho[w], zScores[w], pValues[w], adjP[w],
                    isHotspot, nPairs[w], method, runId));
        }
        return rows.stream();
    }

    /**
     * Benjamini-Hochberg adjusted p-values. Returns an array where
     * {@code adj[i]} is the BH-corrected p-value for the original
     * (un-sorted) input position {@code i}. NaN and out-of-range
     * inputs survive as 1.0.
     */
    static double[] benjaminiHochberg(double[] pValues) {
        int n = pValues.length;
        Integer[] order = new Integer[n];
        for (int i = 0; i < n; i++) order[i] = i;
        java.util.Arrays.sort(order, (a, b) -> Double.compare(
                safe(pValues[a]), safe(pValues[b])));
        double[] adj = new double[n];
        double prev = 1.0;
        for (int rank = n - 1; rank >= 0; rank--) {
            int idx = order[rank];
            double p = safe(pValues[idx]);
            double bh = p * n / (rank + 1.0);
            bh = Math.min(bh, prev);
            adj[idx] = Math.min(1.0, bh);
            prev = adj[idx];
        }
        return adj;
    }

    private static double safe(double p) {
        if (Double.isNaN(p) || p > 1.0) return 1.0;
        if (p < 0.0) return 0.0;
        return p;
    }

    /** One-sided upper-tail p-value for the standard normal at z. */
    static double oneSidedUpperTailNormal(double z) {
        // Φ(z) via the Abramowitz & Stegun erf approximation.
        return 0.5 * (1.0 - erf(z / Math.sqrt(2.0)));
    }

    /** Abramowitz & Stegun 7.1.26: max abs err ≈ 1.5e-7. */
    private static double erf(double x) {
        double sign = (x < 0) ? -1.0 : 1.0;
        x = Math.abs(x);
        double a1 =  0.254829592;
        double a2 = -0.284496736;
        double a3 =  1.421413741;
        double a4 = -1.453152027;
        double a5 =  1.061405429;
        double p =  0.3275911;
        double t = 1.0 / (1.0 + p * x);
        double y = 1.0
                - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1)
                        * t * Math.exp(-x * x);
        return sign * y;
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
