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
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

/**
 * Per-window LD-decay table + Hudson-Kaplan moment estimator
 * (M11 part B).
 *
 * <pre>
 * CALL graphpop.recombination.ld_decay($sampleIds,
 *         {window_size: 10000, step: 5000, min_maf: 0.05,
 *          max_pair_distance: 5000})
 *   YIELD start, end, n_variant_pairs, mean_r2,
 *         mean_pair_distance, rho_per_bp, n_samples,
 *         method, runId
 * </pre>
 *
 * <p>Per window:</p>
 *
 * <ol>
 *   <li>Pull every {@code :Variant} in the window with derived-
 *       allele frequency ≥ {@code min_maf} in the focal sample set.</li>
 *   <li>For every pair of variants at physical distance
 *       {@code d ≤ max_pair_distance}, compute r² from the 2×2
 *       haplotype counts (Hill 1968 closed form).</li>
 *   <li>Aggregate {@code mean(r²)} and {@code mean(d)} across
 *       pairs; invert via {@link HudsonRecombination#solveRhoFromMeanR2}
 *       at the mean distance to get {@code rho_per_bp}.</li>
 * </ol>
 *
 * <p>{@code method = "hudson_moment"} in v1. Slot reserved for
 * LDhat-style composite likelihood in a v2.</p>
 */
public class LdDecayProcedure {

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.recombination.ld_decay", mode = Mode.READ)
    @Description("Per-window LD-decay + Hudson-Kaplan moment "
            + "estimator. Reports mean r², mean pair distance, and "
            + "ρ-per-bp inverted from Hudson 1985's E[r²|n,C].")
    @SuppressWarnings("unchecked")
    public Stream<LdDecayResult> ldDecay(
            @Name("sampleIds") List<String> sampleIds,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        if (options == null) options = new HashMap<>();
        if (sampleIds == null || sampleIds.isEmpty()) return Stream.empty();
        long windowSize = getLong(options, "window_size", 10_000L);
        long step = getLong(options, "step", windowSize);
        double minMaf = getDouble(options, "min_maf", 0.05);
        long maxPairDistance = getLong(options, "max_pair_distance", 5_000L);
        if (windowSize <= 0 || step <= 0) {
            throw new IllegalArgumentException(
                    "window_size and step must be positive");
        }
        if (maxPairDistance <= 0) {
            throw new IllegalArgumentException(
                    "max_pair_distance must be positive");
        }

        // Find genome span (max position) and the focal sample count.
        Set<String> sampleSet = new HashSet<>(sampleIds);
        int nSamples = sampleSet.size();
        if (nSamples > 63) {
            throw new IllegalArgumentException(
                    "ld_decay v1 requires ≤ 63 focal samples (bitmask "
                  + "representation). Got " + nSamples + ".");
        }
        long seqEnd = maxVariantPosition(tx);
        if (seqEnd <= 0) return Stream.empty();
        long[][] windows = SelectionUtils.slidingWindows(
                seqEnd + 1, windowSize, step);
        if (windows.length == 0) return Stream.empty();

        // Resolve runId opportunistically from the first :MUTATED_ON
        // edge — used only for the runId column.
        String runId = peekRunId(tx);

        List<LdDecayResult> rows = new ArrayList<>(windows.length);
        for (long[] win : windows) {
            long lo = win[0];
            long hi = win[1];
            // Variants in window with their carrier sets.
            List<long[]> variants = loadVariantsInWindow(
                    tx, sampleSet, lo, hi, minMaf);
            // Sort by position so we can early-exit on distance.
            variants.sort((a, b) -> Long.compare(a[0], b[0]));

            double sumR2 = 0.0;
            double sumD = 0.0;
            long nPairs = 0L;
            for (int i = 0; i < variants.size(); i++) {
                long posI = variants.get(i)[0];
                long carriersI = variants.get(i)[1];
                int kI = Long.bitCount(carriersI);
                for (int j = i + 1; j < variants.size(); j++) {
                    long posJ = variants.get(j)[0];
                    long d = posJ - posI;
                    if (d > maxPairDistance) break;
                    long carriersJ = variants.get(j)[1];
                    int kJ = Long.bitCount(carriersJ);
                    int kAB = Long.bitCount(carriersI & carriersJ);
                    double r2 = rSquared(kI, kJ, kAB, nSamples);
                    if (Double.isFinite(r2)) {
                        sumR2 += r2;
                        sumD += d;
                        nPairs++;
                    }
                }
            }

            double meanR2 = (nPairs > 0) ? sumR2 / nPairs : Double.NaN;
            double meanD = (nPairs > 0) ? sumD / nPairs : Double.NaN;
            double rho = 0.0;
            if (nPairs > 0 && meanD > 0.0) {
                rho = HudsonRecombination.solveRhoFromMeanR2(meanR2, meanD);
                if (Double.isNaN(rho)) rho = 0.0;
            }
            rows.add(new LdDecayResult(
                    lo, hi, nPairs, meanR2, meanD, rho, nSamples,
                    "hudson_moment", runId));
        }
        return rows.stream();
    }

    /**
     * Returns variants in the window with their carrier bitmask
     * over the focal sample set. Each long-array row is
     * {@code [position, carrier_bitmask]}.
     */
    private static List<long[]> loadVariantsInWindow(
            Transaction tx, Set<String> focalSet,
            long lo, long hi, double minMaf) {
        if (focalSet.size() > 63) {
            throw new IllegalArgumentException(
                    "ld_decay v1 requires ≤ 63 focal samples (bitmask "
                  + "representation). Got " + focalSet.size() + ".");
        }
        // Build a stable focal-sample index (sample_id → bit).
        List<String> ordered = new ArrayList<>(focalSet);
        ordered.sort(null);
        Map<String, Integer> bit = new HashMap<>();
        for (int i = 0; i < ordered.size(); i++) bit.put(ordered.get(i), i);

        // Pull CARRIES edges within the window's variants.
        Result r = tx.execute(
                "MATCH (s:Sample)-[c:CARRIES]->(v:Variant) "
              + "WHERE v.position >= $lo AND v.position < $hi "
              + "  AND s.sampleId IN $sids "
              + "RETURN v.variantId AS vid, v.position AS pos, "
              + "s.sampleId AS sid",
                Map.of("lo", lo, "hi", hi, "sids", ordered));
        Map<String, long[]> byVariant = new HashMap<>();
        try {
            while (r.hasNext()) {
                Map<String, Object> row = r.next();
                String vid = (String) row.get("vid");
                long pos = ((Number) row.get("pos")).longValue();
                String sid = (String) row.get("sid");
                Integer b = bit.get(sid);
                if (b == null) continue;
                long[] entry = byVariant.computeIfAbsent(vid,
                        k -> new long[]{pos, 0L});
                entry[1] |= 1L << b;
            }
        } finally {
            r.close();
        }
        // MAF filter.
        int nSamples = ordered.size();
        double minCarriers = minMaf * nSamples;
        double maxCarriers = (1.0 - minMaf) * nSamples;
        List<long[]> out = new ArrayList<>(byVariant.size());
        for (long[] v : byVariant.values()) {
            int k = Long.bitCount(v[1]);
            if (k <= minCarriers || k >= maxCarriers) continue;
            out.add(v);
        }
        return out;
    }

    /**
     * Hill 1968 closed-form r² for biallelic SNPs given
     * carrier counts kA, kB, kAB and the haploid sample size n.
     */
    static double rSquared(int kA, int kB, int kAB, int n) {
        if (n <= 1) return Double.NaN;
        double pA = (double) kA / n;
        double pB = (double) kB / n;
        double pAB = (double) kAB / n;
        double d = pAB - pA * pB;
        double denom = pA * (1.0 - pA) * pB * (1.0 - pB);
        if (denom <= 0.0) return Double.NaN;
        return d * d / denom;
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
