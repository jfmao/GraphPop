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
 * Per-variant allele-age z-score scan (M8 part A).
 *
 * <p>Sweep signal: a variant whose allele age is "too young for its
 * frequency" is suspect for recent positive selection. The scan
 * stratifies all genome-wide variants by derived-allele frequency
 * (logit-spaced bins by default), computes per-bin
 * (mean, sd) of log-age, and emits a z-score per variant.</p>
 *
 * <p>Highly negative z-scores (younger than expected) flag candidate
 * sweeps; highly positive flag balanced selection or other
 * age-inflating signals.</p>
 *
 * <pre>
 * CALL graphpop.selection.allele_age_scan($runId,
 *         {n_freq_bins: 20, min_freq: 0.05, max_freq: 0.95})
 *   YIELD variant_id, freq, n_carriers, n_samples, age_midpoint,
 *         log_age, bin_index, bin_n, bin_mean_log_age,
 *         bin_sd_log_age, z_score, runId
 * </pre>
 */
public class AlleleAgeScanProcedure {

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.selection.allele_age_scan", mode = Mode.READ)
    @Description("Per-variant allele-age z-score conditional on derived-"
            + "allele frequency. Negative z = younger than expected for "
            + "frequency = candidate sweep.")
    @SuppressWarnings("unchecked")
    public Stream<AlleleAgeScanResult> alleleAgeScan(
            @Name("runId") String runId,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        if (options == null) options = new HashMap<>();
        int nBins = getInt(options, "n_freq_bins", 20);
        double minFreq = getDouble(options, "min_freq", 0.05);
        double maxFreq = getDouble(options, "max_freq", 0.95);

        ARG arg = ARGTraversal.load(tx, runId, 0L, Long.MAX_VALUE);
        if (arg.nNodes == 0) return Stream.empty();
        long nSamples = arg.nSamples();
        if (nSamples < 2) return Stream.empty();

        // Pull every variant + its mutated_on edge for this run, with
        // child/parent times. Carrier counts come from descendant counts
        // at the variant's site in its marginal tree.
        List<MutationRow> muts = loadMutations(tx, runId, arg);
        if (muts.isEmpty()) return Stream.empty();

        // First pass: per-bin Welford accumulator over log_age.
        double[] edges = SelectionUtils.logitFreqEdges(nBins, minFreq, maxFreq);
        SelectionUtils.Welford[] bins = new SelectionUtils.Welford[nBins];
        for (int i = 0; i < nBins; i++) bins[i] = new SelectionUtils.Welford();

        for (MutationRow m : muts) {
            int bi = SelectionUtils.binIndex(edges, m.freq);
            bins[bi].add(m.logAge);
        }

        // Second pass: emit z-score per variant.
        List<AlleleAgeScanResult> rows = new ArrayList<>(muts.size());
        for (MutationRow m : muts) {
            int bi = SelectionUtils.binIndex(edges, m.freq);
            SelectionUtils.Welford bin = bins[bi];
            double z = SelectionUtils.zScore(m.logAge, bin);
            rows.add(new AlleleAgeScanResult(
                    m.variantId, m.freq, m.nCarriers, nSamples,
                    m.ageMidpoint, m.logAge,
                    bin.mean(), bin.sd(),
                    (long) bi, bin.n(),
                    z, runId));
        }
        return rows.stream();
    }

    /** Per-mutation row gathered from the graph. */
    private static final class MutationRow {
        final String variantId;
        final long position;
        final double childTime;
        final double parentTime;
        final double ageMidpoint;
        final double logAge;
        final long nCarriers;
        final double freq;

        MutationRow(String variantId, long position,
                    double childTime, double parentTime,
                    long nCarriers, long nSamples) {
            this.variantId = variantId;
            this.position = position;
            this.childTime = childTime;
            this.parentTime = parentTime;
            this.ageMidpoint = 0.5 * (childTime + parentTime);
            this.logAge = Math.log(Math.max(this.ageMidpoint, 1e-300));
            this.nCarriers = nCarriers;
            this.freq = (double) nCarriers / (double) nSamples;
        }
    }

    private static List<MutationRow> loadMutations(Transaction tx,
                                                    String runId, ARG arg) {
        // Cypher: (:Variant)-[m:MUTATED_ON]->(:TreeNode child) where
        // m.parent_node_id gives the upper end of the bracket. We then
        // count descendants of `child` in the marginal tree at the
        // variant's position to get the carrier count.
        Map<Integer, Integer> idMap = arg.tskitIdToIndexMap();

        Result r = tx.execute(
                "MATCH (v:Variant)-[m:MUTATED_ON {runId: $runId}]->"
              + "(c:TreeNode) "
              + "RETURN v.variantId AS variantId, v.position AS position, "
              + "c.nodeId AS childId, c.time AS childTime, "
              + "m.parent_node_id AS parentId",
                Map.of("runId", runId));
        List<MutationRow> out = new ArrayList<>();
        try {
            // Cache parent times to avoid repeated lookups.
            Map<Long, Double> parentTimeCache = new HashMap<>();
            while (r.hasNext()) {
                Map<String, Object> row = r.next();
                String vid = (String) row.get("variantId");
                long pos = ((Number) row.get("position")).longValue();
                int childTskit = ((Number) row.get("childId")).intValue();
                double childTime = ((Number) row.get("childTime")).doubleValue();
                long parentTskit = ((Number) row.get("parentId")).longValue();

                Double pt = parentTimeCache.get(parentTskit);
                if (pt == null) {
                    Result pr = tx.execute(
                            "MATCH (n:TreeNode {runId: $rid, nodeId: $pid}) "
                          + "RETURN n.time AS t",
                            Map.of("rid", runId, "pid", parentTskit));
                    if (pr.hasNext()) {
                        pt = ((Number) pr.next().get("t")).doubleValue();
                    }
                    pr.close();
                    if (pt != null) parentTimeCache.put(parentTskit, pt);
                }
                if (pt == null) continue;  // missing parent — skip

                // Carrier count: descendants of child in the marginal tree
                // at the variant's position.
                Integer childPacked = idMap.get(childTskit);
                if (childPacked == null) continue;
                int[] parent = ArgTraversalUtils.marginalTree(arg, pos);
                long nCarriers = countSampleDescendants(arg, parent, childPacked);
                if (nCarriers <= 0) continue;

                out.add(new MutationRow(vid, pos, childTime, pt,
                        nCarriers, arg.nSamples()));
            }
        } finally {
            r.close();
        }
        return out;
    }

    private static long countSampleDescendants(ARG arg, int[] parent, int root) {
        long count = 0L;
        for (int leaf : arg.sampleNodes) {
            int cur = leaf;
            while (cur != -1) {
                if (cur == root) { count++; break; }
                cur = parent[cur];
            }
        }
        return count;
    }

    private static int getInt(Map<String, Object> opts, String key, int def) {
        Object v = opts.get(key);
        return (v instanceof Number) ? ((Number) v).intValue() : def;
    }

    private static double getDouble(Map<String, Object> opts, String key,
                                    double def) {
        Object v = opts.get(key);
        return (v instanceof Number) ? ((Number) v).doubleValue() : def;
    }
}
