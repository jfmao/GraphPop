package org.graphpop.procedures.pairwise;

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
 * Algorithm V matrix-vector form of {@code branch_grm}: returns
 * {@code G · v} without materialising {@code G}, in
 * {@code O(T · N log N)} per multiplication. Lifts the v1 full-matrix
 * cap and unblocks UK-Biobank-class cohorts.
 *
 * <pre>
 * CALL graphpop.kinship.branch_grm_apply(
 *   'tsinfer_chr22_v1',
 *   [0.5, -0.3, 0.1, ...],          -- length = n_samples
 *   {restrict_to_pathway: 'GO:0006281'}
 * ) YIELD sample_id, col, value, n_branches, method
 * </pre>
 *
 * <p>Multi-vector input is supported via the {@code vectors} option
 * (a {@code List&lt;List&lt;Double&gt;&gt;} of equal-length vectors).
 * Each input column produces one result column tagged by 0-based
 * {@code col}.</p>
 *
 * <p>Conditional predicates ({@code restrict_to_pathway},
 * {@code mutation_filter}, {@code time_window}) compose unchanged via
 * {@link BranchGrmConditioning#fromOptions}.</p>
 */
public class BranchGrmApplyProcedure {

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.kinship.branch_grm_apply", mode = Mode.READ)
    @Description("Apply the branch GRM to a vector (Algorithm V). Returns "
            + "G·v without materialising G; tractable at biobank scale. "
            + "Conditional predicates from branch_grm compose unchanged.")
    @SuppressWarnings("unchecked")
    public Stream<KinshipVectorResult> branchGrmApply(
            @Name("run_id") String runId,
            @Name(value = "vector", defaultValue = "[]") List<Double> vector,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        if (options == null) options = new HashMap<>();
        long start = getLong(options, "start", 0L);
        long end = getLong(options, "end", Long.MAX_VALUE);

        ARG arg = ARGTraversal.load(tx, runId, start, end);
        if (arg.nNodes == 0 || arg.nSamples() == 0) return Stream.empty();

        BranchWeightFn weight =
                BranchGrmConditioning.fromOptions(tx, runId, options);

        // Resolve sample order before parsing the vector(s) so we can
        // validate length up-front.
        String[] sampleIds = resolveSampleIds(tx, runId, arg);
        int n = sampleIds.length;

        Object multi = options.get("vectors");
        double[][] inputVectors = parseVectors(vector, multi, n);

        double[][] gv = new double[inputVectors.length][];
        for (int k = 0; k < inputVectors.length; k++) {
            gv[k] = BranchGrmMatVec.apply(arg, start, end, weight, inputVectors[k]);
        }

        long nBranches = arg.nEdges;
        List<KinshipVectorResult> rows = new ArrayList<>();
        for (int k = 0; k < gv.length; k++) {
            double[] r = gv[k];
            for (int i = 0; i < n; i++) {
                rows.add(new KinshipVectorResult(
                        sampleIds[i], k, r[i], nBranches,
                        "branch_grm_apply"));
            }
        }
        return rows.stream();
    }

    @SuppressWarnings("unchecked")
    private static double[][] parseVectors(List<Double> single, Object multi, int n) {
        if (multi != null) {
            List<List<Number>> raw = (List<List<Number>>) multi;
            double[][] out = new double[raw.size()][];
            for (int k = 0; k < raw.size(); k++) {
                List<Number> v = raw.get(k);
                if (v.size() != n) {
                    throw new RuntimeException(
                        "vectors[" + k + "] length " + v.size()
                      + " != n_samples " + n);
                }
                double[] d = new double[n];
                for (int i = 0; i < n; i++) d[i] = v.get(i).doubleValue();
                out[k] = d;
            }
            return out;
        }
        if (single == null || single.isEmpty()) {
            throw new RuntimeException(
                "either 'vector' (length " + n + ") or 'vectors' "
              + "(list of length-" + n + " vectors) must be supplied");
        }
        if (single.size() != n) {
            throw new RuntimeException(
                "vector length " + single.size() + " != n_samples " + n);
        }
        double[] d = new double[n];
        for (int i = 0; i < n; i++) d[i] = single.get(i).doubleValue();
        return new double[][]{d};
    }

    private static String[] resolveSampleIds(Transaction tx, String runId, ARG arg) {
        String[] out = new String[arg.nSamples()];
        for (int k = 0; k < arg.nSamples(); k++) {
            int packed = arg.sampleNodes[k];
            out[k] = runId + ":" + arg.tskitNodeId[packed];
        }
        Result r = tx.execute(
                "MATCH (n:TreeNode {runId: $runId})-[rel:REPRESENTS]->(s:Sample) "
              + "RETURN n.nodeId AS nodeId, s.sampleId AS sampleId, rel.haplotype AS haplotype",
                Map.of("runId", runId));
        try {
            Map<Integer, String> tskitToSample = new HashMap<>();
            while (r.hasNext()) {
                Map<String, Object> row = r.next();
                int nodeId = ((Number) row.get("nodeId")).intValue();
                String sid = (String) row.get("sampleId");
                Object hap = row.get("haplotype");
                String label = (hap != null) ? sid + ":h" + ((Number) hap).intValue() : sid;
                tskitToSample.put(nodeId, label);
            }
            for (int k = 0; k < arg.nSamples(); k++) {
                int tskitId = arg.tskitNodeId[arg.sampleNodes[k]];
                String mapped = tskitToSample.get(tskitId);
                if (mapped != null) out[k] = mapped;
            }
        } finally {
            r.close();
        }
        return out;
    }

    private static long getLong(Map<String, Object> opts, String key, long def) {
        Object v = opts.get(key);
        return (v instanceof Number) ? ((Number) v).longValue() : def;
    }
}
