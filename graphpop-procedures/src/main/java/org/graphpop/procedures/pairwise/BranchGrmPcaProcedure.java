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
 * Top-K principal components of the branch GRM via Lanczos
 * iteration (M4.4). Builds on {@link BranchGrmMatVec#apply} so
 * memory is {@code O(K · n)} rather than {@code O(n²)}; tractable
 * at biobank scale.
 *
 * <pre>
 * CALL graphpop.kinship.branch_grm_pca(
 *   'tsinfer_chr22_v1',
 *   10,
 *   {n_iter: 30, seed: 42, restrict_to_pathway: 'GO:0006281'}
 * ) YIELD sample_id, pc, value, eigenvalue, method
 * </pre>
 *
 * <p>Conditional predicates (restrict_to_pathway, mutation_filter,
 * time_window) compose unchanged via
 * {@link BranchGrmConditioning#fromOptions}. Output is one row per
 * (sample, pc) triple.</p>
 */
public class BranchGrmPcaProcedure {

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.kinship.branch_grm_pca", mode = Mode.READ)
    @Description("Top-K principal components of branch GRM via Lanczos. "
            + "Uses Algorithm V mat-vec; tractable at biobank scale. "
            + "Conditional predicates compose.")
    @SuppressWarnings("unchecked")
    public Stream<PcaResult> branchGrmPca(
            @Name("run_id") String runId,
            @Name("k") long k,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        if (options == null) options = new HashMap<>();
        long start = getLong(options, "start", 0L);
        long end = getLong(options, "end", Long.MAX_VALUE);
        long seed = getLong(options, "seed", 42L);
        int nIter = (int) getLong(options, "n_iter", 0L);  // 0 = default 3*k

        ARG arg = ARGTraversal.load(tx, runId, start, end);
        if (arg.nNodes == 0 || arg.nSamples() == 0) return Stream.empty();
        int n = arg.nSamples();
        if (k < 1) {
            throw new RuntimeException("k must be >= 1; got " + k);
        }
        if (k >= n) {
            throw new RuntimeException(
                "k (" + k + ") must be < n_samples (" + n + ")");
        }

        BranchWeightFn weight =
                BranchGrmConditioning.fromOptions(tx, runId, options);
        BranchGrmLanczos.Result r = BranchGrmLanczos.runTopK(
                arg, start, end, weight,
                (int) k, nIter, seed, 1e-12);

        String[] sampleIds = resolveSampleIds(tx, runId, arg);

        List<PcaResult> rows = new ArrayList<>();
        for (int j = 0; j < r.eigenvalues.length; j++) {
            double[] pc = r.eigenvectors[j];
            double eigenvalue = r.eigenvalues[j];
            for (int i = 0; i < n; i++) {
                rows.add(new PcaResult(sampleIds[i], j, pc[i], eigenvalue,
                        "branch_grm_pca"));
            }
        }
        return rows.stream();
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
