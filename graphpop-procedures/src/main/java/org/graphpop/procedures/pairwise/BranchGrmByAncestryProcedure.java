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
 * Ancestry-decomposed branch GRM (M4.1 step 7, closes G3).
 *
 * <p>Buckets each branch's eGRM contribution by the ancestry painting
 * on the branch's child node. Sum across ancestries reproduces the
 * unconditional {@code branch_grm} matrix (modulo unpainted nodes).</p>
 *
 * <p>Conditional predicates from step 4
 * ({@code restrict_to_pathway}, {@code mutation_filter},
 * {@code time_window}) compose unchanged. Composition example:</p>
 *
 * <pre>
 * CALL graphpop.kinship.branch_grm_by_ancestry('singer_chr22_run_001', {
 *   restrict_to_pathway: 'GO:0006281',
 *   painter: 'majority_vote'
 * }) YIELD sample_a, sample_b, ancestry, b_ij_component
 * </pre>
 *
 * <p>The {@code painter} option selects which painting to use when a
 * run has multiple (e.g. RFMix vs Flare). Default: any painting.</p>
 */
public class BranchGrmByAncestryProcedure {

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.kinship.branch_grm_by_ancestry", mode = Mode.READ)
    @Description("Ancestry-decomposed branch GRM. Returns one row per "
            + "(sample_a, sample_b, ancestry) triple with the share of the "
            + "unconditional B_ij attributable to that ancestry. Closes "
            + "literature gap G3 (local-ancestry × ARG topology). "
            + "Conditional predicates from branch_grm compose unchanged.")
    @SuppressWarnings("unchecked")
    public Stream<AncestryKinshipResult> branchGrmByAncestry(
            @Name("run_id") String runId,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        if (options == null) options = new HashMap<>();
        long start = getLong(options, "start", 0L);
        long end = getLong(options, "end", Long.MAX_VALUE);
        boolean includeSelf = getBoolean(options, "include_self", true);
        String painter = (String) options.get("painter");

        ARG arg = ARGTraversal.load(tx, runId, start, end);
        if (arg.nNodes == 0 || arg.nSamples() == 0) return Stream.empty();

        Map<Integer, List<BranchGrmByAncestryComputer.AncestryAssignment>> painting =
                BranchGrmByAncestryComputer.loadPainting(tx, runId, painter);
        if (painting.isEmpty()) return Stream.empty();

        BranchWeightFn weight =
                BranchGrmConditioning.fromOptions(tx, runId, options);
        BranchGrmByAncestryComputer.Result result =
                BranchGrmByAncestryComputer.compute(arg, start, end, weight, painting);

        String[] sampleIds = resolveSampleIds(tx, runId, arg);
        long nBranches = arg.nEdges;

        List<AncestryKinshipResult> rows = new ArrayList<>();
        int n = result.n;
        for (String ancestry : result.ancestries) {
            double[][] m = result.matrices.get(ancestry);
            for (int i = 0; i < n; i++) {
                for (int j = i; j < n; j++) {
                    if (i == j && !includeSelf) continue;
                    rows.add(new AncestryKinshipResult(
                            sampleIds[i], sampleIds[j], ancestry,
                            m[i][j], nBranches, "branch_grm_by_ancestry"));
                }
            }
        }
        return rows.stream();
    }

    private static String[] resolveSampleIds(Transaction tx, String runId, ARG arg) {
        String[] out = new String[arg.nSamples()];
        for (int k = 0; k < arg.nSamples(); k++) {
            int packed = arg.sampleNodes[k];
            int tskitId = arg.tskitNodeId[packed];
            out[k] = runId + ":" + tskitId;
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

    private static boolean getBoolean(Map<String, Object> opts, String key, boolean def) {
        Object v = opts.get(key);
        return (v instanceof Boolean) ? (Boolean) v : def;
    }
}
