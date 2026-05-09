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
 * Unconditional branch GRM via {@code :TreeNode}/{@code :PARENT_OF}
 * traversal (M4.1 step 3).
 *
 * <p>Algorithmic core: {@link BranchGrmComputer}, validated to match
 * {@code egrm.varGRM(ts)} (Fan, Mancuso &amp; Chiang 2022 <i>AJHG</i>) on
 * msprime fixtures to relative error &lt; 10⁻⁶.</p>
 *
 * <pre>
 * CALL graphpop.kinship.branch_grm('tsinfer_chr22_v1', {start: 1, end: 50000000})
 *   YIELD sample_a, sample_b, phi, ibs0, het_het, n_snp, n_aa_min, method
 * </pre>
 *
 * <p>The {@code phi} column carries B<sub>ij</sub>; the unused KingRobust
 * fields ({@code ibs0}, {@code het_het}, {@code n_aa_min}) are reported as
 * 0 since they have no analogue in the branch-GRM formulation. {@code n_snp}
 * is overloaded to mean {@code n_trees} (number of marginal trees that
 * contributed); {@code method} is {@code "branch_grm"}.</p>
 */
public class BranchGrmProcedure {

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.kinship.branch_grm", mode = Mode.READ)
    @Description("Unconditional branch GRM (Fan, Mancuso & Chiang 2022) via "
            + "ARG traversal. Validation baseline: egrm.varGRM Python package.")
    @SuppressWarnings("unchecked")
    public Stream<KinshipResult> branchGrm(
            @Name("run_id") String runId,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        if (options == null) options = new HashMap<>();

        long start = getLong(options, "start", 0L);
        long end = getLong(options, "end", Long.MAX_VALUE);
        boolean includeSelf = getBoolean(options, "include_self", true);

        ARG arg = ARGTraversal.load(tx, runId, start, end);
        if (arg.nNodes == 0 || arg.nSamples() == 0) return Stream.empty();

        BranchWeightFn weight = BranchGrmConditioning.fromOptions(tx, runId, options);
        BranchGrmComputer.Result result =
                BranchGrmComputer.compute(arg, start, end, weight);

        String[] sampleIds = resolveSampleIds(tx, runId, arg);
        long nTrees = countTrees(arg, start, end);

        List<KinshipResult> rows = new ArrayList<>();
        int n = result.n;
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                if (i == j && !includeSelf) continue;
                rows.add(new KinshipResult(
                        sampleIds[i], sampleIds[j],
                        result.matrix[i][j],
                        0L, 0L,
                        nTrees, 0L,
                        "branch_grm"));
            }
        }
        return rows.stream();
    }

    /**
     * For each sample-flagged TreeNode (in {@link ARG#sampleNodes} order),
     * resolve its {@code :Sample.sampleId} via {@code :REPRESENTS}. Falls
     * back to {@code "{runId}:{tskit_node_id}"} when no edge is present
     * (test fixtures, partial ingests).
     */
    private static String[] resolveSampleIds(Transaction tx, String runId, ARG arg) {
        String[] out = new String[arg.nSamples()];
        for (int k = 0; k < arg.nSamples(); k++) {
            int packed = arg.sampleNodes[k];
            int tskitId = arg.tskitNodeId[packed];
            String tnodeId = runId + ":" + tskitId;
            out[k] = tnodeId;  // default fallback
        }
        Result r = tx.execute(
                "MATCH (n:TreeNode {runId: $runId})-[rel:REPRESENTS]->(s:Sample) "
              + "RETURN n.nodeId AS nodeId, s.sampleId AS sampleId, rel.haplotype AS haplotype",
                Map.of("runId", runId));
        try {
            // nodeId -> sampleId map
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

    /** Number of marginal-tree intervals that overlap {@code [start, end)}. */
    private static long countTrees(ARG arg, long start, long end) {
        long lo = Math.max(start, 0L);
        long hi = (end == Long.MAX_VALUE) ? arg.sequenceLength : Math.min(end, arg.sequenceLength);
        long[] bp = arg.breakpointsClipped(lo, hi);
        return Math.max(0, bp.length - 1);
    }

    private static long getLong(Map<String, Object> opts, String key, long def) {
        Object v = opts.get(key);
        if (v instanceof Number) return ((Number) v).longValue();
        return def;
    }

    private static boolean getBoolean(Map<String, Object> opts, String key, boolean def) {
        Object v = opts.get(key);
        if (v instanceof Boolean) return (Boolean) v;
        return def;
    }
}
