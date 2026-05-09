package org.graphpop.procedures.pairwise;

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
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

/**
 * Posterior aggregator for {@code branch_grm}: computes the
 * posterior-mean kinship and standard error of the mean over a set of
 * {@code :ARGRun}s (e.g. M=100 SINGER posterior samples). Closes the
 * G1 (ARG-inference uncertainty) gap in full.
 *
 * <p>Conditional predicates from M4.1 step 4 ({@code restrict_to_pathway},
 * {@code mutation_filter}, {@code time_window}) compose unchanged.
 * Per-run predicate state is rebuilt for each run because the
 * lit-child mask depends on per-run {@code :MUTATED_ON} mappings.</p>
 *
 * <pre>
 * CALL graphpop.kinship.branch_grm_posterior(
 *   ['singer_chr22_run_001', ..., 'singer_chr22_run_100'],
 *   {restrict_to_pathway: 'GO:0006281', time_window: [0, 1000]}
 * ) YIELD sample_a, sample_b, b_ij_mean, b_ij_sd, n_runs, method
 * </pre>
 */
public class PosteriorBranchGrmProcedure {

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.kinship.branch_grm_posterior", mode = Mode.READ)
    @Description("Posterior aggregator for branch_grm: emits posterior-mean "
            + "and standard-error-of-mean kinship over a set of :ARGRun "
            + "(e.g. SINGER posterior samples). Conditional predicates "
            + "(restrict_to_pathway, mutation_filter, time_window) compose.")
    @SuppressWarnings("unchecked")
    public Stream<PosteriorKinshipResult> branchGrmPosterior(
            @Name("posterior_run_ids") List<String> posteriorRunIds,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        if (posteriorRunIds == null || posteriorRunIds.isEmpty()) {
            return Stream.empty();
        }
        if (options == null) options = new HashMap<>();

        long start = getLong(options, "start", 0L);
        long end = getLong(options, "end", Long.MAX_VALUE);
        boolean includeSelf = getBoolean(options, "include_self", true);

        WelfordMatrix welford = null;
        String[] referenceSampleIds = null;

        for (String runId : posteriorRunIds) {
            ARG arg = ARGTraversal.load(tx, runId, start, end);
            if (arg.nNodes == 0 || arg.nSamples() == 0) {
                throw new RuntimeException(
                    "ARGRun '" + runId + "' has no TreeNodes or no samples; "
                  + "cannot include in posterior aggregator");
            }

            String[] sampleIds = resolveSampleIds(tx, runId, arg);

            if (referenceSampleIds == null) {
                referenceSampleIds = sampleIds;
                welford = new WelfordMatrix(sampleIds.length);
            } else {
                if (!Arrays.equals(referenceSampleIds, sampleIds)) {
                    throw new RuntimeException(
                        "Sample alignment mismatch between run '"
                      + posteriorRunIds.get(0) + "' and run '" + runId
                      + "'. Posterior aggregator requires identical "
                      + "ordered sample sets across all runs.");
                }
            }

            BranchWeightFn weight =
                    BranchGrmConditioning.fromOptions(tx, runId, options);
            BranchGrmComputer.Result r =
                    BranchGrmComputer.compute(arg, start, end, weight);
            welford.update(r.matrix);
        }

        double[][] mean = welford.mean();
        double[][] se = welford.stderr();
        long nRuns = posteriorRunIds.size();

        List<PosteriorKinshipResult> rows = new ArrayList<>();
        int n = referenceSampleIds.length;
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                if (i == j && !includeSelf) continue;
                rows.add(new PosteriorKinshipResult(
                        referenceSampleIds[i], referenceSampleIds[j],
                        mean[i][j], se[i][j],
                        nRuns, "branch_grm_posterior"));
            }
        }
        return rows.stream();
    }

    /**
     * Same as {@code BranchGrmProcedure#resolveSampleIds}: maps the
     * sample-flagged TreeNodes (in {@link ARG#sampleNodes} order) to
     * their {@code :Sample.sampleId} via {@code :REPRESENTS}, falling
     * back to the {@code "{runId}:{tskit_node_id}"} format.
     */
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
