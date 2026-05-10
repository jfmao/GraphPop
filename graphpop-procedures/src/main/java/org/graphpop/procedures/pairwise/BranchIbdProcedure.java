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
 * ARG-derived IBD-segment caller (M4.3).
 *
 * <p>Walks {@code :TreeNode}/{@code :PARENT_OF} via {@link ARGTraversal}
 * and uses {@link BranchIbdComputer} to extract maximal IBD segments
 * per pair under an optional TMRCA cap. Equivalent to
 * {@code tskit.TreeSequence.ibd_segments(within=samples, max_time=...,
 * min_span=...)}.</p>
 *
 * <p>Mode.WRITE: emits {@code :IBD_SEGMENT} edges between the
 * corresponding {@code :Sample} nodes (via {@code :REPRESENTS}) AND
 * streams one row per segment back to the caller. Existing
 * {@code (sample_a)-[:IBD_SEGMENT {source: 'arg_derived', runId: $rid}]->(sample_b)}
 * edges for the same {@code runId} are deleted before re-insert
 * (idempotent re-runs).</p>
 *
 * <pre>
 * CALL graphpop.ibd.from_arg('singer_chr22_run_001', {
 *   max_tmrca: 50.0,
 *   min_length_bp: 2000000,
 *   chr: 'chr22'
 * }) YIELD sample_a, sample_b, chr, start, end, length_bp,
 *           mrca_node_id, tmrca, source, runId
 * </pre>
 */
public class BranchIbdProcedure {

    private static final String SOURCE = "arg_derived";

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.ibd.from_arg", mode = Mode.WRITE)
    @Description("ARG-derived IBD segments. Walks :PARENT_OF to extract "
            + "maximal MRCA intervals per pair under an optional TMRCA cap. "
            + "Emits :IBD_SEGMENT edges (idempotent per runId) and streams "
            + "rows back. Equivalent to tskit's ibd_segments.")
    @SuppressWarnings("unchecked")
    public Stream<IbdSegmentResult> fromArg(
            @Name("run_id") String runId,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        if (options == null) options = new HashMap<>();

        long start = getLong(options, "start", 0L);
        long end = getLong(options, "end", Long.MAX_VALUE);
        double maxTmrca = getDouble(options, "max_tmrca",
                BranchIbdComputer.NO_MAX_TMRCA);
        long minLengthBp = getLong(options, "min_length_bp",
                BranchIbdComputer.NO_MIN_LENGTH);
        String chr = (String) options.getOrDefault("chr", "chr1");

        ARG arg = ARGTraversal.load(tx, runId, start, end);
        if (arg.nNodes == 0 || arg.nSamples() == 0) return Stream.empty();

        List<IbdSegment> segments = BranchIbdComputer.compute(
                arg, start, end, maxTmrca, minLengthBp);
        String[] sampleIds = resolveSampleIds(tx, runId, arg);

        // Delete existing arg-derived edges for this runId (idempotent re-run).
        tx.execute(
                "MATCH ()-[r:IBD_SEGMENT {source: $source, runId: $runId}]->() "
              + "DELETE r",
                Map.of("source", SOURCE, "runId", runId));

        // Emit edges in batches.
        List<IbdSegmentResult> rows = new ArrayList<>(segments.size());
        List<Map<String, Object>> batch = new ArrayList<>();
        long now = System.currentTimeMillis();
        for (IbdSegment s : segments) {
            String sa = sampleIds[s.sampleA];
            String sb = sampleIds[s.sampleB];
            long lengthBp = s.lengthBp();
            rows.add(new IbdSegmentResult(
                    sa, sb, chr, s.start, s.end, lengthBp,
                    s.mrcaNodeId, s.tmrca, SOURCE, runId));

            Map<String, Object> row = new HashMap<>();
            row.put("sample_a", sa);
            row.put("sample_b", sb);
            row.put("chr", chr);
            row.put("start", s.start);
            row.put("end", s.end);
            row.put("length_bp", lengthBp);
            row.put("mrca_node_id", (long) s.mrcaNodeId);
            row.put("tmrca", s.tmrca);
            batch.add(row);
            if (batch.size() >= 10_000) {
                writeBatch(tx, batch, runId, now);
                batch.clear();
            }
        }
        if (!batch.isEmpty()) writeBatch(tx, batch, runId, now);

        return rows.stream();
    }

    private static void writeBatch(Transaction tx, List<Map<String, Object>> batch,
                                    String runId, long createdMillis) {
        tx.execute(
                "UNWIND $rows AS r "
              + "MATCH (a:Sample {sampleId: r.sample_a}), "
              + "(b:Sample {sampleId: r.sample_b}) "
              + "CREATE (a)-[:IBD_SEGMENT {"
              + "  chr: r.chr, start: r.start, end: r.end, "
              + "  length_bp: r.length_bp, mrca_node_id: r.mrca_node_id, "
              + "  tmrca: r.tmrca, source: $source, runId: $runId, "
              + "  created_at: datetime({epochMillis: $createdMillis})"
              + "}]->(b)",
                Map.of("rows", batch, "source", SOURCE,
                       "runId", runId, "createdMillis", createdMillis));
    }

    private static String[] resolveSampleIds(Transaction tx, String runId, ARG arg) {
        String[] out = new String[arg.nSamples()];
        for (int k = 0; k < arg.nSamples(); k++) {
            int packed = arg.sampleNodes[k];
            out[k] = runId + ":" + arg.tskitNodeId[packed];
        }
        Result r = tx.execute(
                "MATCH (n:TreeNode {runId: $runId})-[rel:REPRESENTS]->(s:Sample) "
              + "RETURN n.nodeId AS nodeId, s.sampleId AS sampleId",
                Map.of("runId", runId));
        try {
            Map<Integer, String> tskitToSample = new HashMap<>();
            while (r.hasNext()) {
                Map<String, Object> row = r.next();
                int nodeId = ((Number) row.get("nodeId")).intValue();
                String sid = (String) row.get("sampleId");
                tskitToSample.put(nodeId, sid);
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

    private static double getDouble(Map<String, Object> opts, String key, double def) {
        Object v = opts.get(key);
        return (v instanceof Number) ? ((Number) v).doubleValue() : def;
    }
}
