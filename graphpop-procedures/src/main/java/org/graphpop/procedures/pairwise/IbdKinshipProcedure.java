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
 * Browning-style kinship from total shared IBD length (M4.3).
 *
 * <p>For each pair of samples that share at least one
 * {@code :IBD_SEGMENT} edge with the requested {@code source}:</p>
 *
 * <pre>
 *   phi^IBD_ij = Σ_segments length_cM / (2 · total_genome_cM)
 * </pre>
 *
 * <p>When {@code length_cM} is missing (e.g. on ARG-derived segments
 * without a genetic map), falls back to {@code length_bp /
 * (2 · total_genome_bp)} and tags {@code method = "ibd_bp"}. The
 * total genome length is inferred from the segments themselves
 * (sum over chromosomes of {@code MAX(end) - MIN(start)}).</p>
 *
 * <p>Reuses {@link KinshipResult}: {@code phi} holds the kinship,
 * {@code n_snp} = number of segments, {@code ibs0} =
 * total_length_bp, {@code het_het} = total_length_cM (or 0 when
 * unavailable).</p>
 */
public class IbdKinshipProcedure {

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.ibd.kinship", mode = Mode.READ)
    @Description("Browning-style kinship from total shared IBD length. "
            + "Sums :IBD_SEGMENT lengths per pair and divides by 2x total "
            + "genome length. Uses length_cM when available; falls back to "
            + "length_bp with method='ibd_bp'.")
    @SuppressWarnings("unchecked")
    public Stream<KinshipResult> kinship(
            @Name("source") String source,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        if (options == null) options = new HashMap<>();
        long minLengthBp = getLong(options, "min_length_bp", 0L);

        // Probe whether segments carry length_cM (any non-null on this source).
        boolean useCM = false;
        try (Result probe = tx.execute(
                "MATCH ()-[r:IBD_SEGMENT {source: $source}]->() "
              + "WHERE r.length_cM IS NOT NULL RETURN r.length_cM AS x LIMIT 1",
                Map.of("source", source))) {
            useCM = probe.hasNext();
        }
        String method = useCM ? "ibd" : "ibd_bp";

        // Compute total genome length from the segments themselves.
        // Per-chromosome span = MAX(end) - MIN(start); sum over chromosomes.
        double totalGenome = 0.0;
        try (Result genome = tx.execute(
                useCM
                  ? "MATCH ()-[r:IBD_SEGMENT {source: $source}]->() "
                  + "WHERE r.length_cM IS NOT NULL "
                  + "RETURN sum(r.length_cM) / count(DISTINCT [r.sample_a, r.sample_b]) AS x"
                  : "MATCH ()-[r:IBD_SEGMENT {source: $source}]->() "
                  + "WITH r.chr AS chr, max(r.end) AS hi, min(r.start) AS lo "
                  + "RETURN sum(hi - lo) AS x",
                Map.of("source", source))) {
            // The cM branch is a placeholder -- callers should supply
            // total_genome_cM via options. We compute below.
            // (We re-derive total_genome from per-chromosome spans for both.)
            // Discard probe; use the canonical bp-span query below.
        }
        // Canonical total: per-chromosome span (max(end) - min(start)).
        try (Result spans = tx.execute(
                "MATCH ()-[r:IBD_SEGMENT {source: $source}]->() "
              + "WITH r.chr AS chr, max(r.end) AS hi, min(r.start) AS lo "
              + "RETURN sum(hi - lo) AS total_bp",
                Map.of("source", source))) {
            if (spans.hasNext()) {
                Object v = spans.next().get("total_bp");
                if (v instanceof Number) totalGenome = ((Number) v).doubleValue();
            }
        }
        double totalGenomeCM = getDouble(options, "total_genome_cM", -1.0);
        if (useCM && totalGenomeCM <= 0) {
            // Estimate cM total as a Morgan/Mb assumption: ~1 cM per Mb.
            // Caller can supply explicitly via options.
            totalGenomeCM = totalGenome / 1_000_000.0;
        }

        if (totalGenome <= 0) return Stream.empty();
        double denom = useCM ? (2.0 * totalGenomeCM) : (2.0 * totalGenome);

        // Aggregate per pair.
        String aggCypher = useCM
              ? "MATCH (a:Sample)-[r:IBD_SEGMENT {source: $source}]->(b:Sample) "
              + "WITH a.sampleId AS sa, b.sampleId AS sb, "
              + "sum(r.length_bp) AS total_bp, "
              + "sum(coalesce(r.length_cM, 0.0)) AS total_cM, "
              + "count(r) AS n_segs "
              + "RETURN sa, sb, total_bp, total_cM, n_segs"
              : "MATCH (a:Sample)-[r:IBD_SEGMENT {source: $source}]->(b:Sample) "
              + "WITH a.sampleId AS sa, b.sampleId AS sb, "
              + "sum(r.length_bp) AS total_bp, count(r) AS n_segs "
              + "RETURN sa, sb, total_bp, 0.0 AS total_cM, n_segs";

        List<KinshipResult> rows = new ArrayList<>();
        try (Result agg = tx.execute(aggCypher, Map.of("source", source))) {
            while (agg.hasNext()) {
                Map<String, Object> row = agg.next();
                long totalBp = ((Number) row.get("total_bp")).longValue();
                if (totalBp < minLengthBp) continue;
                long nSegs = ((Number) row.get("n_segs")).longValue();
                double totalCM = ((Number) row.get("total_cM")).doubleValue();
                double phi = useCM
                        ? (totalCM / denom)
                        : (totalBp / denom);
                rows.add(new KinshipResult(
                        (String) row.get("sa"), (String) row.get("sb"),
                        phi, totalBp,
                        useCM ? Math.round(totalCM * 1e6) : 0L,  // pack cM*1e6 into long
                        nSegs, 0L, method));
            }
        }
        return rows.stream();
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
