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
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

/**
 * Family clustering via connected components on a thresholded
 * relatedness graph (M5 Part B).
 *
 * <pre>
 * CALL graphpop.relate.families('king', {max_degree: 2})
 *   YIELD sample_id, family_id, family_size, method
 * </pre>
 *
 * <p>Edge predicate (mutually exclusive — pick one):</p>
 * <ul>
 *   <li>{@code max_degree} (default 2) — link two samples when there
 *       is a {@code :RELATIVE {source: $source, degree <= max_degree}}
 *       edge between them.</li>
 *   <li>{@code min_total_ibd_bp} — link two samples when their summed
 *       {@code :IBD_SEGMENT {source: $source}} length is at least the
 *       given threshold.</li>
 *   <li>{@code min_phi} — link two samples when their
 *       {@code :RELATIVE.phi} is at least the given threshold.</li>
 * </ul>
 *
 * <p>With {@code persist: true} (default), writes
 * {@code (:Sample)-[:IN_FAMILY {source}]->(:Family {family_id})}
 * idempotently per source.</p>
 */
public class RelateFamiliesProcedure {

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.relate.families", mode = Mode.WRITE)
    @Description("Family clustering via connected components on a "
            + "thresholded relatedness graph. Edge predicate: max_degree | "
            + "min_total_ibd_bp | min_phi (mutually exclusive).")
    public Stream<FamilyResult> families(
            @Name("source") String source,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        if (options == null) options = new HashMap<>();
        boolean persist = getBoolean(options, "persist", true);

        // Predicate selection.
        Long maxDegree = getNullableLong(options, "max_degree");
        Long minTotalIbdBp = getNullableLong(options, "min_total_ibd_bp");
        Double minPhi = getNullableDouble(options, "min_phi");

        String method;
        List<String[]> edges;
        if (minTotalIbdBp != null) {
            method = "ibd_bp";
            edges = loadEdgesByIbdBp(tx, source, minTotalIbdBp);
        } else if (minPhi != null) {
            method = "phi";
            edges = loadEdgesByPhi(tx, source, minPhi);
        } else {
            // Default predicate: max_degree (2 = "extended family or closer").
            long md = (maxDegree != null) ? maxDegree : 2L;
            method = "degree";
            edges = loadEdgesByDegree(tx, source, md);
        }

        // Pull every :Sample to ensure singletons appear in the output.
        List<String> samples = loadSamples(tx);
        List<FamilyClusterer.Assignment> assignments =
                FamilyClusterer.cluster(samples, edges);

        // Idempotent: clear existing :IN_FAMILY for this source.
        if (persist) {
            tx.execute(
                "MATCH (:Sample)-[r:IN_FAMILY {source: $source}]->(:Family) "
              + "DELETE r",
                Map.of("source", source));
            tx.execute(
                "MATCH (f:Family {source: $source}) "
              + "WHERE NOT (:Sample)-[:IN_FAMILY]->(f) DELETE f",
                Map.of("source", source));
        }

        long now = System.currentTimeMillis();
        if (persist) {
            // Distinct family IDs for this source.
            Set<String> familyIds = new HashSet<>();
            for (FamilyClusterer.Assignment a : assignments) {
                familyIds.add(a.familyId);
            }
            List<Map<String, Object>> familyRows = new ArrayList<>(familyIds.size());
            for (String fid : familyIds) {
                Map<String, Object> r = new HashMap<>();
                r.put("family_id", fid);
                familyRows.add(r);
            }
            if (!familyRows.isEmpty()) {
                tx.execute(
                    "UNWIND $rows AS r "
                  + "MERGE (f:Family {family_id: r.family_id, source: $source}) "
                  + "ON CREATE SET f.created_at = datetime({epochMillis: $now})",
                    Map.of("rows", familyRows, "source", source, "now", now));
            }

            List<Map<String, Object>> edgeRows = new ArrayList<>(assignments.size());
            for (FamilyClusterer.Assignment a : assignments) {
                Map<String, Object> r = new HashMap<>();
                r.put("sample_id", a.sampleId);
                r.put("family_id", a.familyId);
                r.put("family_size", (long) a.familySize);
                edgeRows.add(r);
            }
            tx.execute(
                "UNWIND $rows AS r "
              + "MATCH (s:Sample {sampleId: r.sample_id}), "
              + "(f:Family {family_id: r.family_id, source: $source}) "
              + "CREATE (s)-[:IN_FAMILY {source: $source, "
              + "family_size: r.family_size, "
              + "created_at: datetime({epochMillis: $now})}]->(f)",
                Map.of("rows", edgeRows, "source", source, "now", now));
        }

        List<FamilyResult> rows = new ArrayList<>(assignments.size());
        for (FamilyClusterer.Assignment a : assignments) {
            rows.add(new FamilyResult(a.sampleId, a.familyId,
                    a.familySize, method));
        }
        return rows.stream();
    }

    private static List<String> loadSamples(Transaction tx) {
        List<String> out = new ArrayList<>();
        Result r = tx.execute("MATCH (s:Sample) RETURN s.sampleId AS sid");
        try {
            while (r.hasNext()) out.add((String) r.next().get("sid"));
        } finally {
            r.close();
        }
        return out;
    }

    private static List<String[]> loadEdgesByDegree(Transaction tx,
                                                     String source,
                                                     long maxDegree) {
        List<String[]> out = new ArrayList<>();
        Result r = tx.execute(
                "MATCH (a:Sample)-[r:RELATIVE {source: $source}]->(b:Sample) "
              + "WHERE r.degree <= $maxDegree "
              + "RETURN a.sampleId AS a, b.sampleId AS b",
                Map.of("source", source, "maxDegree", maxDegree));
        try {
            while (r.hasNext()) {
                Map<String, Object> row = r.next();
                out.add(new String[]{(String) row.get("a"), (String) row.get("b")});
            }
        } finally {
            r.close();
        }
        return out;
    }

    private static List<String[]> loadEdgesByIbdBp(Transaction tx,
                                                    String source,
                                                    long minTotalIbdBp) {
        List<String[]> out = new ArrayList<>();
        Result r = tx.execute(
                "MATCH (a:Sample)-[s:IBD_SEGMENT {source: $source}]->(b:Sample) "
              + "WITH a.sampleId AS sa, b.sampleId AS sb, "
              + "sum(s.length_bp) AS total_bp "
              + "WHERE total_bp >= $threshold "
              + "RETURN sa, sb",
                Map.of("source", source, "threshold", minTotalIbdBp));
        try {
            while (r.hasNext()) {
                Map<String, Object> row = r.next();
                out.add(new String[]{(String) row.get("sa"), (String) row.get("sb")});
            }
        } finally {
            r.close();
        }
        return out;
    }

    private static List<String[]> loadEdgesByPhi(Transaction tx,
                                                  String source,
                                                  double minPhi) {
        List<String[]> out = new ArrayList<>();
        Result r = tx.execute(
                "MATCH (a:Sample)-[r:RELATIVE {source: $source}]->(b:Sample) "
              + "WHERE r.phi >= $minPhi "
              + "RETURN a.sampleId AS a, b.sampleId AS b",
                Map.of("source", source, "minPhi", minPhi));
        try {
            while (r.hasNext()) {
                Map<String, Object> row = r.next();
                out.add(new String[]{(String) row.get("a"), (String) row.get("b")});
            }
        } finally {
            r.close();
        }
        return out;
    }

    private static Double getNullableDouble(Map<String, Object> opts, String key) {
        Object v = opts.get(key);
        return (v instanceof Number) ? ((Number) v).doubleValue() : null;
    }

    private static Long getNullableLong(Map<String, Object> opts, String key) {
        Object v = opts.get(key);
        return (v instanceof Number) ? ((Number) v).longValue() : null;
    }

    private static boolean getBoolean(Map<String, Object> opts, String key, boolean def) {
        Object v = opts.get(key);
        return (v instanceof Boolean) ? (Boolean) v : def;
    }
}
