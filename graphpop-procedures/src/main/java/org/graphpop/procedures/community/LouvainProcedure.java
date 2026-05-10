package org.graphpop.procedures.community;

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
 * Louvain modularity community detection on the relatedness graph
 * (M9, deferred from M5).
 *
 * <pre>
 * CALL graphpop.community.louvain('king',
 *         {edge_weight: 'phi', persist: true})
 *   YIELD sample_id, community_id, modularity, n_communities, source
 * </pre>
 *
 * <p>Reads {@code (:Sample)-[:RELATIVE {source}]->(:Sample)} edges,
 * runs Blondel et al. 2008 Louvain on the resulting weighted graph,
 * and emits one row per sample. With {@code persist: true} (the
 * default) also writes
 * {@code (:Sample)-[:IN_COMMUNITY {source}]->(:Community {community_id, source})}
 * edges idempotently per source.</p>
 *
 * <p>Edge-weight options:</p>
 *
 * <ul>
 *   <li>{@code edge_weight: 'phi'}    — use {@code :RELATIVE.phi}</li>
 *   <li>{@code edge_weight: 'unit'}   — every edge weight 1.0 (default)</li>
 *   <li>{@code edge_weight: 'degree'} — use the inverse of the
 *       categorical degree (closer relatives weigh more)</li>
 * </ul>
 *
 * <p>Pure-Java implementation; no GDS plugin required.</p>
 */
public class LouvainProcedure {

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.community.louvain", mode = Mode.WRITE)
    @Description("Louvain modularity community detection on :RELATIVE "
            + "edges. Pure-Java; no GDS plugin needed. Persists "
            + ":IN_COMMUNITY edges by default.")
    @SuppressWarnings("unchecked")
    public Stream<LouvainResult> louvain(
            @Name("source") String source,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        if (options == null) options = new HashMap<>();
        boolean persist = (Boolean) options.getOrDefault("persist", true);
        String weightMode = (String) options.getOrDefault("edge_weight", "unit");
        long seed = ((Number) options.getOrDefault("seed", 42L)).longValue();

        // Pull sample IDs and edges in one Cypher pass each.
        Map<String, Integer> idIndex = new HashMap<>();
        List<String> idList = new ArrayList<>();
        Result r = tx.execute(
                "MATCH (a:Sample)-[r:RELATIVE {source: $src}]->(b:Sample) "
              + "WITH collect(DISTINCT a) + collect(DISTINCT b) AS samples "
              + "UNWIND samples AS s "
              + "RETURN DISTINCT s.sampleId AS sampleId",
                Map.of("src", source));
        try {
            while (r.hasNext()) {
                String sid = (String) r.next().get("sampleId");
                if (!idIndex.containsKey(sid)) {
                    idIndex.put(sid, idList.size());
                    idList.add(sid);
                }
            }
        } finally {
            r.close();
        }
        int n = idList.size();
        if (n == 0) return Stream.empty();

        List<int[]> edgePairs = new ArrayList<>();
        List<Double> edgeWeights = new ArrayList<>();
        Result er = tx.execute(
                "MATCH (a:Sample)-[r:RELATIVE {source: $src}]->(b:Sample) "
              + "RETURN a.sampleId AS sa, b.sampleId AS sb, "
              + "r.phi AS phi, r.degree AS deg",
                Map.of("src", source));
        try {
            while (er.hasNext()) {
                Map<String, Object> row = er.next();
                Integer ai = idIndex.get(row.get("sa"));
                Integer bi = idIndex.get(row.get("sb"));
                if (ai == null || bi == null || ai.intValue() == bi.intValue()) {
                    continue;
                }
                double w = edgeWeight(row, weightMode);
                if (w <= 0) continue;
                edgePairs.add(new int[]{ai, bi});
                edgeWeights.add(w);
            }
        } finally {
            er.close();
        }
        int m = edgePairs.size();
        int[] uA = new int[m];
        int[] vA = new int[m];
        double[] wA = new double[m];
        for (int i = 0; i < m; i++) {
            uA[i] = edgePairs.get(i)[0];
            vA[i] = edgePairs.get(i)[1];
            wA[i] = edgeWeights.get(i);
        }

        LouvainComputer.Result lou = LouvainComputer.run(
                n, uA, vA, wA, 20, 50, 1e-9, seed);

        long now = System.currentTimeMillis();
        if (persist) {
            // Idempotent: clear prior :IN_FAMILY edges and Community nodes
            // for this source before writing the new partition.
            tx.execute(
                "MATCH (:Sample)-[r:IN_COMMUNITY {source: $src}]->(:Community) "
              + "DELETE r",
                Map.of("src", source));
            tx.execute(
                "MATCH (c:Community {source: $src}) "
              + "WHERE NOT (:Sample)-[:IN_COMMUNITY]->(c) DELETE c",
                Map.of("src", source));

            List<Map<String, Object>> communityRows = new ArrayList<>();
            for (int c = 0; c < lou.nCommunities; c++) {
                Map<String, Object> row = new HashMap<>();
                row.put("community_id", (long) c);
                communityRows.add(row);
            }
            tx.execute(
                "UNWIND $rows AS r "
              + "MERGE (c:Community {community_id: r.community_id, "
              + "                    source: $src}) "
              + "ON CREATE SET c.created_at = datetime({epochMillis: $now}), "
              + "              c.modularity = $mod",
                Map.of("rows", communityRows, "src", source,
                       "now", now, "mod", lou.modularity));

            List<Map<String, Object>> edgeRows = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                Map<String, Object> row = new HashMap<>();
                row.put("sample_id", idList.get(i));
                row.put("community_id", (long) lou.community[i]);
                edgeRows.add(row);
            }
            tx.execute(
                "UNWIND $rows AS r "
              + "MATCH (s:Sample {sampleId: r.sample_id}), "
              + "(c:Community {community_id: r.community_id, source: $src}) "
              + "CREATE (s)-[:IN_COMMUNITY {source: $src, "
              + "  created_at: datetime({epochMillis: $now})}]->(c)",
                Map.of("rows", edgeRows, "src", source, "now", now));
        }

        List<LouvainResult> rows = new ArrayList<>(n);
        for (int i = 0; i < n; i++) {
            rows.add(new LouvainResult(
                    idList.get(i), lou.community[i],
                    lou.modularity, lou.nCommunities, source));
        }
        return rows.stream();
    }

    private static double edgeWeight(Map<String, Object> row, String mode) {
        switch (mode) {
            case "phi":
                Object phi = row.get("phi");
                if (phi instanceof Number) {
                    double v = ((Number) phi).doubleValue();
                    return Math.max(v, 0.0);
                }
                return 0.0;
            case "degree":
                Object deg = row.get("deg");
                if (deg instanceof Number) {
                    double d = ((Number) deg).doubleValue();
                    if (d <= 0) return 0.0;
                    return 1.0 / (1.0 + d);  // 1st-degree → 0.5, 2nd → 0.33, …
                }
                return 0.0;
            case "unit":
            default:
                return 1.0;
        }
    }
}
