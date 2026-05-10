package org.graphpop.procedures.embeddings;

import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.neo4j.driver.Driver;
import org.neo4j.driver.GraphDatabase;
import org.neo4j.driver.Record;
import org.neo4j.driver.Session;
import org.neo4j.harness.Neo4j;
import org.neo4j.harness.Neo4jBuilders;

import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Cypher integration test for {@link EmbeddingKnnProcedure}. Hand-
 * authored 6-sample fixture with two well-separated 3-clusters in
 * R^2: kNN must rank within-cluster neighbours above cross-cluster.
 */
class EmbeddingKnnProcedureTest {

    private static Neo4j embeddedNeo4j;
    private static Driver driver;

    @BeforeAll
    static void setUp() {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(EmbeddingKnnProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());
        try (Session session = driver.session()) {
            // Cluster A: roughly along (1, 0)
            createSample(session, "A1", List.of(1.0, 0.0));
            createSample(session, "A2", List.of(0.95, 0.05));
            createSample(session, "A3", List.of(0.90, 0.10));
            // Cluster B: roughly along (0, 1)
            createSample(session, "B1", List.of(0.0, 1.0));
            createSample(session, "B2", List.of(0.05, 0.95));
            createSample(session, "B3", List.of(0.10, 0.90));
            // Sample with no embedding (must be skipped)
            session.run("CREATE (:Sample {sampleId: 'NOEMB'})");
            // Sample with zero-norm embedding (must be skipped)
            createSample(session, "ZERO", List.of(0.0, 0.0));
        }
    }

    private static void createSample(Session session, String sid,
                                      List<Double> emb) {
        session.run(
            "CREATE (:Sample {sampleId: $sid, embedding: $emb})",
            Map.of("sid", sid, "emb", emb));
    }

    @AfterAll
    static void tearDown() {
        if (driver != null) driver.close();
        if (embeddedNeo4j != null) embeddedNeo4j.close();
    }

    @Test
    void knn_returns_top_k_within_cluster_first() {
        try (Session session = driver.session()) {
            List<Record> rs = session.run(
                "CALL graphpop.embedding.knn('A1', 3) "
              + "YIELD query_sample_id, neighbor_sample_id, "
              + "cosine_similarity, rank "
              + "RETURN query_sample_id, neighbor_sample_id, "
              + "cosine_similarity, rank ORDER BY rank").list();
            assertEquals(3, rs.size());
            // Top 2 should be A2 and A3 (within-cluster).
            assertEquals("A1", rs.get(0).get("query_sample_id").asString());
            String r1 = rs.get(0).get("neighbor_sample_id").asString();
            String r2 = rs.get(1).get("neighbor_sample_id").asString();
            assertTrue(r1.startsWith("A"), "rank-1 within-cluster, got " + r1);
            assertTrue(r2.startsWith("A"), "rank-2 within-cluster, got " + r2);
            // Ranks contiguous starting at 1.
            for (int i = 0; i < 3; i++) {
                assertEquals(i + 1L, rs.get(i).get("rank").asLong());
            }
            // Similarities monotonically non-increasing.
            for (int i = 1; i < 3; i++) {
                double prev = rs.get(i - 1).get("cosine_similarity").asDouble();
                double curr = rs.get(i).get("cosine_similarity").asDouble();
                assertTrue(prev >= curr,
                        "rank " + i + " sim " + curr + " > prev " + prev);
            }
        }
    }

    @Test
    void knn_excludes_query_sample() {
        try (Session session = driver.session()) {
            List<Record> rs = session.run(
                "CALL graphpop.embedding.knn('A1', 100) "
              + "YIELD neighbor_sample_id RETURN neighbor_sample_id").list();
            for (Record r : rs) {
                assertNotEquals("A1", r.get("neighbor_sample_id").asString());
            }
        }
    }

    @Test
    void knn_skips_zero_norm_and_no_embedding_samples() {
        try (Session session = driver.session()) {
            List<Record> rs = session.run(
                "CALL graphpop.embedding.knn('A1', 100) "
              + "YIELD neighbor_sample_id RETURN neighbor_sample_id").list();
            for (Record r : rs) {
                String nid = r.get("neighbor_sample_id").asString();
                assertNotEquals("NOEMB", nid);
                assertNotEquals("ZERO", nid);
            }
        }
    }

    @Test
    void knn_unknown_sample_returns_empty() {
        try (Session session = driver.session()) {
            List<Record> rs = session.run(
                "CALL graphpop.embedding.knn('does_not_exist', 3) "
              + "YIELD neighbor_sample_id RETURN neighbor_sample_id").list();
            assertEquals(0, rs.size());
        }
    }

    @Test
    void knn_zero_k_returns_empty() {
        try (Session session = driver.session()) {
            List<Record> rs = session.run(
                "CALL graphpop.embedding.knn('A1', 0) "
              + "YIELD neighbor_sample_id RETURN neighbor_sample_id").list();
            assertEquals(0, rs.size());
        }
    }

    @Test
    void knn_cosine_values_sane() {
        try (Session session = driver.session()) {
            // A1 = (1,0), A2 = (0.95, 0.05). Cosine ≈ 0.998.
            List<Record> rs = session.run(
                "CALL graphpop.embedding.knn('A1', 1) "
              + "YIELD neighbor_sample_id, cosine_similarity "
              + "RETURN neighbor_sample_id, cosine_similarity").list();
            assertEquals("A2", rs.get(0).get("neighbor_sample_id").asString());
            double sim = rs.get(0).get("cosine_similarity").asDouble();
            // 0.95 / sqrt(0.95^2 + 0.05^2) ≈ 0.99861
            assertEquals(0.95 / Math.sqrt(0.9025 + 0.0025), sim, 1e-9);
        }
    }
}
