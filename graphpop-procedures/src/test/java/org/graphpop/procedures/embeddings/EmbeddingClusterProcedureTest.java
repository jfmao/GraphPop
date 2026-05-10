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

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Cypher integration test for {@link EmbeddingClusterProcedure}. Two
 * well-separated 3-clusters in R^2 must each fall into one cluster
 * under k-means with k=2.
 */
class EmbeddingClusterProcedureTest {

    private static Neo4j embeddedNeo4j;
    private static Driver driver;

    @BeforeAll
    static void setUp() {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(EmbeddingClusterProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());
        try (Session session = driver.session()) {
            // Cluster A near (0, 0)
            createSample(session, "A1", List.of(0.0, 0.0));
            createSample(session, "A2", List.of(0.05, 0.0));
            createSample(session, "A3", List.of(0.0, 0.05));
            // Cluster B near (10, 10)
            createSample(session, "B1", List.of(10.0, 10.0));
            createSample(session, "B2", List.of(10.05, 10.0));
            createSample(session, "B3", List.of(10.0, 10.05));
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
    void kmeans_two_clusters_on_well_separated_data() {
        try (Session session = driver.session()) {
            List<Record> rs = session.run(
                "CALL graphpop.embedding.cluster('kmeans', 2, "
              + "{seed: 42}) "
              + "YIELD sample_id, cluster_id, n_clusters, method, "
              + "distance_to_centroid "
              + "RETURN sample_id, cluster_id, n_clusters, method, "
              + "distance_to_centroid").list();
            assertEquals(6, rs.size());

            Map<String, Long> assign = new HashMap<>();
            for (Record r : rs) {
                assign.put(r.get("sample_id").asString(),
                           r.get("cluster_id").asLong());
                assertEquals(2L, r.get("n_clusters").asLong());
                assertEquals("kmeans", r.get("method").asString());
                assertTrue(r.get("distance_to_centroid").asDouble() < 1.0,
                        "distance " + r.get("distance_to_centroid").asDouble()
                            + " too large for tight cluster");
            }
            // A* in one cluster, B* in the other.
            long ca = assign.get("A1");
            for (String s : new String[]{"A2", "A3"}) {
                assertEquals(ca, assign.get(s),
                        "A1 and " + s + " should share a cluster");
            }
            long cb = assign.get("B1");
            for (String s : new String[]{"B2", "B3"}) {
                assertEquals(cb, assign.get(s),
                        "B1 and " + s + " should share a cluster");
            }
            assertNotEquals(ca, cb, "A and B clusters should differ");
        }
    }

    @Test
    void kmeans_deterministic_under_same_seed() {
        try (Session session = driver.session()) {
            List<Record> r1 = session.run(
                "CALL graphpop.embedding.cluster('kmeans', 2, "
              + "{seed: 99}) "
              + "YIELD sample_id, cluster_id "
              + "RETURN sample_id, cluster_id ORDER BY sample_id").list();
            List<Record> r2 = session.run(
                "CALL graphpop.embedding.cluster('kmeans', 2, "
              + "{seed: 99}) "
              + "YIELD sample_id, cluster_id "
              + "RETURN sample_id, cluster_id ORDER BY sample_id").list();
            assertEquals(r1.size(), r2.size());
            for (int i = 0; i < r1.size(); i++) {
                assertEquals(r1.get(i).get("cluster_id").asLong(),
                        r2.get(i).get("cluster_id").asLong(),
                        "deterministic under same seed; row " + i);
            }
        }
    }

    @Test
    void kmeans_unsupported_method_throws() {
        try (Session session = driver.session()) {
            assertThrows(Exception.class, () -> session.run(
                "CALL graphpop.embedding.cluster('hdbscan', 2, {}) "
              + "YIELD sample_id RETURN sample_id").list());
        }
    }

    @Test
    void kmeans_invalid_k_throws() {
        try (Session session = driver.session()) {
            assertThrows(Exception.class, () -> session.run(
                "CALL graphpop.embedding.cluster('kmeans', 0, {}) "
              + "YIELD sample_id RETURN sample_id").list());
            assertThrows(Exception.class, () -> session.run(
                "CALL graphpop.embedding.cluster('kmeans', 100, {}) "
              + "YIELD sample_id RETURN sample_id").list());
        }
    }

    @Test
    void kmeans_k_equals_n_each_sample_own_cluster() {
        try (Session session = driver.session()) {
            List<Record> rs = session.run(
                "CALL graphpop.embedding.cluster('kmeans', 6, "
              + "{seed: 42}) "
              + "YIELD sample_id, cluster_id RETURN sample_id, cluster_id").list();
            Set<Long> clusters = new HashSet<>();
            for (Record r : rs) clusters.add(r.get("cluster_id").asLong());
            assertEquals(6, clusters.size());
        }
    }
}
