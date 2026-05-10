package org.graphpop.procedures.community;

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
 * Cypher integration test for {@link LouvainProcedure}. Hand-built
 * three-community :RELATIVE graph (3 disjoint K4 cliques + 3 weak
 * bridges); Louvain must recover the three communities.
 */
class LouvainProcedureTest {

    private static Neo4j embeddedNeo4j;
    private static Driver driver;

    @BeforeAll
    static void setUp() {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(LouvainProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());

        try (Session session = driver.session()) {
            for (int i = 0; i < 12; i++) {
                session.run("CREATE (:Sample {sampleId: $sid})",
                            Map.of("sid", "S" + i));
            }
            // Three disjoint K4s with weight 0.5 (1st-degree-like).
            for (int base : new int[]{0, 4, 8}) {
                for (int i = 0; i < 4; i++) {
                    for (int j = i + 1; j < 4; j++) {
                        createRel(session, "S" + (base + i),
                                  "S" + (base + j), 0.5, 1, "king");
                    }
                }
            }
            // Three weak bridge edges (2nd-degree-like).
            createRel(session, "S0", "S4",  0.05, 2, "king");
            createRel(session, "S4", "S8",  0.05, 2, "king");
            createRel(session, "S8", "S0",  0.05, 2, "king");
        }
    }

    private static void createRel(Session session, String a, String b,
                                   double phi, long degree, String source) {
        session.run(
            "MATCH (a:Sample {sampleId: $a}), (b:Sample {sampleId: $b}) "
          + "CREATE (a)-[:RELATIVE {phi: $phi, degree: $deg, "
          + "  relationship: 'rel', ibs0_frac: 0.0, source: $src, "
          + "  created_at: datetime()}]->(b)",
            Map.of("a", a, "b", b, "phi", phi, "deg", degree, "src", source));
    }

    @AfterAll
    static void tearDown() {
        if (driver != null) driver.close();
        if (embeddedNeo4j != null) embeddedNeo4j.close();
    }

    @Test
    void louvain_recovers_three_cliques_with_phi_weights() {
        try (Session session = driver.session()) {
            List<Record> rs = session.run(
                "CALL graphpop.community.louvain('king', "
              + "{edge_weight: 'phi', persist: false}) "
              + "YIELD sample_id, community_id, modularity, n_communities "
              + "RETURN sample_id, community_id, modularity, n_communities").list();
            assertEquals(12, rs.size());

            Map<String, Long> byId = new HashMap<>();
            Set<Long> communities = new HashSet<>();
            double mod = -1.0;
            long nComm = -1L;
            for (Record r : rs) {
                byId.put(r.get("sample_id").asString(),
                         r.get("community_id").asLong());
                communities.add(r.get("community_id").asLong());
                mod = r.get("modularity").asDouble();
                nComm = r.get("n_communities").asLong();
            }
            assertEquals(3, communities.size(), "must recover 3 communities");
            assertEquals(3L, nComm);
            assertTrue(mod >= 0.4, "modularity " + mod + " too low");

            // Each clique stays together.
            for (int base : new int[]{0, 4, 8}) {
                long c0 = byId.get("S" + base);
                for (int k = 1; k < 4; k++) {
                    assertEquals(c0, byId.get("S" + (base + k)),
                            "clique " + base + " split: S" + (base + k));
                }
            }
        }
    }

    @Test
    void louvain_persists_in_community_and_community_nodes() {
        try (Session session = driver.session()) {
            session.run(
                "CALL graphpop.community.louvain('king', "
              + "{edge_weight: 'phi'}) "
              + "YIELD sample_id RETURN count(*) AS c").single();

            long edges = session.run(
                "MATCH (:Sample)-[r:IN_COMMUNITY {source:'king'}]->(:Community) "
              + "RETURN count(r) AS c"
            ).single().get("c").asLong();
            assertEquals(12L, edges);

            long communities = session.run(
                "MATCH (c:Community {source:'king'}) RETURN count(c) AS c"
            ).single().get("c").asLong();
            assertEquals(3L, communities);
        }
    }

    @Test
    void louvain_idempotent_re_run() {
        try (Session session = driver.session()) {
            session.run(
                "CALL graphpop.community.louvain('king', {edge_weight: 'phi'}) "
              + "YIELD sample_id RETURN count(*) AS c").single();
            session.run(
                "CALL graphpop.community.louvain('king', {edge_weight: 'phi'}) "
              + "YIELD sample_id RETURN count(*) AS c").single();
            long edges = session.run(
                "MATCH (:Sample)-[r:IN_COMMUNITY {source:'king'}]->(:Community) "
              + "RETURN count(r) AS c"
            ).single().get("c").asLong();
            assertEquals(12L, edges, "re-run must not duplicate");
        }
    }

    @Test
    void louvain_unknown_source_returns_empty() {
        try (Session session = driver.session()) {
            List<Record> rs = session.run(
                "CALL graphpop.community.louvain('does_not_exist', "
              + "{persist: false}) "
              + "YIELD sample_id RETURN sample_id").list();
            assertEquals(0, rs.size());
        }
    }

    @Test
    void louvain_unit_weights_still_resolve_communities() {
        try (Session session = driver.session()) {
            // edge_weight='unit' = unweighted; on this fixture it should
            // still recover 3 communities (the 6 within-clique edges
            // dwarf the 1 bridge).
            List<Record> rs = session.run(
                "CALL graphpop.community.louvain('king', "
              + "{edge_weight: 'unit', persist: false}) "
              + "YIELD community_id RETURN DISTINCT community_id").list();
            assertEquals(3, rs.size());
        }
    }
}
