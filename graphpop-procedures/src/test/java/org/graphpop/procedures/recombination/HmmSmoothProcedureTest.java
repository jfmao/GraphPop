package org.graphpop.procedures.recombination;

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

class HmmSmoothProcedureTest {

    private static Neo4j embeddedNeo4j;
    private static Driver driver;

    @BeforeAll
    static void setUp() {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(HmmSmoothProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());
        try (Session session = driver.session()) {
            for (int i = 0; i < 10; i++) {
                session.run("CREATE (:Sample {sampleId: $sid})",
                    Map.of("sid", "S" + i));
            }
            for (long pos = 100; pos <= 5_000; pos += 500) {
                String vid = "v_" + pos;
                session.run("CREATE (:Variant {variantId: $vid, position: $pos})",
                    Map.of("vid", vid, "pos", pos));
                int nCar = 3 + ((int)(pos / 500)) % 5;
                for (int i = 0; i < nCar; i++) {
                    session.run(
                        "MATCH (s:Sample {sampleId: $sid}), "
                      + "(v:Variant {variantId: $vid}) "
                      + "CREATE (s)-[:CARRIES {gt: 1}]->(v)",
                        Map.of("sid", "S" + i, "vid", vid));
                }
            }
        }
    }

    @AfterAll
    static void tearDown() {
        if (driver != null) driver.close();
        if (embeddedNeo4j != null) embeddedNeo4j.close();
    }

    private static List<String> sids10() {
        return List.of("S0", "S1", "S2", "S3", "S4",
                       "S5", "S6", "S7", "S8", "S9");
    }

    @Test
    void hmm_smooth_emits_one_row_per_window() {
        try (Session session = driver.session()) {
            List<Record> rs = session.run(
                "CALL graphpop.recombination.hmm_smooth($sids, "
              + "{window_size: 1000, step: 1000, "
              + " max_pair_distance: 2000, n_states: 15}) "
              + "YIELD start, end, rho_per_bp, rho_smoothed, "
              + "hmm_state, n_variant_pairs "
              + "RETURN start, end, rho_per_bp, rho_smoothed, "
              + "hmm_state, n_variant_pairs ORDER BY start",
                Map.of("sids", sids10())).list();
            assertTrue(rs.size() >= 4,
                    "expected ≥ 4 windows; got " + rs.size());
            for (Record r : rs) {
                double rawRho = r.get("rho_per_bp").asDouble();
                double smoothed = r.get("rho_smoothed").asDouble();
                long state = r.get("hmm_state").asLong();
                assertTrue(rawRho >= 0.0 && Double.isFinite(rawRho));
                assertTrue(smoothed >= 0.0 && Double.isFinite(smoothed));
                assertTrue(state >= 0 && state < 15);
            }
        }
    }

    @Test
    void hmm_smooth_invalid_options_throw() {
        try (Session session = driver.session()) {
            assertThrows(Exception.class, () -> session.run(
                "CALL graphpop.recombination.hmm_smooth($sids, "
              + "{window_size: 0}) YIELD start RETURN start",
                Map.of("sids", sids10())).list());
        }
    }

    @Test
    void hmm_smooth_empty_sample_list_returns_empty() {
        try (Session session = driver.session()) {
            List<Record> rs = session.run(
                "CALL graphpop.recombination.hmm_smooth([], {}) "
              + "YIELD start RETURN start").list();
            assertEquals(0, rs.size());
        }
    }
}
