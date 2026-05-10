package org.graphpop.procedures.arg;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.neo4j.driver.Driver;
import org.neo4j.driver.GraphDatabase;
import org.neo4j.driver.Record;
import org.neo4j.driver.Session;
import org.neo4j.harness.Neo4j;
import org.neo4j.harness.Neo4jBuilders;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Integration test for {@link ArgBranchDiversityProcedure}. Validates
 * pi-mode against {@code tskit.TreeSequence.diversity(..., mode='branch')}
 * stored in the M6 reference fixture.
 */
class ArgBranchDiversityProcedureTest {

    private static final String RUN_ID = "egrm_fixture_20samples";

    private static Neo4j embeddedNeo4j;
    private static Driver driver;
    private static JsonNode statsRef;

    @BeforeAll
    static void setUp() throws Exception {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(ArgBranchDiversityProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());

        ObjectMapper mapper = new ObjectMapper();
        JsonNode arg;
        try (InputStream is = ArgBranchDiversityProcedureTest.class
                .getResourceAsStream("/egrm_fixture_20samples_arg.json")) {
            arg = mapper.readTree(is);
        }
        try (InputStream is = ArgBranchDiversityProcedureTest.class
                .getResourceAsStream("/egrm_fixture_20samples_arg_stats.json")) {
            statsRef = mapper.readTree(is);
        }
        ArgFixtureLoader.loadInto(driver, arg, RUN_ID);
    }

    @AfterAll
    static void tearDown() {
        if (driver != null) driver.close();
        if (embeddedNeo4j != null) embeddedNeo4j.close();
    }

    private static List<String> sampleIdsFrom(JsonNode arr) {
        List<String> out = new ArrayList<>(arr.size());
        for (JsonNode n : arr) out.add("hap_" + n.asInt());
        return out;
    }

    @Test
    void branch_diversity_pi_full_matches_tskit() {
        try (Session session = driver.session()) {
            double expected = statsRef.get("branch_diversity_pi_full").asDouble();
            List<String> sids = sampleIdsFrom(statsRef.get("samples_full"));
            List<Record> rs = session.run(
                "CALL graphpop.arg.branch_diversity($rid, $sids, {mode: 'pi'}) "
              + "YIELD branch_pi, n_samples RETURN branch_pi, n_samples",
                Map.of("rid", RUN_ID, "sids", sids)).list();
            assertEquals(1, rs.size());
            assertEquals(20L, rs.get(0).get("n_samples").asLong());
            assertEquals(expected, rs.get(0).get("branch_pi").asDouble(),
                    Math.max(1e-12, 1e-9 * Math.abs(expected)),
                    "pi_branch (full) mismatch");
        }
    }

    @Test
    void branch_diversity_pi_half_matches_tskit() {
        try (Session session = driver.session()) {
            double expected = statsRef.get("branch_diversity_pi_half").asDouble();
            List<String> sids = sampleIdsFrom(statsRef.get("samples_half"));
            List<Record> rs = session.run(
                "CALL graphpop.arg.branch_diversity($rid, $sids, {}) "
              + "YIELD branch_pi RETURN branch_pi",
                Map.of("rid", RUN_ID, "sids", sids)).list();
            assertEquals(1, rs.size());
            assertEquals(expected, rs.get(0).get("branch_pi").asDouble(),
                    Math.max(1e-12, 1e-9 * Math.abs(expected)),
                    "pi_branch (10-sample subset) mismatch");
        }
    }

    @Test
    void branch_diversity_with_windows_emits_one_row_per_window() {
        try (Session session = driver.session()) {
            List<String> sids = sampleIdsFrom(statsRef.get("samples_full"));
            List<Long> wins = List.of(0L, 25_000L, 50_000L);
            List<Record> rs = session.run(
                "CALL graphpop.arg.branch_diversity($rid, $sids, "
              + "{windows: $wins}) "
              + "YIELD start, end, branch_pi RETURN start, end, branch_pi "
              + "ORDER BY start",
                Map.of("rid", RUN_ID, "sids", sids, "wins", wins)).list();
            assertEquals(2, rs.size());
            assertEquals(0L, rs.get(0).get("start").asLong());
            assertEquals(25_000L, rs.get(0).get("end").asLong());
            assertEquals(25_000L, rs.get(1).get("start").asLong());
            assertEquals(50_000L, rs.get(1).get("end").asLong());
            // Both windows must be positive and finite.
            for (Record r : rs) {
                assertTrue(Double.isFinite(r.get("branch_pi").asDouble()));
                assertTrue(r.get("branch_pi").asDouble() > 0.0);
            }
        }
    }

    @Test
    void branch_diversity_unsupported_mode_throws() {
        try (Session session = driver.session()) {
            List<String> sids = sampleIdsFrom(statsRef.get("samples_full"));
            assertThrows(Exception.class, () -> session.run(
                "CALL graphpop.arg.branch_diversity($rid, $sids, "
              + "{mode: 'theta_T'}) YIELD branch_pi RETURN branch_pi",
                Map.of("rid", RUN_ID, "sids", sids)).list());
        }
    }

    @Test
    void branch_diversity_singleton_returns_empty() {
        try (Session session = driver.session()) {
            List<Record> rs = session.run(
                "CALL graphpop.arg.branch_diversity($rid, $sids, {}) "
              + "YIELD branch_pi RETURN branch_pi",
                Map.of("rid", RUN_ID, "sids", List.of("hap_0"))).list();
            assertEquals(0, rs.size());
        }
    }
}
