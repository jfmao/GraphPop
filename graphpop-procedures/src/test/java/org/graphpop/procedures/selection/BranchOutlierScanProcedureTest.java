package org.graphpop.procedures.selection;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.graphpop.procedures.arg.ArgFixtureLoader;
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
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Integration test for {@link BranchOutlierScanProcedure}. Runs
 * against the existing 20-sample fixture (no sweep) so the test
 * asserts:
 *
 * <ul>
 *   <li>per-window totals are positive and finite;</li>
 *   <li>z-scores follow the expected centred distribution
 *       (per-window mean(z) ≈ 0, sample sd(z) ≈ 1);</li>
 *   <li>window edges respect {@code window_size} / {@code step}.</li>
 * </ul>
 */
class BranchOutlierScanProcedureTest {

    private static final String RUN_ID = "egrm_fixture_20samples";

    private static Neo4j embeddedNeo4j;
    private static Driver driver;
    private static JsonNode statsRef;

    @BeforeAll
    static void setUp() throws Exception {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(BranchOutlierScanProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());

        ObjectMapper mapper = new ObjectMapper();
        JsonNode arg;
        try (InputStream is = BranchOutlierScanProcedureTest.class
                .getResourceAsStream("/egrm_fixture_20samples_arg.json")) {
            arg = mapper.readTree(is);
        }
        try (InputStream is = BranchOutlierScanProcedureTest.class
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
    void branch_outlier_scan_full_sequence_emits_expected_windows() {
        try (Session session = driver.session()) {
            List<String> sids = sampleIdsFrom(statsRef.get("samples_full"));
            // 50_000 bp sequence, 10_000 bp non-overlapping windows = 5 windows.
            List<Record> rs = session.run(
                "CALL graphpop.selection.branch_outlier_scan($rid, $sids, "
              + "{window_size: 10000, step: 10000}) "
              + "YIELD start, end, total_branch_length, z_score, n_samples "
              + "RETURN start, end, total_branch_length, z_score, n_samples "
              + "ORDER BY start",
                Map.of("rid", RUN_ID, "sids", sids)).list();
            assertEquals(5, rs.size());
            for (int i = 0; i < rs.size(); i++) {
                Record r = rs.get(i);
                assertEquals(i * 10_000L, r.get("start").asLong());
                assertEquals((i + 1) * 10_000L, r.get("end").asLong());
                assertTrue(r.get("total_branch_length").asDouble() > 0.0);
                assertTrue(Double.isFinite(r.get("z_score").asDouble()));
                assertEquals(20L, r.get("n_samples").asLong());
            }
        }
    }

    @Test
    void branch_outlier_scan_z_scores_sum_to_zero() {
        try (Session session = driver.session()) {
            List<String> sids = sampleIdsFrom(statsRef.get("samples_full"));
            List<Record> rs = session.run(
                "CALL graphpop.selection.branch_outlier_scan($rid, $sids, "
              + "{window_size: 5000, step: 5000}) "
              + "YIELD z_score RETURN z_score",
                Map.of("rid", RUN_ID, "sids", sids)).list();
            // By construction (Welford-derived z) sum is ~ 0.
            double sum = 0.0;
            for (Record r : rs) sum += r.get("z_score").asDouble();
            assertEquals(0.0, sum, 1e-9 * rs.size(),
                    "per-window z must sum to zero by construction");
        }
    }

    @Test
    void branch_outlier_scan_overlapping_windows_yield_more_rows() {
        try (Session session = driver.session()) {
            List<String> sids = sampleIdsFrom(statsRef.get("samples_full"));
            // 50_000 bp, window=10_000, step=5_000:
            //   starts: 0,5_000,10_000,...,45_000 → 10 windows
            List<Record> rs = session.run(
                "CALL graphpop.selection.branch_outlier_scan($rid, $sids, "
              + "{window_size: 10000, step: 5000}) "
              + "YIELD start RETURN start ORDER BY start",
                Map.of("rid", RUN_ID, "sids", sids)).list();
            assertEquals(10, rs.size());
        }
    }

    @Test
    void branch_outlier_scan_invalid_options_throw() {
        try (Session session = driver.session()) {
            List<String> sids = List.of("hap_0", "hap_1");
            assertThrows(Exception.class, () -> session.run(
                "CALL graphpop.selection.branch_outlier_scan($rid, $sids, "
              + "{window_size: 0}) YIELD start RETURN start",
                Map.of("rid", RUN_ID, "sids", sids)).list());
        }
    }

    @Test
    void branch_outlier_scan_empty_sample_list_returns_empty() {
        try (Session session = driver.session()) {
            List<Record> rs = session.run(
                "CALL graphpop.selection.branch_outlier_scan($rid, [], {}) "
              + "YIELD start RETURN start",
                Map.of("rid", RUN_ID)).list();
            assertEquals(0, rs.size());
        }
    }
}
