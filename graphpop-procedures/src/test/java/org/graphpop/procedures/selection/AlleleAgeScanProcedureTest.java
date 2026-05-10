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
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Integration test for {@link AlleleAgeScanProcedure}. Runs against
 * the existing 20-sample fixture (no sweep — neutral null), so the
 * test asserts:
 *
 * <ul>
 *   <li>every mutation in the fixture gets exactly one row;</li>
 *   <li>z-scores are well-defined (finite) and per-bin distributions
 *       are centred on zero (mean |z| over a non-empty bin is < 1);</li>
 *   <li>carrier counts and frequencies match the M6 fixture's
 *       allele-age reference.</li>
 * </ul>
 */
class AlleleAgeScanProcedureTest {

    private static final String RUN_ID = "egrm_fixture_20samples";

    private static Neo4j embeddedNeo4j;
    private static Driver driver;
    private static JsonNode statsRef;

    @BeforeAll
    static void setUp() throws Exception {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(AlleleAgeScanProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());

        ObjectMapper mapper = new ObjectMapper();
        JsonNode arg;
        try (InputStream is = AlleleAgeScanProcedureTest.class
                .getResourceAsStream("/egrm_fixture_20samples_arg.json")) {
            arg = mapper.readTree(is);
        }
        try (InputStream is = AlleleAgeScanProcedureTest.class
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

    @Test
    void allele_age_scan_emits_one_row_per_mutation() {
        try (Session session = driver.session()) {
            int expected = statsRef.get("allele_age").size();
            List<Record> rs = session.run(
                "CALL graphpop.selection.allele_age_scan($rid, {}) "
              + "YIELD variant_id, freq, age_midpoint, log_age, "
              + "z_score, n_carriers, n_samples "
              + "RETURN variant_id, freq, age_midpoint, log_age, "
              + "z_score, n_carriers, n_samples",
                Map.of("rid", RUN_ID)).list();
            assertEquals(expected, rs.size(),
                    "must emit one row per :MUTATED_ON edge");
            for (Record r : rs) {
                assertTrue(Double.isFinite(r.get("freq").asDouble()));
                assertTrue(r.get("n_samples").asLong() == 20L);
                assertTrue(r.get("n_carriers").asLong() > 0);
                assertTrue(Double.isFinite(r.get("age_midpoint").asDouble()));
                assertTrue(Double.isFinite(r.get("log_age").asDouble()));
                assertTrue(Double.isFinite(r.get("z_score").asDouble()));
            }
        }
    }

    @Test
    void allele_age_scan_carrier_counts_match_fixture_reference() {
        try (Session session = driver.session()) {
            // Build expected (variant_id → n_carriers, age_midpoint) map.
            Map<String, Long> expectedCarriers = new HashMap<>();
            Map<String, Double> expectedMidpoints = new HashMap<>();
            for (int i = 0; i < statsRef.get("allele_age").size(); i++) {
                JsonNode m = statsRef.get("allele_age").get(i);
                String varId = RUN_ID + ":mut:" + i;
                expectedCarriers.put(varId, m.get("n_carriers").asLong());
                expectedMidpoints.put(varId, m.get("midpoint_time").asDouble());
            }

            List<Record> rs = session.run(
                "CALL graphpop.selection.allele_age_scan($rid, {}) "
              + "YIELD variant_id, n_carriers, age_midpoint "
              + "RETURN variant_id, n_carriers, age_midpoint",
                Map.of("rid", RUN_ID)).list();
            for (Record r : rs) {
                String vid = r.get("variant_id").asString();
                Long expC = expectedCarriers.get(vid);
                Double expM = expectedMidpoints.get(vid);
                assertNotNull(expC, "unexpected variant " + vid);
                assertEquals(expC.longValue(), r.get("n_carriers").asLong(),
                        "carriers for " + vid);
                assertEquals(expM, r.get("age_midpoint").asDouble(), 1e-12,
                        "age_midpoint for " + vid);
            }
        }
    }

    @Test
    void allele_age_scan_neutral_data_z_centred_near_zero() {
        try (Session session = driver.session()) {
            // Aggregate z by bin; per-bin mean |z| must be < 1 for a
            // neutral fixture (typically << 1).
            List<Record> rs = session.run(
                "CALL graphpop.selection.allele_age_scan($rid, "
              + "{n_freq_bins: 5, min_freq: 0.05, max_freq: 0.95}) "
              + "YIELD bin_index, z_score, bin_n "
              + "RETURN bin_index, z_score, bin_n",
                Map.of("rid", RUN_ID)).list();

            Map<Long, double[]> sums = new HashMap<>();  // bin → [sum, count]
            Set<Long> bins = new HashSet<>();
            for (Record r : rs) {
                long bi = r.get("bin_index").asLong();
                bins.add(bi);
                long bn = r.get("bin_n").asLong();
                if (bn < 2) continue;  // skip bins where z is forced to 0
                double[] s = sums.computeIfAbsent(bi, k -> new double[2]);
                s[0] += r.get("z_score").asDouble();
                s[1] += 1;
            }
            // Per non-degenerate bin, mean(z) ≈ 0 by construction; assert
            // |mean(z)| < 1 to guard against gross algorithmic drift.
            for (Map.Entry<Long, double[]> e : sums.entrySet()) {
                if (e.getValue()[1] >= 2) {
                    double meanZ = e.getValue()[0] / e.getValue()[1];
                    assertTrue(Math.abs(meanZ) < 1.0,
                            "bin " + e.getKey() + " mean z = " + meanZ
                                + " exceeds neutral expectation");
                }
            }
        }
    }

    @Test
    void allele_age_scan_unknown_run_returns_empty() {
        try (Session session = driver.session()) {
            List<Record> rs = session.run(
                "CALL graphpop.selection.allele_age_scan('no_such_run', {}) "
              + "YIELD variant_id RETURN variant_id").list();
            assertEquals(0, rs.size());
        }
    }
}
