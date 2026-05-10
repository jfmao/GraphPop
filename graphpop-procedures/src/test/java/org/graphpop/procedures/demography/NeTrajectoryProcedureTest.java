package org.graphpop.procedures.demography;

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
 * Cypher integration test for {@link NeTrajectoryProcedure}.
 *
 * <p>Validates the closed-form Ne(t) inversion of the M6 coalescence-
 * rate output. The expected per-bin Ne and SE are derived in-test
 * from the existing {@code coalescence_rate_full} reference (which
 * is itself dumped by the Python builder and validated against
 * tskit), so no fixture regeneration is needed.</p>
 */
class NeTrajectoryProcedureTest {

    private static final String RUN_ID = "egrm_fixture_20samples";

    private static Neo4j embeddedNeo4j;
    private static Driver driver;
    private static JsonNode statsRef;

    @BeforeAll
    static void setUp() throws Exception {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(NeTrajectoryProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());

        ObjectMapper mapper = new ObjectMapper();
        JsonNode arg;
        try (InputStream is = NeTrajectoryProcedureTest.class
                .getResourceAsStream("/egrm_fixture_20samples_arg.json")) {
            arg = mapper.readTree(is);
        }
        try (InputStream is = NeTrajectoryProcedureTest.class
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
    void ne_trajectory_matches_closed_form_inversion_with_default_ploidy() {
        try (Session session = driver.session()) {
            List<String> sids = sampleIdsFrom(statsRef.get("samples_full"));
            List<Double> bins = new ArrayList<>();
            for (JsonNode n : statsRef.get("coalescence_rate_bins")) {
                bins.add(n.asDouble());
            }
            JsonNode coalRef = statsRef.get("coalescence_rate_full");

            List<Record> rs = session.run(
                "CALL graphpop.demography.ne_trajectory($rid, $sids, "
              + "{time_bins: $bins}) "
              + "YIELD time_lo, time_hi, n_coalescent_events, "
              + "lineage_pair_time, rate, ne, ne_se, flag "
              + "RETURN time_lo, time_hi, n_coalescent_events, "
              + "lineage_pair_time, rate, ne, ne_se, flag "
              + "ORDER BY time_lo",
                Map.of("rid", RUN_ID, "sids", sids, "bins", bins)).list();
            assertEquals(coalRef.size(), rs.size());

            for (int i = 0; i < rs.size(); i++) {
                Record r = rs.get(i);
                JsonNode e = coalRef.get(i);
                double rate = e.get("rate").asDouble();
                double pairTime = e.get("lineage_pair_time").asDouble();

                // Underlying coalescence-rate values must round-trip.
                assertEquals(rate, r.get("rate").asDouble(), 1e-9,
                        "rate mismatch bin " + i);
                assertEquals(pairTime,
                        r.get("lineage_pair_time").asDouble(), 1e-9,
                        "pair_time mismatch bin " + i);

                if (rate > 0.0 && pairTime > 0.0) {
                    double expectedNe = 1.0 / (2.0 * rate);  // ploidy=2
                    double expectedSe = (1.0 / (2.0 * rate * rate))
                                       * Math.sqrt(rate / pairTime);
                    assertEquals(expectedNe, r.get("ne").asDouble(),
                            Math.max(1e-9, 1e-9 * Math.abs(expectedNe)),
                            "Ne mismatch bin " + i);
                    assertEquals(expectedSe, r.get("ne_se").asDouble(),
                            Math.max(1e-9, 1e-9 * Math.abs(expectedSe)),
                            "Ne SE mismatch bin " + i);
                    assertEquals("ok", r.get("flag").asString(),
                            "flag should be 'ok' for non-empty bin " + i);
                } else {
                    assertEquals(Double.POSITIVE_INFINITY,
                            r.get("ne").asDouble(),
                            "empty bin " + i + " should have Ne=+Inf");
                    assertTrue(Double.isNaN(r.get("ne_se").asDouble()),
                            "empty bin " + i + " should have NaN SE");
                    assertEquals("no_events", r.get("flag").asString());
                }
            }
        }
    }

    @Test
    void ne_trajectory_haploid_is_double_diploid() {
        try (Session session = driver.session()) {
            List<String> sids = sampleIdsFrom(statsRef.get("samples_full"));
            List<Double> bins = new ArrayList<>();
            for (JsonNode n : statsRef.get("coalescence_rate_bins")) {
                bins.add(n.asDouble());
            }

            List<Record> diploid = session.run(
                "CALL graphpop.demography.ne_trajectory($rid, $sids, "
              + "{time_bins: $bins, ploidy: 2}) "
              + "YIELD time_lo, ne RETURN time_lo, ne ORDER BY time_lo",
                Map.of("rid", RUN_ID, "sids", sids, "bins", bins)).list();
            List<Record> haploid = session.run(
                "CALL graphpop.demography.ne_trajectory($rid, $sids, "
              + "{time_bins: $bins, ploidy: 1}) "
              + "YIELD time_lo, ne RETURN time_lo, ne ORDER BY time_lo",
                Map.of("rid", RUN_ID, "sids", sids, "bins", bins)).list();
            assertEquals(diploid.size(), haploid.size());

            for (int i = 0; i < diploid.size(); i++) {
                double dNe = diploid.get(i).get("ne").asDouble();
                double hNe = haploid.get(i).get("ne").asDouble();
                if (Double.isInfinite(dNe)) {
                    assertEquals(Double.POSITIVE_INFINITY, hNe);
                } else {
                    assertEquals(2.0 * dNe, hNe,
                            Math.max(1e-9, 1e-9 * Math.abs(dNe)),
                            "Haploid Ne should be 2× diploid Ne (bin " + i + ")");
                }
            }
        }
    }

    @Test
    void ne_trajectory_zero_event_bin_returns_no_events_flag() {
        try (Session session = driver.session()) {
            List<String> sids = sampleIdsFrom(statsRef.get("samples_full"));
            // Place a tiny bin in a regime far past the deepest internal-
            // node time on the 20-sample fixture (deepest ~3-4 generations
            // for the seed=42 sim). [10, 20] generations contains zero
            // events.
            List<Double> bins = List.of(0.0, 10.0, 20.0);
            List<Record> rs = session.run(
                "CALL graphpop.demography.ne_trajectory($rid, $sids, "
              + "{time_bins: $bins}) "
              + "YIELD time_lo, ne, ne_se, flag, n_coalescent_events "
              + "RETURN time_lo, ne, ne_se, flag, n_coalescent_events "
              + "ORDER BY time_lo",
                Map.of("rid", RUN_ID, "sids", sids, "bins", bins)).list();
            assertEquals(2, rs.size());
            // Bin 1 (10–20) must be empty.
            Record empty = rs.get(1);
            assertEquals(0.0, empty.get("n_coalescent_events").asDouble(), 1e-12);
            assertEquals(Double.POSITIVE_INFINITY, empty.get("ne").asDouble());
            assertTrue(Double.isNaN(empty.get("ne_se").asDouble()));
            assertEquals("no_events", empty.get("flag").asString());
        }
    }

    @Test
    void ne_trajectory_missing_time_bins_throws() {
        try (Session session = driver.session()) {
            List<String> sids = List.of("hap_0", "hap_1");
            assertThrows(Exception.class, () -> session.run(
                "CALL graphpop.demography.ne_trajectory($rid, $sids, {}) "
              + "YIELD ne RETURN ne",
                Map.of("rid", RUN_ID, "sids", sids)).list());
        }
    }

    @Test
    void ne_trajectory_invalid_ploidy_throws() {
        try (Session session = driver.session()) {
            List<String> sids = List.of("hap_0", "hap_1");
            assertThrows(Exception.class, () -> session.run(
                "CALL graphpop.demography.ne_trajectory($rid, $sids, "
              + "{time_bins: [0, 1], ploidy: 0}) "
              + "YIELD ne RETURN ne",
                Map.of("rid", RUN_ID, "sids", sids)).list());
        }
    }

    @Test
    void ne_trajectory_unknown_population_returns_empty() {
        try (Session session = driver.session()) {
            List<Record> rs = session.run(
                "CALL graphpop.demography.ne_trajectory($rid, [], "
              + "{time_bins: [0, 1], population: 'no_such_pop'}) "
              + "YIELD ne RETURN ne",
                Map.of("rid", RUN_ID)).list();
            assertEquals(0, rs.size());
        }
    }
}
