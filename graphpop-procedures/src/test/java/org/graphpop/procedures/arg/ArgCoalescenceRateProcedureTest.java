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
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Integration test for {@link ArgCoalescenceRateProcedure}. Validates
 * per-bin estimator against the reference computed by
 * {@code build_egrm_fixture.py:coalescence_rate_reference}.
 */
class ArgCoalescenceRateProcedureTest {

    private static final String RUN_ID = "egrm_fixture_20samples";

    private static Neo4j embeddedNeo4j;
    private static Driver driver;
    private static JsonNode statsRef;

    @BeforeAll
    static void setUp() throws Exception {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(ArgCoalescenceRateProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());

        ObjectMapper mapper = new ObjectMapper();
        JsonNode arg;
        try (InputStream is = ArgCoalescenceRateProcedureTest.class
                .getResourceAsStream("/egrm_fixture_20samples_arg.json")) {
            arg = mapper.readTree(is);
        }
        try (InputStream is = ArgCoalescenceRateProcedureTest.class
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
    void coalescence_rate_per_bin_matches_reference() {
        try (Session session = driver.session()) {
            List<String> sids = sampleIdsFrom(statsRef.get("samples_full"));
            List<Double> bins = new ArrayList<>();
            for (JsonNode n : statsRef.get("coalescence_rate_bins")) {
                bins.add(n.asDouble());
            }
            JsonNode expected = statsRef.get("coalescence_rate_full");

            List<Record> rs = session.run(
                "CALL graphpop.arg.coalescence_rate($rid, $sids, "
              + "{time_bins: $bins}) "
              + "YIELD time_lo, time_hi, n_coalescent_events, "
              + "lineage_pair_time, rate "
              + "RETURN time_lo, time_hi, n_coalescent_events, "
              + "lineage_pair_time, rate ORDER BY time_lo",
                Map.of("rid", RUN_ID, "sids", sids, "bins", bins)).list();
            assertEquals(expected.size(), rs.size());

            for (int i = 0; i < rs.size(); i++) {
                Record r = rs.get(i);
                JsonNode e = expected.get(i);
                assertEquals(e.get("time_lo").asDouble(),
                        r.get("time_lo").asDouble(), 1e-12);
                assertEquals(e.get("time_hi").asDouble(),
                        r.get("time_hi").asDouble(), 1e-12);
                assertEquals(e.get("n_coalescent_events").asDouble(),
                        r.get("n_coalescent_events").asDouble(), 1e-9,
                        "events bin " + i);
                assertEquals(e.get("lineage_pair_time").asDouble(),
                        r.get("lineage_pair_time").asDouble(), 1e-9,
                        "pair_time bin " + i);
                double expectedRate = e.get("rate").asDouble();
                double actualRate = r.get("rate").asDouble();
                if (expectedRate == 0.0) {
                    assertEquals(0.0, actualRate, 1e-12);
                } else {
                    assertEquals(expectedRate, actualRate,
                            Math.max(1e-12, 1e-9 * Math.abs(expectedRate)),
                            "rate bin " + i);
                }
            }
        }
    }

    @Test
    void coalescence_rate_missing_time_bins_throws() {
        try (Session session = driver.session()) {
            List<String> sids = List.of("hap_0", "hap_1");
            assertThrows(Exception.class, () -> session.run(
                "CALL graphpop.arg.coalescence_rate($rid, $sids, {}) "
              + "YIELD rate RETURN rate",
                Map.of("rid", RUN_ID, "sids", sids)).list());
        }
    }

    @Test
    void coalescence_rate_non_monotonic_bins_throws() {
        try (Session session = driver.session()) {
            List<String> sids = List.of("hap_0", "hap_1");
            assertThrows(Exception.class, () -> session.run(
                "CALL graphpop.arg.coalescence_rate($rid, $sids, "
              + "{time_bins: [0.0, 1.0, 0.5]}) "
              + "YIELD rate RETURN rate",
                Map.of("rid", RUN_ID, "sids", sids)).list());
        }
    }
}
