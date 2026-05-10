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
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Integration test for {@link ArgAlleleAgeProcedure}. Validates per-
 * mutation time bracket + carrier count against the reference dumped
 * by {@code build_egrm_fixture.py}.
 */
class ArgAlleleAgeProcedureTest {

    private static final String RUN_ID = "egrm_fixture_20samples";

    private static Neo4j embeddedNeo4j;
    private static Driver driver;
    private static JsonNode statsRef;

    @BeforeAll
    static void setUp() throws Exception {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(ArgAlleleAgeProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());

        ObjectMapper mapper = new ObjectMapper();
        JsonNode arg;
        try (InputStream is = ArgAlleleAgeProcedureTest.class
                .getResourceAsStream("/egrm_fixture_20samples_arg.json")) {
            arg = mapper.readTree(is);
        }
        try (InputStream is = ArgAlleleAgeProcedureTest.class
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
    void allele_age_matches_tskit_for_every_mutation() {
        try (Session session = driver.session()) {
            JsonNode allele = statsRef.get("allele_age");
            for (int i = 0; i < allele.size(); i++) {
                JsonNode m = allele.get(i);
                String varId = RUN_ID + ":mut:" + i;
                List<Record> rs = session.run(
                    "CALL graphpop.arg.allele_age($rid, $vid) "
                  + "YIELD child_node_id, parent_node_id, child_time, "
                  + "parent_time, midpoint_time, n_carriers "
                  + "RETURN child_node_id, parent_node_id, child_time, "
                  + "parent_time, midpoint_time, n_carriers",
                    Map.of("rid", RUN_ID, "vid", varId)).list();
                assertEquals(1, rs.size(),
                        "expected 1 row for mutation " + i);
                Record r = rs.get(0);
                assertEquals(m.get("child_node_id").asLong(),
                        r.get("child_node_id").asLong(),
                        "child_node_id for mutation " + i);
                assertEquals(m.get("parent_node_id").asLong(),
                        r.get("parent_node_id").asLong(),
                        "parent_node_id for mutation " + i);
                assertEquals(m.get("child_time").asDouble(),
                        r.get("child_time").asDouble(), 1e-12,
                        "child_time for mutation " + i);
                assertEquals(m.get("parent_time").asDouble(),
                        r.get("parent_time").asDouble(), 1e-12,
                        "parent_time for mutation " + i);
                assertEquals(m.get("midpoint_time").asDouble(),
                        r.get("midpoint_time").asDouble(), 1e-12,
                        "midpoint for mutation " + i);
                assertEquals(m.get("n_carriers").asLong(),
                        r.get("n_carriers").asLong(),
                        "n_carriers for mutation " + i);
            }
        }
    }

    @Test
    void allele_age_unknown_variant_returns_empty() {
        try (Session session = driver.session()) {
            List<Record> rs = session.run(
                "CALL graphpop.arg.allele_age($rid, 'no_such_variant') "
              + "YIELD variant_id RETURN variant_id",
                Map.of("rid", RUN_ID)).list();
            assertEquals(0, rs.size());
        }
    }
}
