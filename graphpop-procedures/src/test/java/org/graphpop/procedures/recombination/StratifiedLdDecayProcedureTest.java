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

import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

class StratifiedLdDecayProcedureTest {

    private static Neo4j embeddedNeo4j;
    private static Driver driver;

    @BeforeAll
    static void setUp() {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(StratifiedLdDecayProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());
        try (Session session = driver.session()) {
            // 10 samples: 5 EUR, 5 AFR; sex: 5 M, 5 F (orthogonal).
            for (int i = 0; i < 10; i++) {
                String pop = (i < 5) ? "EUR" : "AFR";
                String sex = (i % 2 == 0) ? "M" : "F";
                session.run(
                    "CREATE (:Sample {sampleId: $sid, population: $pop, "
                  + "sex: $sex})",
                    Map.of("sid", "S" + i, "pop", pop, "sex", sex));
            }
            // Two variants in window: v1 (pos 100) carried by S0-S4 (all EUR),
            // v2 (pos 200) carried by S0-S4 → perfect r² within EUR;
            // AFR has no carriers (monomorphic) → no pairs in AFR.
            session.run("CREATE (:Variant {variantId: 'v1', position: 100})");
            session.run("CREATE (:Variant {variantId: 'v2', position: 200})");
            for (int i = 0; i < 5; i++) {
                session.run(
                    "MATCH (s:Sample {sampleId: $sid}), "
                  + "(v:Variant {variantId: 'v1'}) "
                  + "CREATE (s)-[:CARRIES {gt: 1}]->(v)",
                    Map.of("sid", "S" + i));
                session.run(
                    "MATCH (s:Sample {sampleId: $sid}), "
                  + "(v:Variant {variantId: 'v2'}) "
                  + "CREATE (s)-[:CARRIES {gt: 1}]->(v)",
                    Map.of("sid", "S" + i));
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
    void stratified_by_population_emits_per_stratum_rows() {
        try (Session session = driver.session()) {
            List<Record> rs = session.run(
                "CALL graphpop.recombination.stratified_ld_decay($sids, "
              + "{stratify_by: 'population', "
              + " window_size: 250, step: 250, "
              + " max_pair_distance: 200, min_maf: 0.1}) "
              + "YIELD stratum, n_variant_pairs, mean_r2, "
              + "n_samples, stratify_by "
              + "RETURN stratum, n_variant_pairs, mean_r2, n_samples, "
              + "stratify_by",
                Map.of("sids", sids10())).list();
            Set<String> seen = new HashSet<>();
            for (Record r : rs) {
                seen.add(r.get("stratum").asString());
                assertEquals("population", r.get("stratify_by").asString());
            }
            assertTrue(seen.contains("EUR"));
            assertTrue(seen.contains("AFR"));
        }
    }

    @Test
    void stratified_by_sex_emits_per_stratum_rows() {
        try (Session session = driver.session()) {
            List<Record> rs = session.run(
                "CALL graphpop.recombination.stratified_ld_decay($sids, "
              + "{stratify_by: 'sex', "
              + " window_size: 250, step: 250, "
              + " max_pair_distance: 200, min_maf: 0.1}) "
              + "YIELD stratum, stratify_by "
              + "RETURN stratum, stratify_by",
                Map.of("sids", sids10())).list();
            Set<String> seen = new HashSet<>();
            for (Record r : rs) {
                seen.add(r.get("stratum").asString());
                assertEquals("sex", r.get("stratify_by").asString());
            }
            assertEquals(2, seen.size());
            assertTrue(seen.contains("M"));
            assertTrue(seen.contains("F"));
        }
    }

    @Test
    void stratified_unknown_mode_throws() {
        try (Session session = driver.session()) {
            assertThrows(Exception.class, () -> session.run(
                "CALL graphpop.recombination.stratified_ld_decay($sids, "
              + "{stratify_by: 'ancestry_block'}) "
              + "YIELD start RETURN start",
                Map.of("sids", sids10())).list());
        }
    }

    @Test
    void stratified_empty_sample_list_returns_empty() {
        try (Session session = driver.session()) {
            List<Record> rs = session.run(
                "CALL graphpop.recombination.stratified_ld_decay([], {}) "
              + "YIELD start RETURN start").list();
            assertEquals(0, rs.size());
        }
    }
}
