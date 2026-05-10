package org.graphpop.procedures.pairwise;

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
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Cypher integration test for {@link RelateClassifyProcedure}.
 *
 * <p>Hand-authored {@code :KINSHIP {method:'king-robust'}} edges
 * cover every Manichaikul 2010 relationship category; the classifier
 * must label each correctly and persist {@code :RELATIVE} edges.</p>
 */
class RelateClassifyProcedureTest {

    private static Neo4j embeddedNeo4j;
    private static Driver driver;

    @BeforeAll
    static void setUp() {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(RelateClassifyProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());

        try (Session session = driver.session()) {
            for (int i = 0; i < 10; i++) {
                session.run("CREATE (:Sample {sampleId: $sid})",
                            Map.of("sid", "S" + i));
            }
            // Hand-authored KING-robust edges covering every category.
            createKing(session, "S0", "S1", 0.45, 0.0);    // identical
            createKing(session, "S2", "S3", 0.25, 0.001);  // parent-child
            createKing(session, "S4", "S5", 0.25, 0.06);   // full-sibling
            createKing(session, "S6", "S7", 0.125, 0.05);  // 2nd-degree
            createKing(session, "S0", "S8", 0.0625, 0.05); // 3rd-degree
            createKing(session, "S0", "S9", 0.02, 0.10);   // unrelated

            // IBD fixture: 4 pairs on 4 different chromosomes so the
            // genome span aggregates to 91.5 Mb (denom = 183 Mb), and
            // each pair's total IBD length lands in a known phi bucket.
            for (int i = 0; i < 8; i++) {
                session.run("CREATE (:Sample {sampleId: $sid})",
                            Map.of("sid", "T" + i));
            }
            createIbd(session, "T0", "T1", "1", 0L, 50_000_000L, "test_ibd");
            createIbd(session, "T2", "T3", "2", 0L, 25_000_000L, "test_ibd");
            createIbd(session, "T4", "T5", "3", 0L, 12_500_000L, "test_ibd");
            createIbd(session, "T6", "T7", "4", 0L,  4_000_000L, "test_ibd");
        }
    }

    private static void createIbd(Session session, String a, String b,
                                   String chr, long start, long end,
                                   String source) {
        session.run(
            "MATCH (a:Sample {sampleId: $a}), (b:Sample {sampleId: $b}) "
          + "CREATE (a)-[:IBD_SEGMENT {chr: $chr, start: $start, "
          + "end: $end, length_bp: $length_bp, source: $source}]->(b)",
            Map.of("a", a, "b", b, "chr", chr, "start", start, "end", end,
                   "length_bp", end - start, "source", source));
    }

    private static void createKing(Session session, String a, String b,
                                    double phi, double ibs0Frac) {
        session.run(
            "MATCH (a:Sample {sampleId: $a}), (b:Sample {sampleId: $b}) "
          + "CREATE (a)-[:KINSHIP {method: 'king-robust', phi: $phi, "
          + "ibs0_frac: $ibs0_frac}]->(b)",
            Map.of("a", a, "b", b, "phi", phi, "ibs0_frac", ibs0Frac));
    }

    @AfterAll
    static void tearDown() {
        if (driver != null) driver.close();
        if (embeddedNeo4j != null) embeddedNeo4j.close();
    }

    @Test
    void classify_king_full_category_matrix() {
        try (Session session = driver.session()) {
            // min_phi=0 so even unrelated pair shows up.
            List<Record> rows = session.run(
                "CALL graphpop.relate.classify('king', {min_phi: 0.0}) "
              + "YIELD sample_a, sample_b, relationship, degree, phi, ibs0_frac "
              + "RETURN sample_a, sample_b, relationship, degree "
              + "ORDER BY sample_a, sample_b").list();
            assertEquals(6, rows.size());

            Map<String, String> byPair = new HashMap<>();
            for (Record r : rows) {
                String key = r.get("sample_a").asString() + "," + r.get("sample_b").asString();
                byPair.put(key, r.get("relationship").asString());
            }
            assertEquals("identical", byPair.get("S0,S1"));
            assertEquals("parent_child", byPair.get("S2,S3"));
            assertEquals("full_sibling", byPair.get("S4,S5"));
            assertEquals("second_degree", byPair.get("S6,S7"));
            assertEquals("third_degree", byPair.get("S0,S8"));
            assertEquals("unrelated", byPair.get("S0,S9"));
        }
    }

    @Test
    void classify_default_min_phi_drops_unrelated() {
        try (Session session = driver.session()) {
            // Default min_phi = 0.0442; unrelated (phi=0.02) must drop.
            long c = session.run(
                "CALL graphpop.relate.classify('king', {persist: false}) "
              + "YIELD relationship "
              + "WHERE relationship = 'unrelated' RETURN count(*) AS c"
            ).single().get("c").asLong();
            assertEquals(0L, c);
        }
    }

    @Test
    void classify_persists_relative_edges() {
        try (Session session = driver.session()) {
            session.run(
                "CALL graphpop.relate.classify('king', {min_phi: 0.0}) "
              + "YIELD sample_a RETURN count(*) AS c").single();
            long edgeCount = session.run(
                "MATCH ()-[r:RELATIVE {source: 'king'}]->() RETURN count(r) AS c"
            ).single().get("c").asLong();
            assertEquals(6L, edgeCount);
        }
    }

    @Test
    void classify_idempotent_re_run() {
        try (Session session = driver.session()) {
            session.run(
                "CALL graphpop.relate.classify('king', {min_phi: 0.0}) "
              + "YIELD sample_a RETURN count(*) AS c").single();
            session.run(
                "CALL graphpop.relate.classify('king', {min_phi: 0.0}) "
              + "YIELD sample_a RETURN count(*) AS c").single();
            long edgeCount = session.run(
                "MATCH ()-[r:RELATIVE {source: 'king'}]->() RETURN count(r) AS c"
            ).single().get("c").asLong();
            assertEquals(6L, edgeCount, "re-run must not duplicate");
        }
    }

    @Test
    void classify_custom_cutoffs_override_defaults() {
        try (Session session = driver.session()) {
            // Lower the identical bar to 0.10; now (S0,S1)/(S2,S3)/(S4,S5)
            // all qualify as "identical".
            List<Record> rows = session.run(
                "CALL graphpop.relate.classify('king', {min_phi: 0.0, "
              + "identical: 0.10, first_degree: 0.05, "
              + "second_degree: 0.025, third_degree: 0.01, "
              + "persist: false}) "
              + "YIELD sample_a, sample_b, relationship "
              + "WHERE relationship = 'identical' "
              + "RETURN sample_a, sample_b").list();
            assertTrue(rows.size() >= 3,
                "custom thresholds must apply; got " + rows.size() + " identical pairs");
        }
    }

    @Test
    void classify_ibd_source_phi_to_relationship() {
        try (Session session = driver.session()) {
            // Genome span = 50M + 25M + 12.5M + 4M = 91.5 Mb across 4 chrs;
            // denom = 183 Mb. Per-pair phi:
            //   T0-T1: 50/183 = 0.273  → first_degree (NaN ibs0)
            //   T2-T3: 25/183 = 0.137  → second_degree
            //   T4-T5: 12.5/183 = 0.068 → third_degree
            //   T6-T7: 4/183 = 0.022   → unrelated
            List<Record> rows = session.run(
                "CALL graphpop.relate.classify('test_ibd', {min_phi: 0.0}) "
              + "YIELD sample_a, sample_b, relationship, phi, ibs0_frac "
              + "RETURN sample_a, sample_b, relationship, phi, ibs0_frac "
              + "ORDER BY sample_a, sample_b").list();
            assertEquals(4, rows.size());

            Map<String, String> byPair = new HashMap<>();
            for (Record r : rows) {
                String key = r.get("sample_a").asString() + ","
                           + r.get("sample_b").asString();
                byPair.put(key, r.get("relationship").asString());
                // IBD source must report NaN ibs0_frac (no IBS info).
                assertTrue(Double.isNaN(r.get("ibs0_frac").asDouble()),
                        "IBD source must report NaN ibs0_frac for " + key);
            }
            assertEquals("first_degree", byPair.get("T0,T1"),
                "IBD source cannot disambiguate PC vs FS → first_degree");
            assertEquals("second_degree", byPair.get("T2,T3"));
            assertEquals("third_degree", byPair.get("T4,T5"));
            assertEquals("unrelated", byPair.get("T6,T7"));

            // Persisted edges tagged with the IBD source.
            long edgeCount = session.run(
                "MATCH ()-[r:RELATIVE {source: 'test_ibd'}]->() "
              + "RETURN count(r) AS c"
            ).single().get("c").asLong();
            assertEquals(4L, edgeCount);
        }
    }

    @Test
    void classify_unknown_source_returns_empty() {
        try (Session session = driver.session()) {
            long c = session.run(
                "CALL graphpop.relate.classify('does_not_exist') "
              + "YIELD sample_a RETURN count(*) AS c"
            ).single().get("c").asLong();
            assertEquals(0L, c);
        }
    }
}
