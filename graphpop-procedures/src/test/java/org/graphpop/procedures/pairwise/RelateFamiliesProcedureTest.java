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
 * Cypher integration test for {@link RelateFamiliesProcedure}.
 *
 * <p>Hand-authored {@code :RELATIVE} and {@code :IBD_SEGMENT} edges
 * cover three predicate modes: degree, IBD-bp, and phi.</p>
 */
class RelateFamiliesProcedureTest {

    private static Neo4j embeddedNeo4j;
    private static Driver driver;

    @BeforeAll
    static void setUp() {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(RelateFamiliesProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());

        try (Session session = driver.session()) {
            // 12 samples: three trios + three singletons.
            for (int i = 0; i < 12; i++) {
                session.run("CREATE (:Sample {sampleId: $sid})",
                            Map.of("sid", "S" + i));
            }
            // Trio 1: S0-S1-S2 via first-degree :RELATIVE edges.
            createRel(session, "S0", "S1", 1, 0.25, "king");
            createRel(session, "S1", "S2", 1, 0.25, "king");
            // Trio 2: S3-S4-S5 via mixed first-/second-degree edges.
            createRel(session, "S3", "S4", 1, 0.25, "king");
            createRel(session, "S4", "S5", 2, 0.125, "king");
            // Trio 3: S6-S7-S8 via second-degree edges only.
            createRel(session, "S6", "S7", 2, 0.125, "king");
            createRel(session, "S7", "S8", 2, 0.125, "king");
            // S9, S10, S11: singletons (no edges).

            // IBD-bp predicate fixture (separate sample range).
            for (int i = 0; i < 4; i++) {
                session.run("CREATE (:Sample {sampleId: $sid})",
                            Map.of("sid", "T" + i));
            }
            // T0-T1: 8 Mb IBD (above 5 Mb threshold)
            // T2-T3: 3 Mb IBD (below threshold)
            createIbd(session, "T0", "T1", 8_000_000L, "test_ibd");
            createIbd(session, "T2", "T3", 3_000_000L, "test_ibd");
        }
    }

    private static void createRel(Session session, String a, String b,
                                   long degree, double phi, String source) {
        session.run(
            "MATCH (a:Sample {sampleId: $a}), (b:Sample {sampleId: $b}) "
          + "CREATE (a)-[:RELATIVE {relationship: 'rel', degree: $degree, "
          + "phi: $phi, ibs0_frac: 0.0, source: $source, "
          + "created_at: datetime()}]->(b)",
            Map.of("a", a, "b", b, "degree", degree,
                   "phi", phi, "source", source));
    }

    private static void createIbd(Session session, String a, String b,
                                   long lengthBp, String source) {
        session.run(
            "MATCH (a:Sample {sampleId: $a}), (b:Sample {sampleId: $b}) "
          + "CREATE (a)-[:IBD_SEGMENT {chr: '1', start: 0, end: $end, "
          + "length_bp: $length_bp, source: $source}]->(b)",
            Map.of("a", a, "b", b, "end", lengthBp,
                   "length_bp", lengthBp, "source", source));
    }

    @AfterAll
    static void tearDown() {
        if (driver != null) driver.close();
        if (embeddedNeo4j != null) embeddedNeo4j.close();
    }

    private static Map<String, Long> familySizesById(List<Record> rows) {
        // Per-sample size collapsed to per-family-id size.
        Map<String, Long> bySample = new HashMap<>();
        for (Record r : rows) {
            bySample.put(r.get("family_id").asString(),
                         r.get("family_size").asLong());
        }
        return bySample;
    }

    @Test
    void families_default_max_degree_2_three_trios_plus_singletons() {
        try (Session session = driver.session()) {
            // Default max_degree=2: all three trios (some via 2nd-degree
            // edges) cluster, plus 3 :Sample singletons (S9/S10/S11) and
            // 4 from the IBD fixture (T0..T3 — no :RELATIVE edges) and
            // their pairs are not connected via :RELATIVE.
            List<Record> rows = session.run(
                "CALL graphpop.relate.families('king', {persist: false}) "
              + "YIELD sample_id, family_id, family_size, method "
              + "RETURN sample_id, family_id, family_size, method "
              + "ORDER BY sample_id").list();
            // 12 :Sample S* + 4 :Sample T* = 16 rows.
            assertEquals(16, rows.size());

            Map<String, String> sampleToFamily = new HashMap<>();
            Map<String, Long> sampleToSize = new HashMap<>();
            for (Record r : rows) {
                sampleToFamily.put(r.get("sample_id").asString(),
                                   r.get("family_id").asString());
                sampleToSize.put(r.get("sample_id").asString(),
                                 r.get("family_size").asLong());
                assertEquals("degree", r.get("method").asString());
            }

            // Trio 1: S0/S1/S2 share family_id, size 3.
            assertEquals(sampleToFamily.get("S0"), sampleToFamily.get("S1"));
            assertEquals(sampleToFamily.get("S1"), sampleToFamily.get("S2"));
            assertEquals(3L, sampleToSize.get("S0"));

            // Trio 2 (S3/S4/S5) and Trio 3 (S6/S7/S8) likewise.
            assertEquals(sampleToFamily.get("S3"), sampleToFamily.get("S5"));
            assertEquals(3L, sampleToSize.get("S3"));
            assertEquals(sampleToFamily.get("S6"), sampleToFamily.get("S8"));
            assertEquals(3L, sampleToSize.get("S6"));

            // Singletons: S9/S10/S11 each form their own family.
            assertEquals("S9", sampleToFamily.get("S9"));
            assertEquals(1L, sampleToSize.get("S9"));
            assertEquals("S10", sampleToFamily.get("S10"));
            assertEquals(1L, sampleToSize.get("S10"));
            assertEquals("S11", sampleToFamily.get("S11"));
        }
    }

    @Test
    void families_max_degree_1_drops_second_degree_edges() {
        try (Session session = driver.session()) {
            // With max_degree=1, only first-degree :RELATIVE edges count.
            // Trio 2 (S3-S4 first-deg, S4-S5 second-deg) shrinks: {S3,S4} pair
            // and {S5} singleton. Trio 3 (all second-degree) splits into
            // three singletons {S6}, {S7}, {S8}.
            List<Record> rows = session.run(
                "CALL graphpop.relate.families('king', "
              + "{max_degree: 1, persist: false}) "
              + "YIELD sample_id, family_id, family_size "
              + "RETURN sample_id, family_id, family_size").list();
            Map<String, Long> sizeBySample = new HashMap<>();
            for (Record r : rows) {
                sizeBySample.put(r.get("sample_id").asString(),
                                 r.get("family_size").asLong());
            }
            // Trio 1 unchanged.
            assertEquals(3L, sizeBySample.get("S0"));
            // Trio 2 shrinks.
            assertEquals(2L, sizeBySample.get("S3"));
            assertEquals(2L, sizeBySample.get("S4"));
            assertEquals(1L, sizeBySample.get("S5"));
            // Trio 3 splits.
            assertEquals(1L, sizeBySample.get("S6"));
            assertEquals(1L, sizeBySample.get("S7"));
            assertEquals(1L, sizeBySample.get("S8"));
        }
    }

    @Test
    void families_ibd_bp_predicate_links_above_threshold() {
        try (Session session = driver.session()) {
            List<Record> rows = session.run(
                "CALL graphpop.relate.families('test_ibd', "
              + "{min_total_ibd_bp: 5000000, persist: false}) "
              + "YIELD sample_id, family_id, family_size, method "
              + "WHERE sample_id IN ['T0','T1','T2','T3'] "
              + "RETURN sample_id, family_id, family_size, method").list();
            assertEquals(4, rows.size());
            Map<String, Long> sizeBySample = new HashMap<>();
            Map<String, String> familyBySample = new HashMap<>();
            for (Record r : rows) {
                sizeBySample.put(r.get("sample_id").asString(),
                                 r.get("family_size").asLong());
                familyBySample.put(r.get("sample_id").asString(),
                                   r.get("family_id").asString());
                assertEquals("ibd_bp", r.get("method").asString());
            }
            // T0-T1 above threshold → cluster of size 2.
            assertEquals(2L, sizeBySample.get("T0"));
            assertEquals(familyBySample.get("T0"), familyBySample.get("T1"));
            // T2-T3 below threshold → singletons.
            assertEquals(1L, sizeBySample.get("T2"));
            assertEquals(1L, sizeBySample.get("T3"));
            assertNotEquals(familyBySample.get("T2"), familyBySample.get("T3"));
        }
    }

    @Test
    void families_min_phi_predicate_filters_by_phi() {
        try (Session session = driver.session()) {
            // min_phi=0.20 keeps only first-degree edges (phi=0.25),
            // dropping all second-degree ones (phi=0.125).
            List<Record> rows = session.run(
                "CALL graphpop.relate.families('king', "
              + "{min_phi: 0.20, persist: false}) "
              + "YIELD sample_id, family_id, family_size, method "
              + "WHERE sample_id IN ['S0','S1','S2','S3','S4','S5','S6','S7','S8'] "
              + "RETURN sample_id, family_id, family_size, method").list();
            Map<String, Long> sizeBySample = new HashMap<>();
            for (Record r : rows) {
                sizeBySample.put(r.get("sample_id").asString(),
                                 r.get("family_size").asLong());
                assertEquals("phi", r.get("method").asString());
            }
            // Trio 1 (both phi=0.25): size 3.
            assertEquals(3L, sizeBySample.get("S0"));
            // Trio 2 (S3-S4 phi=0.25, S4-S5 phi=0.125): {S3,S4}, {S5}.
            assertEquals(2L, sizeBySample.get("S3"));
            assertEquals(1L, sizeBySample.get("S5"));
            // Trio 3 (all phi=0.125): three singletons.
            assertEquals(1L, sizeBySample.get("S6"));
            assertEquals(1L, sizeBySample.get("S7"));
            assertEquals(1L, sizeBySample.get("S8"));
        }
    }

    @Test
    void families_persist_writes_in_family_and_family_nodes() {
        try (Session session = driver.session()) {
            // Idempotent re-run: call twice, count edges + Family nodes.
            session.run(
                "CALL graphpop.relate.families('king', {max_degree: 2}) "
              + "YIELD sample_id RETURN count(*) AS c").single();
            session.run(
                "CALL graphpop.relate.families('king', {max_degree: 2}) "
              + "YIELD sample_id RETURN count(*) AS c").single();

            long edgeCount = session.run(
                "MATCH (:Sample)-[r:IN_FAMILY {source:'king'}]->(:Family) "
              + "RETURN count(r) AS c"
            ).single().get("c").asLong();
            // Every Sample lands in a family — including the 4 T* samples
            // which are singletons under source='king'. 12 S* + 4 T* = 16.
            assertEquals(16L, edgeCount);

            // 3 trios (3 families) + 3 S-singletons + 4 T-singletons = 10.
            long familyCount = session.run(
                "MATCH (f:Family {source:'king'}) RETURN count(f) AS c"
            ).single().get("c").asLong();
            assertEquals(10L, familyCount);
        }
    }

    @Test
    void families_empty_source_each_sample_is_singleton() {
        try (Session session = driver.session()) {
            List<Record> rows = session.run(
                "CALL graphpop.relate.families('does_not_exist', "
              + "{persist: false}) "
              + "YIELD sample_id, family_id, family_size "
              + "RETURN sample_id, family_id, family_size").list();
            // Every :Sample appears as its own family.
            assertEquals(16, rows.size());
            for (Record r : rows) {
                assertEquals(r.get("sample_id").asString(),
                             r.get("family_id").asString());
                assertEquals(1L, r.get("family_size").asLong());
            }
        }
    }
}
