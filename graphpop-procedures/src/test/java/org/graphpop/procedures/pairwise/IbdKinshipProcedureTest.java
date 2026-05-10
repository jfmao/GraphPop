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
 * Cypher integration test for {@link IbdKinshipProcedure}. Uses
 * hand-authored {@code :IBD_SEGMENT} edges so the expected aggregate
 * is fully controlled.
 */
class IbdKinshipProcedureTest {

    private static Neo4j embeddedNeo4j;
    private static Driver driver;

    @BeforeAll
    static void setUp() {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(IbdKinshipProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());

        try (Session session = driver.session()) {
            // 4 samples; segments hand-authored for deterministic totals.
            for (int i = 0; i < 4; i++) {
                session.run(
                    "CREATE (:Sample {sampleId: $sid})",
                    Map.of("sid", "S" + i));
            }

            // Segments under source "test_bp" with only length_bp.
            // Total genome span: chr1 [0, 100_000_000] -> 100 Mb.
            //   pair (S0, S1): one 30-Mb segment
            //   pair (S0, S2): two segments, total 50 Mb
            //   pair (S2, S3): one 10-Mb segment
            createSeg(session, "S0", "S1", "chr1", 0, 30_000_000L, "test_bp", null);
            createSeg(session, "S0", "S1", "chr1", 99_000_000L, 100_000_000L, "test_bp", null);  // boundary so chr span = 100 Mb
            createSeg(session, "S0", "S2", "chr1", 0, 30_000_000L, "test_bp", null);
            createSeg(session, "S0", "S2", "chr1", 50_000_000L, 70_000_000L, "test_bp", null);
            createSeg(session, "S2", "S3", "chr1", 80_000_000L, 90_000_000L, "test_bp", null);

            // Segments under source "test_cm" with explicit length_cM.
            // Pairs replicate the same lengths.
            createSeg(session, "S0", "S1", "chr1", 0, 30_000_000L, "test_cm", 30.0);
            createSeg(session, "S0", "S1", "chr1", 99_000_000L, 100_000_000L, "test_cm", 1.0);
            createSeg(session, "S0", "S2", "chr1", 0, 30_000_000L, "test_cm", 30.0);
            createSeg(session, "S0", "S2", "chr1", 50_000_000L, 70_000_000L, "test_cm", 20.0);
            createSeg(session, "S2", "S3", "chr1", 80_000_000L, 90_000_000L, "test_cm", 10.0);
        }
    }

    private static void createSeg(Session session, String a, String b,
                                   String chr, long start, long end,
                                   String source, Double cM) {
        Map<String, Object> p = new HashMap<>();
        p.put("a", a); p.put("b", b);
        p.put("chr", chr); p.put("start", start); p.put("end", end);
        p.put("len_bp", end - start);
        p.put("source", source);
        p.put("cm", cM);
        session.run(
            "MATCH (a:Sample {sampleId: $a}), (b:Sample {sampleId: $b}) "
          + "CREATE (a)-[:IBD_SEGMENT {chr: $chr, start: $start, end: $end, "
          + "length_bp: $len_bp, length_cM: $cm, source: $source}]->(b)",
            p);
    }

    @AfterAll
    static void tearDown() {
        if (driver != null) driver.close();
        if (embeddedNeo4j != null) embeddedNeo4j.close();
    }

    @Test
    void bpFallback_methodIsIbdBp() {
        try (Session session = driver.session()) {
            // total genome span = 100 Mb, denom = 200 Mb.
            // (S0, S1) total bp = 30M + 1M = 31M -> phi = 31/200 = 0.155
            // (S0, S2) total bp = 30M + 20M = 50M -> phi = 50/200 = 0.25
            // (S2, S3) total bp = 10M -> phi = 10/200 = 0.05
            List<Record> rows = session.run(
                "CALL graphpop.ibd.kinship('test_bp') "
              + "YIELD sample_a, sample_b, phi, n_snp, ibs0, method "
              + "RETURN sample_a, sample_b, phi, n_snp, ibs0, method "
              + "ORDER BY sample_a, sample_b").list();
            assertEquals(3, rows.size());
            for (Record r : rows) {
                assertEquals("ibd_bp", r.get("method").asString());
            }
            assertEquals(0.155, rows.get(0).get("phi").asDouble(), 1e-9);
            assertEquals(31_000_000L, rows.get(0).get("ibs0").asLong());
            assertEquals(2L, rows.get(0).get("n_snp").asLong());
            assertEquals(0.25, rows.get(1).get("phi").asDouble(), 1e-9);
            assertEquals(0.05, rows.get(2).get("phi").asDouble(), 1e-9);
        }
    }

    @Test
    void cMPath_methodIsIbd() {
        try (Session session = driver.session()) {
            // Default: total_genome_cM estimated from bp/1e6 = 100 cM,
            // denom = 200 cM.
            // (S0, S1) total cM = 31 -> phi = 31/200 = 0.155
            List<Record> rows = session.run(
                "CALL graphpop.ibd.kinship('test_cm') "
              + "YIELD sample_a, sample_b, phi, method "
              + "RETURN sample_a, sample_b, phi, method "
              + "ORDER BY sample_a, sample_b").list();
            for (Record r : rows) assertEquals("ibd", r.get("method").asString());
            assertEquals(0.155, rows.get(0).get("phi").asDouble(), 1e-9);
            assertEquals(0.25, rows.get(1).get("phi").asDouble(), 1e-9);
            assertEquals(0.05, rows.get(2).get("phi").asDouble(), 1e-9);
        }
    }

    @Test
    void cMPath_explicitTotalGenomeCM() {
        try (Session session = driver.session()) {
            // total_genome_cM = 50 cM (override) -> denom = 100 cM.
            // (S0, S1) phi = 31/100 = 0.31
            List<Record> rows = session.run(
                "CALL graphpop.ibd.kinship('test_cm', {total_genome_cM: 50.0}) "
              + "YIELD sample_a, sample_b, phi "
              + "RETURN sample_a, sample_b, phi "
              + "ORDER BY sample_a, sample_b").list();
            assertEquals(0.31, rows.get(0).get("phi").asDouble(), 1e-9);
        }
    }

    @Test
    void minLengthBpFilter() {
        try (Session session = driver.session()) {
            // min_length_bp = 32M -> drops (S0, S1) (31M total) and
            // (S2, S3) (10M total); keeps (S0, S2) (50M total).
            List<Record> rows = session.run(
                "CALL graphpop.ibd.kinship('test_bp', {min_length_bp: 32000000}) "
              + "YIELD sample_a RETURN sample_a").list();
            assertEquals(1, rows.size());
        }
    }

    @Test
    void unknownSource_returnsEmpty() {
        try (Session session = driver.session()) {
            long c = session.run(
                "CALL graphpop.ibd.kinship('does_not_exist') "
              + "YIELD sample_a RETURN count(*) AS c"
            ).single().get("c").asLong();
            assertEquals(0L, c);
        }
    }
}
