package org.graphpop.procedures.pairwise;

import org.graphpop.procedures.PackedGenotypeReader;
import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.neo4j.driver.Driver;
import org.neo4j.driver.GraphDatabase;
import org.neo4j.driver.Record;
import org.neo4j.driver.Session;
import org.neo4j.driver.Value;
import org.neo4j.harness.Neo4j;
import org.neo4j.harness.Neo4jBuilders;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Cypher integration test for {@code graphpop.kinship.ibs}. Reuses the
 * 4-sample × 8-variant fixture from {@link KingRobustProcedureTest}.
 */
class IbsProcedureTest {

    private static final double EPS = 1e-12;

    private static Neo4j embeddedNeo4j;
    private static Driver driver;

    @BeforeAll
    static void setUp() {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(IbsProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());

        int[][] gt = {
            { 1, 0, 2, 1, 0, 1, 2, 1 },
            { 1, 1, 1, 1, 1, 0, 2, 0 },
            { 2, 2, 0, 0, 2, 0, 1, 2 },
            { 0, 0, 0, 0, 0, 0, 0, 0 },
        };
        int nS = gt.length;
        int nV = gt[0].length;

        try (Session session = driver.session()) {
            session.run("CREATE (:Population {populationId: 'POP1', n_samples: 4})");
            for (int s = 0; s < nS; s++) {
                session.run(
                    "MATCH (p:Population {populationId: 'POP1'}) "
                  + "CREATE (:Sample {sampleId: $sid, packed_index: $pi})-[:IN_POPULATION]->(p)",
                    Map.of("sid", "S" + (s + 1), "pi", s));
            }
            for (int v = 0; v < nV; v++) {
                byte[] gtPacked = new byte[PackedGenotypeReader.gtPackedLength(nS)];
                for (int s = 0; s < nS; s++) {
                    PackedGenotypeReader.setGenotype(gtPacked, s, gt[s][v]);
                }
                session.run(
                    "CREATE (:Variant {variantId: $vid, chr: 'chr1', pos: $pos, "
                  + "ref: 'A', alt: 'T', gt_packed: $gt})",
                    Map.of("vid", "chr1:" + (v * 100 + 100),
                           "pos", (long) (v * 100 + 100),
                           "gt", gtPacked));
            }
        }
    }

    @AfterAll
    static void tearDown() {
        if (driver != null) driver.close();
        if (embeddedNeo4j != null) embeddedNeo4j.close();
    }

    @Test
    void ibs_allInformativePairs() {
        try (Session session = driver.session()) {
            List<Record> rows = session.run(
                "CALL graphpop.kinship.ibs('chr1', 'POP1', {min_snp: 0}) "
              + "YIELD sample_a, sample_b, phi, ibs0, het_het, n_snp, method "
              + "RETURN sample_a, sample_b, phi, ibs0, het_het, n_snp"
            ).list();

            // 6 informative non-self pairs.
            assertEquals(6, rows.size());

            Map<String, Map<String, Value>> byPair = new HashMap<>();
            for (Record r : rows) {
                String key = r.get("sample_a").asString() + "," + r.get("sample_b").asString();
                byPair.put(key, r.asMap(v -> v));
            }
            // Hand-computed values from IbsComputerTest.
            assertEquals(11.0 / 16.0, byPair.get("S1,S2").get("phi").asDouble(), EPS);
            assertEquals(5.0 / 16.0, byPair.get("S1,S3").get("phi").asDouble(), EPS);
            assertEquals(8.0 / 16.0, byPair.get("S1,S4").get("phi").asDouble(), EPS);
            assertEquals(8.0 / 16.0, byPair.get("S2,S3").get("phi").asDouble(), EPS);
            assertEquals(9.0 / 16.0, byPair.get("S2,S4").get("phi").asDouble(), EPS);
            assertEquals(7.0 / 16.0, byPair.get("S3,S4").get("phi").asDouble(), EPS);

            // Counters spot-check on (S1, S2): ibs0=0, het_het=ibs2=3, n_snp=8.
            assertEquals(0L, byPair.get("S1,S2").get("ibs0").asLong());
            assertEquals(3L, byPair.get("S1,S2").get("het_het").asLong());
            assertEquals(8L, byPair.get("S1,S2").get("n_snp").asLong());
        }
    }

    @Test
    void ibs_minSnpFilterDropsAllPairs() {
        try (Session session = driver.session()) {
            // Default min_snp=1000; fixture has 8 variants -> all filtered.
            long c = session.run(
                "CALL graphpop.kinship.ibs('chr1', 'POP1', {}) "
              + "YIELD sample_a RETURN count(*) AS c"
            ).single().get("c").asLong();
            assertEquals(0L, c);
        }
    }

    @Test
    void ibs_minIbsFilter() {
        try (Session session = driver.session()) {
            // Only pairs with IBS >= 0.6 -> just (S1, S2) at 0.6875.
            List<Record> rows = session.run(
                "CALL graphpop.kinship.ibs('chr1', 'POP1', {min_snp: 0, min_ibs: 0.6}) "
              + "YIELD sample_a, sample_b, phi RETURN sample_a, sample_b, phi"
            ).list();
            assertEquals(1, rows.size());
            assertEquals("S1", rows.get(0).get("sample_a").asString());
            assertEquals("S2", rows.get(0).get("sample_b").asString());
        }
    }

    @Test
    void ibs_includeSelfReturnsDiagonalOne() {
        try (Session session = driver.session()) {
            List<Record> rows = session.run(
                "CALL graphpop.kinship.ibs('chr1', 'POP1', "
              + "{min_snp: 0, include_self: true}) "
              + "YIELD sample_a, sample_b, phi WHERE sample_a = sample_b "
              + "RETURN sample_a, phi"
            ).list();
            assertEquals(4, rows.size());
            for (Record r : rows) assertEquals(1.0, r.get("phi").asDouble(), EPS);
        }
    }

    @Test
    void ibs_methodTagIsIbs() {
        try (Session session = driver.session()) {
            String method = session.run(
                "CALL graphpop.kinship.ibs('chr1', 'POP1', {min_snp: 0}) "
              + "YIELD method RETURN method LIMIT 1"
            ).single().get("method").asString();
            assertEquals("ibs", method);
        }
    }
}
