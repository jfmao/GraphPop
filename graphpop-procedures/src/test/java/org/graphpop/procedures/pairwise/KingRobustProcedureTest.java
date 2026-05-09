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
 * Cypher-level integration test for {@code graphpop.kinship.king}.
 *
 * <p>Uses the same 4-sample x 8-variant fixture as
 * {@link KingRobustComputerTest}; expected kinship values are hand-computed.</p>
 */
class KingRobustProcedureTest {

    private static final double EPS = 1e-12;

    private static Neo4j embeddedNeo4j;
    private static Driver driver;

    @BeforeAll
    static void setUp() {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(KingRobustProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());

        // 4-sample x 8-variant fixture (matches KingRobustComputerTest).
        int[][] gt = {
            // V1, V2, V3, V4, V5, V6, V7, V8
            { 1, 0, 2, 1, 0, 1, 2, 1 },  // S1 (packed_index 0)
            { 1, 1, 1, 1, 1, 0, 2, 0 },  // S2 (packed_index 1)
            { 2, 2, 0, 0, 2, 0, 1, 2 },  // S3 (packed_index 2)
            { 0, 0, 0, 0, 0, 0, 0, 0 },  // S4 (packed_index 3)
        };
        int nSamples = gt.length;
        int nVariants = gt[0].length;

        try (Session session = driver.session()) {
            session.run("CREATE (:Population {populationId: 'POP1', n_samples: 4})");

            for (int s = 0; s < nSamples; s++) {
                session.run(
                    "MATCH (p:Population {populationId: 'POP1'}) "
                  + "CREATE (s:Sample {sampleId: $sid, packed_index: $pi})-[:IN_POPULATION]->(p)",
                    Map.of("sid", "S" + (s + 1), "pi", s));
            }

            for (int v = 0; v < nVariants; v++) {
                byte[] gtPacked = new byte[PackedGenotypeReader.gtPackedLength(nSamples)];
                for (int s = 0; s < nSamples; s++) {
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
    void king_allInformativePairs() {
        try (Session session = driver.session()) {
            // Tests pass small fixture so set min_snp=0 and min_phi well below smallest expected.
            List<Record> rows = session.run(
                "CALL graphpop.kinship.king('chr1', 'POP1', "
              + "{min_snp: 0, min_phi: -10.0}) "
              + "YIELD sample_a, sample_b, phi, ibs0, het_het, n_snp, n_aa_min, method "
              + "RETURN sample_a, sample_b, phi, ibs0, het_het, n_snp, n_aa_min"
            ).list();

            // Three pairs are informative (S4 produces NaN for every pair => filtered out):
            //   (S1,S2)=0.25, (S1,S3)=-3.0, (S2,S3)=-1.0
            assertEquals(3, rows.size(),
                "Expected 3 informative non-self pairs; got " + rows.size());

            Map<String, Map<String, Value>> byPair = new HashMap<>();
            for (Record r : rows) {
                String key = r.get("sample_a").asString() + "," + r.get("sample_b").asString();
                byPair.put(key, r.asMap(v -> v));
            }

            assertEquals(0.25,
                byPair.get("S1,S2").get("phi").asDouble(), EPS);
            assertEquals(2L,
                byPair.get("S1,S2").get("het_het").asLong());
            assertEquals(0L,
                byPair.get("S1,S2").get("ibs0").asLong());
            assertEquals(8L,
                byPair.get("S1,S2").get("n_snp").asLong());
            assertEquals(4L,
                byPair.get("S1,S2").get("n_aa_min").asLong());

            assertEquals(-3.0,
                byPair.get("S1,S3").get("phi").asDouble(), EPS);
            assertEquals(3L,
                byPair.get("S1,S3").get("ibs0").asLong());

            assertEquals(-1.0,
                byPair.get("S2,S3").get("phi").asDouble(), EPS);
        }
    }

    @Test
    void king_minPhiFilterDropsLowPairs() {
        try (Session session = driver.session()) {
            // Default min_phi = -0.5 (KING convention). Only (S1,S2) passes.
            List<Record> rows = session.run(
                "CALL graphpop.kinship.king('chr1', 'POP1', {min_snp: 0}) "
              + "YIELD sample_a, sample_b, phi RETURN sample_a, sample_b, phi"
            ).list();
            assertEquals(1, rows.size(), "Only the parent-child pair clears default min_phi");
            assertEquals("S1", rows.get(0).get("sample_a").asString());
            assertEquals("S2", rows.get(0).get("sample_b").asString());
            assertEquals(0.25, rows.get(0).get("phi").asDouble(), EPS);
        }
    }

    @Test
    void king_minSnpFilterDropsAllPairs() {
        try (Session session = driver.session()) {
            // Only 8 variants in the fixture; default min_snp=1000 => no rows.
            List<Record> rows = session.run(
                "CALL graphpop.kinship.king('chr1', 'POP1', {}) "
              + "YIELD sample_a RETURN sample_a"
            ).list();
            assertEquals(0, rows.size(),
                "Default min_snp=1000 must filter out the 8-variant fixture");
        }
    }

    @Test
    void king_includeSelfFlagEmitsDiagonal() {
        try (Session session = driver.session()) {
            // include_self=true with min_phi=-10 emits self-pairs for the 3 samples
            // that have at least one het site. S4 has zero het sites and produces NaN.
            List<Record> rows = session.run(
                "CALL graphpop.kinship.king('chr1', 'POP1', "
              + "{min_snp: 0, min_phi: -10.0, include_self: true}) "
              + "YIELD sample_a, sample_b, phi RETURN sample_a, sample_b, phi"
            ).list();

            long selfCount = rows.stream()
                .filter(r -> r.get("sample_a").asString().equals(r.get("sample_b").asString()))
                .count();
            assertEquals(3, selfCount, "S1, S2, S3 self-pairs emit; S4 is undefined");

            for (Record r : rows) {
                if (r.get("sample_a").asString().equals(r.get("sample_b").asString())) {
                    assertEquals(0.5, r.get("phi").asDouble(), EPS,
                        "Self-pair phi by KING formula is 0.5 when N_Aa > 0");
                }
            }
        }
    }

    @Test
    void king_explicitSamplesSubset() {
        try (Session session = driver.session()) {
            // Pass the samples option to override pop-based lookup.
            List<Record> rows = session.run(
                "CALL graphpop.kinship.king('chr1', 'POP1', "
              + "{min_snp: 0, min_phi: -10.0, samples: ['S1','S2']}) "
              + "YIELD sample_a, sample_b, phi RETURN sample_a, sample_b, phi"
            ).list();
            assertEquals(1, rows.size(), "Subset of 2 samples => 1 pair");
            assertEquals(0.25, rows.get(0).get("phi").asDouble(), EPS);
        }
    }
}
