package org.graphpop.procedures.pairwise;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
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
import java.util.List;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Cypher integration test for {@link BranchGrmPcaProcedure}.
 *
 * <p>Headline correctness gate: the procedure-emitted top-K
 * eigenvalues must match Apache Commons Math's full-matrix
 * eigendecomposition of the {@code BranchGrmComputer} output to
 * relative error &lt; 10⁻⁴ on the 20-sample fixture.</p>
 */
class BranchGrmPcaProcedureTest {

    private static final String RUN_ID = "egrm_fixture_20samples";
    private static final int K = 3;

    private static Neo4j embeddedNeo4j;
    private static Driver driver;
    private static int nSamples;
    private static double[] expectedTopEvals;

    @BeforeAll
    static void setUp() throws Exception {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(BranchGrmPcaProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());

        ObjectMapper mapper = new ObjectMapper();
        JsonNode arg;
        try (InputStream is = BranchGrmPcaProcedureTest.class.getResourceAsStream(
                "/egrm_fixture_20samples_arg.json")) {
            arg = mapper.readTree(is);
        }
        JsonNode expected;
        try (InputStream is = BranchGrmPcaProcedureTest.class.getResourceAsStream(
                "/egrm_expected_20samples.json")) {
            expected = mapper.readTree(is);
        }
        nSamples = expected.get("n_samples").asInt();

        // Build full matrix and compute reference eigenvalues.
        double[][] full = new double[nSamples][nSamples];
        JsonNode mat = expected.get("matrix");
        for (int i = 0; i < nSamples; i++)
            for (int j = 0; j < nSamples; j++)
                full[i][j] = mat.get(i).get(j).asDouble();
        RealMatrix mr = new Array2DRowRealMatrix(full, false);
        double[] e = new EigenDecomposition(mr).getRealEigenvalues();
        java.util.Arrays.sort(e);
        expectedTopEvals = new double[K];
        for (int i = 0; i < K; i++) expectedTopEvals[i] = e[e.length - 1 - i];

        try (Session session = driver.session()) {
            session.run(
                "CREATE (:ARGRun {runId: $rid, source: 'msprime', "
              + "n_samples: $ns, sequence_length: $sl})",
                Map.of("rid", RUN_ID,
                       "ns", arg.get("n_samples").asInt(),
                       "sl", arg.get("sequence_length").asLong()));
            for (JsonNode s : arg.get("samples")) {
                int sn = s.asInt();
                session.run(
                    "CREATE (:Sample {sampleId: $sid, packed_index: $pi})",
                    Map.of("sid", "hap_" + sn, "pi", sn));
            }
            for (JsonNode n : arg.get("nodes")) {
                Map<String, Object> p = new HashMap<>();
                p.put("tid", RUN_ID + ":" + n.get("id").asInt());
                p.put("rid", RUN_ID);
                p.put("nid", n.get("id").asInt());
                p.put("t", n.get("time").asDouble());
                p.put("is", n.get("is_sample").asBoolean());
                p.put("f", n.get("flags").asLong());
                session.run(
                    "CREATE (:TreeNode {treeNodeId: $tid, runId: $rid, "
                  + "nodeId: $nid, time: $t, is_sample: $is, flags: $f})", p);
            }
            for (JsonNode e2 : arg.get("edges")) {
                Map<String, Object> p = new HashMap<>();
                p.put("pid", RUN_ID + ":" + e2.get("parent").asInt());
                p.put("cid", RUN_ID + ":" + e2.get("child").asInt());
                p.put("rid", RUN_ID);
                p.put("s", e2.get("start").asLong());
                p.put("en", e2.get("end").asLong());
                session.run(
                    "MATCH (p:TreeNode {treeNodeId: $pid}), "
                  + "(c:TreeNode {treeNodeId: $cid}) "
                  + "CREATE (p)-[:PARENT_OF {runId: $rid, start: $s, end: $en}]->(c)",
                    p);
            }
            for (JsonNode s : arg.get("samples")) {
                int sn = s.asInt();
                session.run(
                    "MATCH (n:TreeNode {treeNodeId: $tid}), (s:Sample {sampleId: $sid}) "
                  + "CREATE (n)-[:REPRESENTS {runId: $rid, haplotype: 0}]->(s)",
                    Map.of("tid", RUN_ID + ":" + sn,
                           "sid", "hap_" + sn,
                           "rid", RUN_ID));
            }
        }
    }

    @AfterAll
    static void tearDown() {
        if (driver != null) driver.close();
        if (embeddedNeo4j != null) embeddedNeo4j.close();
    }

    private static int parseHap(String sid) {
        int colon = sid.indexOf(':');
        String head = (colon >= 0) ? sid.substring(0, colon) : sid;
        return Integer.parseInt(head.substring("hap_".length()));
    }

    @Test
    void topEigenvalues_matchCommonsMath() {
        try (Session session = driver.session()) {
            // n_iter = 3*K ≈ 9 is enough for biobank-style geometric
            // convergence; for the small fixture push to 15 for safety.
            List<Record> rows = session.run(
                "CALL graphpop.kinship.branch_grm_pca($rid, $k, {n_iter: 15}) "
              + "YIELD pc, eigenvalue RETURN DISTINCT pc, eigenvalue ORDER BY pc",
                Map.of("rid", RUN_ID, "k", (long) K)).list();
            assertEquals(K, rows.size());
            for (int i = 0; i < K; i++) {
                double expected = expectedTopEvals[i];
                double actual = rows.get(i).get("eigenvalue").asDouble();
                double rel = Math.abs(expected - actual)
                        / Math.max(Math.abs(expected), 1e-12);
                assertTrue(rel < 1e-4,
                    String.format("eigenvalue %d: expected=%.6e actual=%.6e rel=%.3e",
                        i, expected, actual, rel));
            }
        }
    }

    @Test
    void emitsKxNRowsTotal() {
        try (Session session = driver.session()) {
            long c = session.run(
                "CALL graphpop.kinship.branch_grm_pca($rid, $k, {}) "
              + "YIELD sample_id RETURN count(*) AS c",
                Map.of("rid", RUN_ID, "k", (long) K)).single().get("c").asLong();
            assertEquals((long) K * nSamples, c);
        }
    }

    @Test
    void everyPCSumsToZero() {
        try (Session session = driver.session()) {
            List<Record> rows = session.run(
                "CALL graphpop.kinship.branch_grm_pca($rid, $k, {}) "
              + "YIELD pc, value WITH pc, sum(value) AS s "
              + "RETURN pc, s ORDER BY pc",
                Map.of("rid", RUN_ID, "k", (long) K)).list();
            for (Record r : rows) {
                assertEquals(0.0, r.get("s").asDouble(), 1e-9,
                    "PC " + r.get("pc").asLong() + " sum non-zero");
            }
        }
    }

    @Test
    void deterministicSeed() {
        try (Session session = driver.session()) {
            List<Record> r1 = session.run(
                "CALL graphpop.kinship.branch_grm_pca($rid, 2, {seed: 42}) "
              + "YIELD sample_id, pc, value RETURN sample_id, pc, value "
              + "ORDER BY pc, sample_id",
                Map.of("rid", RUN_ID)).list();
            List<Record> r2 = session.run(
                "CALL graphpop.kinship.branch_grm_pca($rid, 2, {seed: 42}) "
              + "YIELD sample_id, pc, value RETURN sample_id, pc, value "
              + "ORDER BY pc, sample_id",
                Map.of("rid", RUN_ID)).list();
            assertEquals(r1.size(), r2.size());
            for (int i = 0; i < r1.size(); i++) {
                assertEquals(r1.get(i).get("value").asDouble(),
                              r2.get(i).get("value").asDouble(), 1e-15);
            }
        }
    }

    @Test
    void rejectsKAtOrAboveN() {
        try (Session session = driver.session()) {
            assertThrows(Exception.class, () ->
                session.run(
                    "CALL graphpop.kinship.branch_grm_pca($rid, $k, {}) "
                  + "YIELD sample_id RETURN sample_id LIMIT 1",
                    Map.of("rid", RUN_ID, "k", (long) nSamples)).list()
            );
        }
    }
}
