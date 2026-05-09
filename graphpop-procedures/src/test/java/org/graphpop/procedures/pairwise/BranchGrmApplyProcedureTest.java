package org.graphpop.procedures.pairwise;

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
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Cypher integration test for {@link BranchGrmApplyProcedure}.
 *
 * <p>Headline correctness gate: the procedure-emitted {@code G · v}
 * must match the Java in-memory {@code BranchGrmComputer · v}
 * element-wise on the same 20-sample msprime fixture.</p>
 */
class BranchGrmApplyProcedureTest {

    private static final String RUN_ID = "egrm_fixture_20samples";

    private static Neo4j embeddedNeo4j;
    private static Driver driver;
    private static int nSamples;
    private static double[][] fullMatrix;

    @BeforeAll
    static void setUp() throws Exception {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(BranchGrmApplyProcedure.class)
                .withProcedure(BranchGrmProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());

        ObjectMapper mapper = new ObjectMapper();
        JsonNode arg;
        try (InputStream is = BranchGrmApplyProcedureTest.class.getResourceAsStream(
                "/egrm_fixture_20samples_arg.json")) {
            arg = mapper.readTree(is);
        }
        JsonNode expected;
        try (InputStream is = BranchGrmApplyProcedureTest.class.getResourceAsStream(
                "/egrm_expected_20samples.json")) {
            expected = mapper.readTree(is);
        }
        nSamples = expected.get("n_samples").asInt();
        fullMatrix = new double[nSamples][nSamples];
        JsonNode mat = expected.get("matrix");
        for (int i = 0; i < nSamples; i++)
            for (int j = 0; j < nSamples; j++)
                fullMatrix[i][j] = mat.get(i).get(j).asDouble();

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
            for (JsonNode e : arg.get("edges")) {
                Map<String, Object> p = new HashMap<>();
                p.put("pid", RUN_ID + ":" + e.get("parent").asInt());
                p.put("cid", RUN_ID + ":" + e.get("child").asInt());
                p.put("rid", RUN_ID);
                p.put("s", e.get("start").asLong());
                p.put("en", e.get("end").asLong());
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

    private static double[] matMulVec(double[] v) {
        double[] out = new double[nSamples];
        for (int i = 0; i < nSamples; i++) {
            double s = 0.0;
            for (int j = 0; j < nSamples; j++) s += fullMatrix[i][j] * v[j];
            out[i] = s;
        }
        return out;
    }

    private static double[] runApply(List<Double> v, Map<String, Object> opts) {
        try (Session session = driver.session()) {
            List<Record> rows = session.run(
                "CALL graphpop.kinship.branch_grm_apply($rid, $v, $opts) "
              + "YIELD sample_id, value RETURN sample_id, value",
                Map.of("rid", RUN_ID, "v", v, "opts", opts)).list();
            double[] out = new double[nSamples];
            for (Record r : rows) {
                int idx = parseHap(r.get("sample_id").asString());
                out[idx] = r.get("value").asDouble();
            }
            return out;
        }
    }

    private static int parseHap(String sid) {
        int colon = sid.indexOf(':');
        String head = (colon >= 0) ? sid.substring(0, colon) : sid;
        return Integer.parseInt(head.substring("hap_".length()));
    }

    @Test
    void apply_matchesFullMatrixVecMul_randomVector() {
        Random rng = new Random(42);
        List<Double> v = new java.util.ArrayList<>();
        double[] vd = new double[nSamples];
        for (int i = 0; i < nSamples; i++) {
            vd[i] = rng.nextGaussian();
            v.add(vd[i]);
        }
        double[] expected = matMulVec(vd);
        double[] actual = runApply(v, Map.of());

        for (int i = 0; i < nSamples; i++) {
            double diff = Math.abs(expected[i] - actual[i]);
            double denom = Math.max(Math.abs(expected[i]), 1e-12);
            assertTrue(diff / denom < 1e-6 && diff < 1e-9,
                String.format("apply[%d]: expected=%.10f actual=%.10f diff=%.3e",
                    i, expected[i], actual[i], diff));
        }
    }

    @Test
    void apply_zeroVectorReturnsZero() {
        List<Double> v = new java.util.ArrayList<>();
        for (int i = 0; i < nSamples; i++) v.add(0.0);
        double[] gv = runApply(v, Map.of());
        for (int i = 0; i < nSamples; i++)
            assertEquals(0.0, gv[i], 1e-12);
    }

    @Test
    void apply_constantVectorReturnsZero() {
        List<Double> v = new java.util.ArrayList<>();
        for (int i = 0; i < nSamples; i++) v.add(2.5);
        double[] gv = runApply(v, Map.of());
        for (int i = 0; i < nSamples; i++)
            assertEquals(0.0, gv[i], 1e-9, "constant v -> Gv = 0 (centring)");
    }

    @Test
    void apply_rejectsWrongLengthVector() {
        try (Session session = driver.session()) {
            assertThrows(Exception.class, () ->
                session.run(
                    "CALL graphpop.kinship.branch_grm_apply($rid, $v) "
                  + "YIELD sample_id RETURN sample_id",
                    Map.of("rid", RUN_ID, "v", List.of(1.0, 2.0))).list()
            );
        }
    }

    @Test
    void apply_composesWithRestrictToPathway() throws Exception {
        // Procedure must run end-to-end with the conditional predicate;
        // we do not validate exact numbers vs a Python reference here
        // (the conditional matrix's egrm_pathway_half.json doesn't have
        // an annotation layer ingested in this test setUp) -- the
        // BranchGrmComputer/Procedure tests already cover the
        // conditional path. We just verify the call succeeds.
        Random rng = new Random(7);
        List<Double> v = new java.util.ArrayList<>();
        for (int i = 0; i < nSamples; i++) v.add(rng.nextGaussian());

        // Without annotation graph, restrict_to_pathway returns the
        // empty lit-set so all branches get weight 0. The procedure
        // then returns a zero vector (no error).
        double[] gv = runApply(v, Map.of("restrict_to_pathway", "P_test"));
        for (int i = 0; i < nSamples; i++)
            assertEquals(0.0, gv[i], 1e-12,
                "no annotation layer -> all branches dark -> G·v = 0");
    }

    @Test
    void apply_block_kVectorsSameAsKApplyCalls() {
        try (Session session = driver.session()) {
            Random rng = new Random(99);
            List<List<Double>> vectors = new java.util.ArrayList<>();
            for (int k = 0; k < 3; k++) {
                List<Double> v = new java.util.ArrayList<>();
                for (int i = 0; i < nSamples; i++) v.add(rng.nextGaussian());
                vectors.add(v);
            }
            // Block call.
            List<Record> blockRows = session.run(
                "CALL graphpop.kinship.branch_grm_apply($rid, [], $opts) "
              + "YIELD sample_id, col, value RETURN sample_id, col, value",
                Map.of("rid", RUN_ID, "opts", Map.of("vectors", vectors))).list();
            assertEquals(3 * nSamples, blockRows.size());

            // Per-vector calls and compare element-wise.
            for (int k = 0; k < 3; k++) {
                double[] perCall = runApply(vectors.get(k), Map.of());
                double[] fromBlock = new double[nSamples];
                for (Record r : blockRows) {
                    if (r.get("col").asLong() == k) {
                        fromBlock[parseHap(r.get("sample_id").asString())] =
                                r.get("value").asDouble();
                    }
                }
                for (int i = 0; i < nSamples; i++)
                    assertEquals(perCall[i], fromBlock[i], 1e-15,
                        "block[" + k + "][" + i + "]");
            }
        }
    }
}
