package org.graphpop.procedures.pairwise;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.CholeskyDecomposition;
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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Cypher integration test for {@link BranchGrmHeProcedure}.
 *
 * <p>Headline correctness gate: simulated heritability {@code h² =
 * 0.5} from {@code y ~ N(0, h²·B + (1−h²)·I)} where {@code B} is
 * the 20-sample fixture's branch GRM, recovered within ±0.25 (the
 * tolerance is large because n=20).</p>
 */
class BranchGrmHeProcedureTest {

    private static final String RUN_ID = "egrm_fixture_20samples";

    private static Neo4j embeddedNeo4j;
    private static Driver driver;
    private static int nSamples;
    private static double[][] fullMatrix;

    @BeforeAll
    static void setUp() throws Exception {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(BranchGrmHeProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());

        ObjectMapper mapper = new ObjectMapper();
        JsonNode arg;
        try (InputStream is = BranchGrmHeProcedureTest.class.getResourceAsStream(
                "/egrm_fixture_20samples_arg.json")) {
            arg = mapper.readTree(is);
        }
        JsonNode expected;
        try (InputStream is = BranchGrmHeProcedureTest.class.getResourceAsStream(
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

    /** Simulate y ~ N(0, h²·B + (1−h²)·I). */
    private static List<Double> simulateY(double h2, long seed) {
        // covariance = h² · B + (1−h²) · I
        double[][] cov = new double[nSamples][nSamples];
        for (int i = 0; i < nSamples; i++) {
            for (int j = 0; j < nSamples; j++) {
                cov[i][j] = h2 * fullMatrix[i][j];
                if (i == j) cov[i][j] += (1.0 - h2);
            }
        }
        // The branch GRM may have small negative eigenvalues numerically;
        // add a tiny ridge to ensure PSD.
        for (int i = 0; i < nSamples; i++) cov[i][i] += 1e-9;
        RealMatrix C = new Array2DRowRealMatrix(cov, false);
        CholeskyDecomposition chol = new CholeskyDecomposition(
                C, 1e-10, 1e-10);
        Random rng = new Random(seed);
        double[] z = new double[nSamples];
        for (int i = 0; i < nSamples; i++) z[i] = rng.nextGaussian();
        double[] y = chol.getL().operate(z);
        List<Double> out = new ArrayList<>(nSamples);
        for (double v : y) out.add(v);
        return out;
    }

    @Test
    void recovers_h2_from_simulation_within_tolerance() {
        List<Double> y = simulateY(0.5, 42L);
        try (Session session = driver.session()) {
            Record r = session.run(
                "CALL graphpop.kinship.branch_grm_he($rid, $y, "
              + "{n_hutchinson: 100, seed: 7}) "
              + "YIELD h2, se, num, tr_g_sq, n_samples, n_hutchinson "
              + "RETURN h2, se, num, tr_g_sq, n_samples, n_hutchinson",
                Map.of("rid", RUN_ID, "y", y)).single();
            double h2 = r.get("h2").asDouble();
            double se = r.get("se").asDouble();
            assertEquals(0.5, h2, 0.25,
                "h² recovery within ±0.25 (large for n=20): got " + h2);
            assertTrue(se > 0 && Double.isFinite(se), "SE must be positive finite");
            assertEquals(nSamples, r.get("n_samples").asLong());
            assertEquals(100L, r.get("n_hutchinson").asLong());
        }
    }

    @Test
    void constant_phenotype_returns_h2_zero() {
        List<Double> y = new ArrayList<>();
        for (int i = 0; i < nSamples; i++) y.add(1.0);
        try (Session session = driver.session()) {
            Record r = session.run(
                "CALL graphpop.kinship.branch_grm_he($rid, $y) "
              + "YIELD h2, num RETURN h2, num",
                Map.of("rid", RUN_ID, "y", y)).single();
            assertEquals(0.0, r.get("h2").asDouble(), 1e-12);
            assertEquals(0.0, r.get("num").asDouble(), 1e-12);
        }
    }

    @Test
    void wrong_length_phenotype_throws() {
        try (Session session = driver.session()) {
            assertThrows(Exception.class, () ->
                session.run(
                    "CALL graphpop.kinship.branch_grm_he($rid, [1.0, 2.0]) "
                  + "YIELD h2 RETURN h2",
                    Map.of("rid", RUN_ID)).single()
            );
        }
    }

    @Test
    void hutchinson_convergence_variance_decreases_with_M() {
        // For larger M, repeated runs should give more concentrated h2
        // estimates. This is not a strict test but a sanity check.
        List<Double> y = simulateY(0.5, 13L);
        try (Session session = driver.session()) {
            double[] h2_M10 = new double[5];
            double[] h2_M200 = new double[5];
            for (int i = 0; i < 5; i++) {
                h2_M10[i] = session.run(
                    "CALL graphpop.kinship.branch_grm_he($rid, $y, "
                  + "{n_hutchinson: 10, seed: $seed}) YIELD h2 RETURN h2",
                    Map.of("rid", RUN_ID, "y", y, "seed", (long) i)
                ).single().get("h2").asDouble();
                h2_M200[i] = session.run(
                    "CALL graphpop.kinship.branch_grm_he($rid, $y, "
                  + "{n_hutchinson: 200, seed: $seed}) YIELD h2 RETURN h2",
                    Map.of("rid", RUN_ID, "y", y, "seed", (long) i)
                ).single().get("h2").asDouble();
            }
            assertTrue(stddev(h2_M200) <= stddev(h2_M10) + 1e-3,
                "Larger M should give less variance: "
              + "std(M=200)=" + stddev(h2_M200)
              + " std(M=10)=" + stddev(h2_M10));
        }
    }

    private static double stddev(double[] xs) {
        double mean = 0;
        for (double x : xs) mean += x;
        mean /= xs.length;
        double sq = 0;
        for (double x : xs) sq += (x - mean) * (x - mean);
        return Math.sqrt(sq / xs.length);
    }

    @Test
    void unknown_run_returns_empty() {
        try (Session session = driver.session()) {
            List<Double> y = new ArrayList<>();
            for (int i = 0; i < nSamples; i++) y.add(0.1);
            long c = session.run(
                "CALL graphpop.kinship.branch_grm_he('does_not_exist', $y) "
              + "YIELD h2 RETURN count(*) AS c",
                Map.of("y", y)).single().get("c").asLong();
            assertEquals(0L, c);
        }
    }
}
