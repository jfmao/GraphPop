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

import static org.junit.jupiter.api.Assertions.*;

/**
 * Integration test for {@link PosteriorBranchGrmProcedure}.
 *
 * <p>Headline correctness gates: posterior-mean and posterior-SE
 * matrices over five independent msprime-simulated ARGs match the
 * Python Welford reference (built by
 * {@code build_egrm_fixture_posterior.py}) to relative error
 * &lt; 10⁻⁶.</p>
 */
class PosteriorBranchGrmProcedureTest {

    private static final int[] SEEDS = {42, 43, 44, 45, 46};
    private static final int MID = 20;

    private static Neo4j embeddedNeo4j;
    private static Driver driver;
    private static int nSamples;
    private static List<String> RUN_IDS = new java.util.ArrayList<>();

    @BeforeAll
    static void setUp() throws Exception {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(PosteriorBranchGrmProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());

        // Shared sample / annotation nodes.
        try (Session session = driver.session()) {
            // Look up sample count from the first fixture so we don't hard-code.
            JsonNode firstArg = loadJson(
                "/egrm_posterior_fixture_arg_seed" + SEEDS[0] + ".json");
            nSamples = firstArg.get("n_samples").asInt();

            for (JsonNode s : firstArg.get("samples")) {
                int sampleNodeId = s.asInt();
                session.run(
                    "CREATE (:Sample {sampleId: $sid, packed_index: $pi})",
                    Map.of("sid", "hap_" + sampleNodeId, "pi", sampleNodeId));
            }

            session.run("CREATE (:Pathway {pathwayId: 'P_test', name: 'test pathway'})");
            session.run(
                "MATCH (p:Pathway {pathwayId: 'P_test'}) "
              + "CREATE (g:Gene {geneId: 'G_test', symbol: 'GTEST'})-[:IN_PATHWAY]->(p)");
            session.run(
                "CREATE (g:Gene {geneId: 'G_other', symbol: 'GOTHER'})");
        }

        // Per-run topology + per-run :Variant + :MUTATED_ON.
        for (int seed : SEEDS) {
            String runId = "egrm_post_seed" + seed;
            RUN_IDS.add(runId);
            ingestRun(runId, loadJson(
                "/egrm_posterior_fixture_arg_seed" + seed + ".json"));
        }
    }

    @AfterAll
    static void tearDown() {
        if (driver != null) driver.close();
        if (embeddedNeo4j != null) embeddedNeo4j.close();
    }

    // ---- fixtures -------------------------------------------------------

    private static JsonNode loadJson(String resource) throws Exception {
        ObjectMapper mapper = new ObjectMapper();
        try (InputStream is = PosteriorBranchGrmProcedureTest.class
                .getResourceAsStream(resource)) {
            assertNotNull(is, "missing resource: " + resource);
            return mapper.readTree(is);
        }
    }

    private static void ingestRun(String runId, JsonNode arg) {
        try (Session session = driver.session()) {
            session.run(
                "CREATE (:ARGRun {runId: $rid, source: 'msprime', "
              + "n_samples: $ns, sequence_length: $sl})",
                Map.of("rid", runId,
                       "ns", arg.get("n_samples").asInt(),
                       "sl", arg.get("sequence_length").asLong()));

            for (JsonNode n : arg.get("nodes")) {
                Map<String, Object> p = new HashMap<>();
                p.put("tid", runId + ":" + n.get("id").asInt());
                p.put("rid", runId);
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
                p.put("pid", runId + ":" + e.get("parent").asInt());
                p.put("cid", runId + ":" + e.get("child").asInt());
                p.put("rid", runId);
                p.put("s", e.get("start").asLong());
                p.put("en", e.get("end").asLong());
                session.run(
                    "MATCH (p:TreeNode {treeNodeId: $pid}), "
                  + "(c:TreeNode {treeNodeId: $cid}) "
                  + "CREATE (p)-[:PARENT_OF {runId: $rid, start: $s, end: $en}]->(c)",
                    p);
            }
            for (JsonNode s : arg.get("samples")) {
                int sampleNodeId = s.asInt();
                session.run(
                    "MATCH (n:TreeNode {treeNodeId: $tid}), (s:Sample {sampleId: $sid}) "
                  + "CREATE (n)-[:REPRESENTS {runId: $rid, haplotype: 0}]->(s)",
                    Map.of("tid", runId + ":" + sampleNodeId,
                           "sid", "hap_" + sampleNodeId,
                           "rid", runId));
            }

            JsonNode muts = arg.get("mutations");
            for (int i = 0; i < muts.size(); i++) {
                JsonNode m = muts.get(i);
                long pos = m.get("position").asLong();
                int childNodeId = m.get("child_node_id").asInt();
                String vid = runId + ":chr1:" + pos + ":A:T";
                String gene = (i < MID) ? "G_test" : "G_other";
                String consequence = (i < MID) ? "missense_variant"
                                                : "synonymous_variant";
                Map<String, Object> p = new HashMap<>();
                p.put("vid", vid);
                p.put("pos", pos);
                p.put("rid", runId);
                p.put("tid", runId + ":" + childNodeId);
                p.put("gene", gene);
                p.put("conseq", consequence);
                p.put("derived", m.get("derived_state").asText());
                session.run(
                    "MERGE (v:Variant {variantId: $vid}) "
                  + "ON CREATE SET v.chr = 'chr1', v.pos = $pos, v.ref = 'A', v.alt = 'T' "
                  + "WITH v "
                  + "MATCH (tn:TreeNode {treeNodeId: $tid}) "
                  + "CREATE (v)-[:MUTATED_ON {runId: $rid, parent_node_id: -1, "
                  + "derived_state: $derived}]->(tn) "
                  + "WITH v "
                  + "MATCH (g:Gene {geneId: $gene}) "
                  + "CREATE (v)-[:HAS_CONSEQUENCE {consequence: $conseq, impact: 'MODERATE'}]->(g)",
                    p);
            }
        }
    }

    // ---- helpers --------------------------------------------------------

    private static double[][] runMeanMatrix(List<String> runIds, Map<String, Object> options) {
        try (Session session = driver.session()) {
            List<Record> rows = session.run(
                "CALL graphpop.kinship.branch_grm_posterior($rids, $opts) "
              + "YIELD sample_a, sample_b, b_ij_mean RETURN sample_a, sample_b, b_ij_mean",
                Map.of("rids", runIds, "opts", options)).list();
            double[][] m = new double[nSamples][nSamples];
            for (Record r : rows) {
                int a = parseHap(r.get("sample_a").asString());
                int b = parseHap(r.get("sample_b").asString());
                double v = r.get("b_ij_mean").asDouble();
                m[a][b] = v;
                m[b][a] = v;
            }
            return m;
        }
    }

    private static double[][] runSeMatrix(List<String> runIds, Map<String, Object> options) {
        try (Session session = driver.session()) {
            List<Record> rows = session.run(
                "CALL graphpop.kinship.branch_grm_posterior($rids, $opts) "
              + "YIELD sample_a, sample_b, b_ij_sd RETURN sample_a, sample_b, b_ij_sd",
                Map.of("rids", runIds, "opts", options)).list();
            double[][] m = new double[nSamples][nSamples];
            for (Record r : rows) {
                int a = parseHap(r.get("sample_a").asString());
                int b = parseHap(r.get("sample_b").asString());
                double v = r.get("b_ij_sd").asDouble();
                m[a][b] = v;
                m[b][a] = v;
            }
            return m;
        }
    }

    private static int parseHap(String sid) {
        int colon = sid.indexOf(':');
        String head = (colon >= 0) ? sid.substring(0, colon) : sid;
        return Integer.parseInt(head.substring("hap_".length()));
    }

    private static double[][] loadMatrix(String resourceName) throws Exception {
        JsonNode root = loadJson("/" + resourceName);
        JsonNode mat = root.get("matrix");
        int n = mat.size();
        double[][] m = new double[n][n];
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                m[i][j] = mat.get(i).get(j).asDouble();
        return m;
    }

    private static void assertMatrixAgrees(double[][] expected, double[][] actual,
                                            double relTol, double absTol,
                                            String label) {
        double maxAbs = 0, maxRel = 0;
        int mi = -1, mj = -1;
        for (int i = 0; i < expected.length; i++) {
            for (int j = 0; j < expected.length; j++) {
                double diff = Math.abs(expected[i][j] - actual[i][j]);
                double denom = Math.max(Math.abs(expected[i][j]), 1e-12);
                double rel = diff / denom;
                if (diff > maxAbs) { maxAbs = diff; maxRel = rel; mi = i; mj = j; }
            }
        }
        assertTrue(maxRel < relTol && maxAbs < absTol,
            String.format("%s mismatch (%d,%d): exp=%.10f act=%.10f rel=%.3e abs=%.3e",
                label, mi, mj, expected[mi][mj], actual[mi][mj], maxRel, maxAbs));
    }

    // ---- tests ----------------------------------------------------------

    @Test
    void emptyRunIds_returnsEmpty() {
        try (Session session = driver.session()) {
            long c = session.run(
                "CALL graphpop.kinship.branch_grm_posterior([]) "
              + "YIELD sample_a RETURN count(*) AS c"
            ).single().get("c").asLong();
            assertEquals(0L, c);
        }
    }

    @Test
    void singleRun_seIsNaN_meanIsValid() {
        try (Session session = driver.session()) {
            List<Record> rows = session.run(
                "CALL graphpop.kinship.branch_grm_posterior($rids) "
              + "YIELD sample_a, sample_b, b_ij_mean, b_ij_sd, n_runs",
                Map.of("rids", List.of(RUN_IDS.get(0)))).list();
            assertEquals((long) nSamples * (nSamples + 1) / 2, rows.size());
            for (Record r : rows) {
                assertEquals(1L, r.get("n_runs").asLong());
                assertTrue(Double.isNaN(r.get("b_ij_sd").asDouble()),
                    "single-run sd must be NaN");
            }
        }
    }

    @Test
    void identicalRuns_seIsZero() {
        try (Session session = driver.session()) {
            // The 'identical run' is the SAME runId twice: the procedure
            // sees two runs with identical matrices, so SE should be 0.
            List<Record> rows = session.run(
                "CALL graphpop.kinship.branch_grm_posterior($rids) "
              + "YIELD b_ij_sd, n_runs RETURN b_ij_sd, n_runs",
                Map.of("rids", List.of(RUN_IDS.get(0), RUN_IDS.get(0)))).list();
            for (Record r : rows) {
                assertEquals(2L, r.get("n_runs").asLong());
                assertEquals(0.0, r.get("b_ij_sd").asDouble(), 1e-15);
            }
        }
    }

    @Test
    void multiRunMean_matchesPythonReference() throws Exception {
        double[][] expected = loadMatrix("egrm_posterior_expected_mean.json");
        double[][] actual = runMeanMatrix(RUN_IDS, Map.of());
        assertMatrixAgrees(expected, actual, 1e-6, 1e-9, "posterior-mean unconditional");
    }

    @Test
    void multiRunSe_matchesPythonReference() throws Exception {
        double[][] expected = loadMatrix("egrm_posterior_expected_sd.json");
        double[][] actual = runSeMatrix(RUN_IDS, Map.of());
        assertMatrixAgrees(expected, actual, 1e-6, 1e-9, "posterior-se unconditional");
    }

    @Test
    void multiRunPathwayMean_matchesPythonReference() throws Exception {
        double[][] expected = loadMatrix("egrm_posterior_pathway_expected_mean.json");
        double[][] actual = runMeanMatrix(RUN_IDS,
            Map.of("restrict_to_pathway", "P_test"));
        assertMatrixAgrees(expected, actual, 1e-6, 1e-9,
            "posterior-mean conditional");
    }

    @Test
    void multiRunPathwaySe_matchesPythonReference() throws Exception {
        double[][] expected = loadMatrix("egrm_posterior_pathway_expected_sd.json");
        double[][] actual = runSeMatrix(RUN_IDS,
            Map.of("restrict_to_pathway", "P_test"));
        assertMatrixAgrees(expected, actual, 1e-6, 1e-9,
            "posterior-se conditional");
    }

    @Test
    void sampleAlignmentMismatch_throws() {
        // Ingest a small "rogue" run with a different sample set.
        try (Session session = driver.session()) {
            session.run("CREATE (:ARGRun {runId: 'rogue', source: 'test', "
                      + "n_samples: 1, sequence_length: 100})");
            session.run("CREATE (:Sample {sampleId: 'rogue_hap_0'})");
            session.run(
                "CREATE (:TreeNode {treeNodeId: 'rogue:0', runId: 'rogue', "
              + "nodeId: 0, time: 0.0, is_sample: true, flags: 1}), "
              + "(:TreeNode {treeNodeId: 'rogue:1', runId: 'rogue', "
              + "nodeId: 1, time: 1.0, is_sample: false, flags: 0})");
            session.run(
                "MATCH (p:TreeNode {treeNodeId: 'rogue:1'}), "
              + "(c:TreeNode {treeNodeId: 'rogue:0'}) "
              + "CREATE (p)-[:PARENT_OF {runId: 'rogue', start: 0, end: 100}]->(c)");
            session.run(
                "MATCH (n:TreeNode {treeNodeId: 'rogue:0'}), (s:Sample {sampleId: 'rogue_hap_0'}) "
              + "CREATE (n)-[:REPRESENTS {runId: 'rogue', haplotype: 0}]->(s)");
        }

        try (Session session = driver.session()) {
            Exception ex = assertThrows(Exception.class, () ->
                session.run(
                    "CALL graphpop.kinship.branch_grm_posterior($rids) "
                  + "YIELD sample_a RETURN sample_a",
                    Map.of("rids", List.of(RUN_IDS.get(0), "rogue"))).list()
            );
            String msg = ex.getMessage();
            assertTrue(msg.contains("alignment") || msg.contains("Sample"),
                "exception must mention alignment: " + msg);
        }
    }
}
