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
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Integration test for {@link BranchGrmByAncestryProcedure}.
 *
 * <p>Headline correctness gates: the per-ancestry sub-matrices match
 * the Python reference (built by {@code build_egrm_fixture.py}) to
 * relative error &lt; 10⁻⁶, and their sum reproduces the
 * unconditional {@code branch_grm} matrix to relative error &lt;
 * 10⁻⁹.</p>
 */
class BranchGrmByAncestryProcedureTest {

    private static final String RUN_ID = "egrm_fixture_20samples";
    private static final int MID = 20;

    private static Neo4j embeddedNeo4j;
    private static Driver driver;
    private static int nSamples;

    @BeforeAll
    static void setUp() throws Exception {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(BranchGrmByAncestryProcedure.class)
                .withProcedure(BranchGrmProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());

        ObjectMapper mapper = new ObjectMapper();
        JsonNode arg = loadJson("/egrm_fixture_20samples_arg.json");
        JsonNode painting = loadJson("/egrm_by_ancestry_fixture_painting.json");
        nSamples = arg.get("n_samples").asInt();

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

            // Annotation layer for pathway-restricted decomposition test.
            session.run("CREATE (:Pathway {pathwayId: 'P_test', name: 'test'})");
            session.run(
                "MATCH (p:Pathway {pathwayId: 'P_test'}) "
              + "CREATE (g:Gene {geneId: 'G_test'})-[:IN_PATHWAY]->(p)");
            session.run("CREATE (g:Gene {geneId: 'G_other'})");

            JsonNode muts = arg.get("mutations");
            for (int i = 0; i < muts.size(); i++) {
                JsonNode m = muts.get(i);
                long pos = m.get("position").asLong();
                int childNodeId = m.get("child_node_id").asInt();
                String vid = "chr1:" + pos + ":A:T";
                String gene = (i < MID) ? "G_test" : "G_other";
                String consequence = (i < MID) ? "missense_variant" : "synonymous_variant";
                Map<String, Object> p = new HashMap<>();
                p.put("vid", vid); p.put("pos", pos); p.put("rid", RUN_ID);
                p.put("tid", RUN_ID + ":" + childNodeId);
                p.put("gene", gene); p.put("conseq", consequence);
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

            // Ingest the painting from the JSON fixture.
            Set<String> populations = new HashSet<>();
            for (JsonNode row : painting.get("rows")) {
                populations.add(row.get("population_id").asText());
            }
            for (String pop : populations) {
                session.run("CREATE (:Population {populationId: $pop})",
                    Map.of("pop", pop));
            }
            for (JsonNode row : painting.get("rows")) {
                Map<String, Object> p = new HashMap<>();
                p.put("tid", RUN_ID + ":" + row.get("tskit_node_id").asInt());
                p.put("pop", row.get("population_id").asText());
                p.put("rid", RUN_ID);
                p.put("prob", row.get("posterior_prob").asDouble());
                p.put("painter", "majority_vote");
                session.run(
                    "MATCH (n:TreeNode {treeNodeId: $tid}), (pop:Population {populationId: $pop}) "
                  + "CREATE (n)-[:HAS_ANCESTRY {runId: $rid, posterior_prob: $prob, "
                  + "painter: $painter}]->(pop)", p);
            }
        }
    }

    @AfterAll
    static void tearDown() {
        if (driver != null) driver.close();
        if (embeddedNeo4j != null) embeddedNeo4j.close();
    }

    private static JsonNode loadJson(String resource) throws Exception {
        ObjectMapper mapper = new ObjectMapper();
        try (InputStream is = BranchGrmByAncestryProcedureTest.class
                .getResourceAsStream(resource)) {
            assertNotNull(is, "missing resource: " + resource);
            return mapper.readTree(is);
        }
    }

    private static double[][] runMatrixForAncestry(Map<String, Object> options,
                                                    String ancestry) {
        try (Session session = driver.session()) {
            List<Record> rows = session.run(
                "CALL graphpop.kinship.branch_grm_by_ancestry($rid, $opts) "
              + "YIELD sample_a, sample_b, ancestry, b_ij_component "
              + "WHERE ancestry = $anc "
              + "RETURN sample_a, sample_b, b_ij_component",
                Map.of("rid", RUN_ID, "opts", options, "anc", ancestry)).list();
            double[][] m = new double[nSamples][nSamples];
            for (Record r : rows) {
                int a = parseHap(r.get("sample_a").asString());
                int b = parseHap(r.get("sample_b").asString());
                double v = r.get("b_ij_component").asDouble();
                m[a][b] = v;
                m[b][a] = v;
            }
            return m;
        }
    }

    private static double[][] runUnconditionalMatrix() {
        try (Session session = driver.session()) {
            List<Record> rows = session.run(
                "CALL graphpop.kinship.branch_grm($rid) "
              + "YIELD sample_a, sample_b, phi RETURN sample_a, sample_b, phi",
                Map.of("rid", RUN_ID)).list();
            double[][] m = new double[nSamples][nSamples];
            for (Record r : rows) {
                int a = parseHap(r.get("sample_a").asString());
                int b = parseHap(r.get("sample_b").asString());
                double v = r.get("phi").asDouble();
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

    @Test
    void unconditional_eur_matchesPythonReference() throws Exception {
        double[][] expected = loadMatrix("egrm_by_ancestry_expected_EUR.json");
        double[][] actual = runMatrixForAncestry(Map.of(), "EUR");
        assertMatrixAgrees(expected, actual, 1e-6, 1e-9, "decompose-EUR");
    }

    @Test
    void unconditional_afr_matchesPythonReference() throws Exception {
        double[][] expected = loadMatrix("egrm_by_ancestry_expected_AFR.json");
        double[][] actual = runMatrixForAncestry(Map.of(), "AFR");
        assertMatrixAgrees(expected, actual, 1e-6, 1e-9, "decompose-AFR");
    }

    @Test
    void partition_sumOfAncestriesEqualsUnconditional() {
        double[][] eur = runMatrixForAncestry(Map.of(), "EUR");
        double[][] afr = runMatrixForAncestry(Map.of(), "AFR");
        double[][] uncond = runUnconditionalMatrix();
        double[][] sum = new double[nSamples][nSamples];
        for (int i = 0; i < nSamples; i++)
            for (int j = 0; j < nSamples; j++)
                sum[i][j] = eur[i][j] + afr[i][j];
        assertMatrixAgrees(uncond, sum, 1e-6, 1e-9, "Σ_a B_ij^a == B_ij^uncond");
    }

    @Test
    void pathway_eur_matchesPythonReference() throws Exception {
        double[][] expected = loadMatrix("egrm_by_ancestry_pathway_expected_EUR.json");
        double[][] actual = runMatrixForAncestry(
            Map.of("restrict_to_pathway", "P_test"), "EUR");
        assertMatrixAgrees(expected, actual, 1e-6, 1e-9, "pathway-decompose-EUR");
    }

    @Test
    void pathway_afr_matchesPythonReference() throws Exception {
        double[][] expected = loadMatrix("egrm_by_ancestry_pathway_expected_AFR.json");
        double[][] actual = runMatrixForAncestry(
            Map.of("restrict_to_pathway", "P_test"), "AFR");
        assertMatrixAgrees(expected, actual, 1e-6, 1e-9, "pathway-decompose-AFR");
    }

    @Test
    void unpaintedRun_returnsEmpty() {
        try (Session session = driver.session()) {
            // Create a separate run without painting.
            session.run("CREATE (:ARGRun {runId: 'unpainted', source: 'test', "
                      + "n_samples: 1, sequence_length: 10}), "
                      + "(:TreeNode {treeNodeId: 'unpainted:0', runId: 'unpainted', "
                      + "nodeId: 0, time: 0.0, is_sample: true, flags: 1}), "
                      + "(:TreeNode {treeNodeId: 'unpainted:1', runId: 'unpainted', "
                      + "nodeId: 1, time: 1.0, is_sample: false, flags: 0})");
            session.run(
                "MATCH (p:TreeNode {treeNodeId: 'unpainted:1'}), "
              + "(c:TreeNode {treeNodeId: 'unpainted:0'}) "
              + "CREATE (p)-[:PARENT_OF {runId: 'unpainted', start: 0, end: 10}]->(c)");

            long c = session.run(
                "CALL graphpop.kinship.branch_grm_by_ancestry('unpainted') "
              + "YIELD sample_a RETURN count(*) AS c"
            ).single().get("c").asLong();
            assertEquals(0L, c);
        }
    }

    @Test
    void emitsRowsForEveryAncestryAndPair() {
        try (Session session = driver.session()) {
            // n*(n+1)/2 pairs × 2 ancestries (EUR, AFR)
            long c = session.run(
                "CALL graphpop.kinship.branch_grm_by_ancestry($rid) "
              + "YIELD sample_a RETURN count(*) AS c",
                Map.of("rid", RUN_ID)).single().get("c").asLong();
            long expected = (long) nSamples * (nSamples + 1) / 2 * 2;
            assertEquals(expected, c);
        }
    }
}
