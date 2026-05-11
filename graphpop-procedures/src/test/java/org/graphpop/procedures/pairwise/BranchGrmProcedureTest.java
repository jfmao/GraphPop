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

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Cypher integration test for {@link BranchGrmProcedure}.
 *
 * <p><b>Headline correctness gate</b>: the matrix returned by
 * {@code graphpop.kinship.branch_grm} must match the {@code egrm.varGRM}
 * Python reference (checked-in JSON) to relative error &lt; 1×10⁻⁶
 * element-wise.</p>
 *
 * <p>Fixture sources (regenerate via
 * {@code src/test/python/build_egrm_fixture.py}):</p>
 * <ul>
 *   <li>{@code src/test/resources/egrm_fixture_20samples_arg.json} — ARG
 *       topology (nodes, edges, samples, sequence_length).</li>
 *   <li>{@code src/test/resources/egrm_expected_20samples.json} —
 *       reference matrix from {@code egrm.varGRM(ts)}.</li>
 * </ul>
 */
class BranchGrmProcedureTest {

    private static final String RUN_ID = "egrm_fixture_20samples";

    private static Neo4j embeddedNeo4j;
    private static Driver driver;
    private static double[][] expectedMatrix;
    private static int nSamples;

    @BeforeAll
    static void setUp() throws Exception {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(BranchGrmProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());

        // Load the ARG topology and CREATE Neo4j nodes/edges.
        ObjectMapper mapper = new ObjectMapper();
        JsonNode arg;
        try (InputStream is = BranchGrmProcedureTest.class.getResourceAsStream(
                "/egrm_fixture_20samples_arg.json")) {
            assertNotNull(is, "ARG fixture JSON missing on classpath");
            arg = mapper.readTree(is);
        }
        JsonNode expected;
        try (InputStream is = BranchGrmProcedureTest.class.getResourceAsStream(
                "/egrm_expected_20samples.json")) {
            assertNotNull(is, "Expected matrix JSON missing on classpath");
            expected = mapper.readTree(is);
        }

        nSamples = expected.get("n_samples").asInt();
        expectedMatrix = new double[nSamples][nSamples];
        JsonNode mat = expected.get("matrix");
        for (int i = 0; i < nSamples; i++) {
            for (int j = 0; j < nSamples; j++) {
                expectedMatrix[i][j] = mat.get(i).get(j).asDouble();
            }
        }

        try (Session session = driver.session()) {
            session.run(
                "CREATE (:ARGRun {runId: $rid, source: 'msprime', n_samples: $ns, "
              + "sequence_length: $sl})",
                Map.of("rid", RUN_ID,
                       "ns", arg.get("n_samples").asInt(),
                       "sl", arg.get("sequence_length").asLong()));

            // Sample nodes (haplotype level: one :Sample per sample-flagged TreeNode).
            for (JsonNode s : arg.get("samples")) {
                int sampleNodeId = s.asInt();
                session.run(
                    "CREATE (:Sample {sampleId: $sid, packed_index: $pi})",
                    Map.of("sid", sampleId(sampleNodeId), "pi", sampleNodeId));
            }

            // TreeNodes.
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

            // PARENT_OF edges.
            for (JsonNode e : arg.get("edges")) {
                Map<String, Object> p = new HashMap<>();
                p.put("pid", RUN_ID + ":" + e.get("parent").asInt());
                p.put("cid", RUN_ID + ":" + e.get("child").asInt());
                p.put("rid", RUN_ID);
                p.put("s", e.get("start").asLong());
                p.put("en", e.get("end").asLong());
                session.run(
                    "MATCH (p:TreeNode {treeNodeId: $pid}), (c:TreeNode {treeNodeId: $cid}) "
                  + "CREATE (p)-[:PARENT_OF {runId: $rid, start: $s, end: $en}]->(c)", p);
            }

            // REPRESENTS edges (sample TreeNode -> :Sample).
            for (JsonNode s : arg.get("samples")) {
                int sampleNodeId = s.asInt();
                session.run(
                    "MATCH (n:TreeNode {treeNodeId: $tid}), (s:Sample {sampleId: $sid}) "
                  + "CREATE (n)-[:REPRESENTS {runId: $rid, haplotype: 0}]->(s)",
                    Map.of("tid", RUN_ID + ":" + sampleNodeId,
                           "sid", sampleId(sampleNodeId),
                           "rid", RUN_ID));
            }

            // ---- Synthetic annotation layer for step-4 conditional tests ----
            // First MID=20 mutations are "lit" (in pathway P_test, consequence
            // missense). The remaining 21 are "dark" (in G_other, synonymous).
            session.run("CREATE (:Pathway {pathwayId: 'P_test', name: 'test pathway'})");
            session.run("CREATE (:Pathway {pathwayId: 'P_other', name: 'other pathway'})");
            session.run(
                "MATCH (p:Pathway {pathwayId: 'P_test'}) "
              + "CREATE (g:Gene {geneId: 'G_test', symbol: 'GTEST'})-[:IN_PATHWAY]->(p)");
            session.run(
                "MATCH (p:Pathway {pathwayId: 'P_other'}) "
              + "CREATE (g:Gene {geneId: 'G_other', symbol: 'GOTHER'})-[:IN_PATHWAY]->(p)");

            JsonNode muts = arg.get("mutations");
            for (int i = 0; i < muts.size(); i++) {
                JsonNode m = muts.get(i);
                long pos = m.get("position").asLong();
                int childNodeId = m.get("child_node_id").asInt();
                String vid = "chr1:" + pos + ":A:T";
                String gene = (i < MID) ? "G_test" : "G_other";
                String consequence = (i < MID) ? "missense_variant"
                                                : "synonymous_variant";
                Map<String, Object> p = new HashMap<>();
                p.put("vid", vid);
                p.put("pos", pos);
                p.put("rid", RUN_ID);
                p.put("tid", RUN_ID + ":" + childNodeId);
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

    /**
     * Index of the first {@code MID} mutations that are tagged as
     * &quot;lit&quot; (in pathway {@code P_test} and consequence
     * {@code missense_variant}).  Must match the constant in
     * {@code build_egrm_fixture.py}.
     */
    private static final int MID = 20;

    @AfterAll
    static void tearDown() {
        if (driver != null) driver.close();
        if (embeddedNeo4j != null) embeddedNeo4j.close();
    }

    private static String sampleId(int tskitSampleNodeId) {
        return "hap_" + tskitSampleNodeId;
    }

    private static double[][] runProcedureMatrix() {
        return runProcedureMatrix(Map.of());
    }

    private static double[][] runProcedureMatrix(Map<String, Object> options) {
        try (Session session = driver.session()) {
            List<Record> rows = session.run(
                "CALL graphpop.kinship.branch_grm($rid, $opts) "
              + "YIELD sample_a, sample_b, phi RETURN sample_a, sample_b, phi",
                Map.of("rid", RUN_ID, "opts", options)).list();

            double[][] m = new double[nSamples][nSamples];
            for (Record r : rows) {
                int a = parseHapIndex(r.get("sample_a").asString());
                int b = parseHapIndex(r.get("sample_b").asString());
                double v = r.get("phi").asDouble();
                m[a][b] = v;
                m[b][a] = v;  // procedure emits upper-triangular only
            }
            return m;
        }
    }

    private static double[][] loadExpectedMatrix(String resourceName) throws Exception {
        ObjectMapper mapper = new ObjectMapper();
        try (InputStream is = BranchGrmProcedureTest.class.getResourceAsStream(
                "/" + resourceName)) {
            assertNotNull(is, "missing resource: " + resourceName);
            JsonNode root = mapper.readTree(is);
            JsonNode mat = root.get("matrix");
            int n = mat.size();
            double[][] m = new double[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    m[i][j] = mat.get(i).get(j).asDouble();
                }
            }
            return m;
        }
    }

    private static void assertMatricesAgree(double[][] expected, double[][] actual,
                                            double relTol, double absTol,
                                            String label) {
        double maxAbs = 0;
        double maxRel = 0;
        int maxI = -1, maxJ = -1;
        for (int i = 0; i < expected.length; i++) {
            for (int j = 0; j < expected.length; j++) {
                double diff = Math.abs(expected[i][j] - actual[i][j]);
                double denom = Math.max(Math.abs(expected[i][j]), 1e-12);
                double rel = diff / denom;
                if (diff > maxAbs) {
                    maxAbs = diff; maxRel = rel; maxI = i; maxJ = j;
                }
            }
        }
        assertTrue(maxRel < relTol && maxAbs < absTol,
            String.format(
                "%s mismatch at (%d,%d): expected=%.10f actual=%.10f rel=%.3e abs=%.3e",
                label, maxI, maxJ,
                expected[maxI][maxJ], actual[maxI][maxJ], maxRel, maxAbs));
    }

    private static double readJsonDouble(String resourceName, String key) throws Exception {
        ObjectMapper mapper = new ObjectMapper();
        try (InputStream is = BranchGrmProcedureTest.class.getResourceAsStream(
                "/" + resourceName)) {
            return mapper.readTree(is).get(key).asDouble();
        }
    }

    private static int parseHapIndex(String sampleId) {
        // sampleId format: "hap_<i>:h<haplotype>" or "hap_<i>"
        int colon = sampleId.indexOf(':');
        String head = (colon >= 0) ? sampleId.substring(0, colon) : sampleId;
        return Integer.parseInt(head.substring("hap_".length()));
    }

    @Test
    void branchGrm_matchesEgrmReference() {
        double[][] actual = runProcedureMatrix();

        // Element-wise rel-err comparison.
        double maxAbsDiff = 0.0;
        double maxRelErr = 0.0;
        int maxI = -1, maxJ = -1;
        for (int i = 0; i < nSamples; i++) {
            for (int j = 0; j < nSamples; j++) {
                double exp = expectedMatrix[i][j];
                double got = actual[i][j];
                double absDiff = Math.abs(exp - got);
                double denom = Math.max(Math.abs(exp), 1e-12);
                double relErr = absDiff / denom;
                if (absDiff > maxAbsDiff) {
                    maxAbsDiff = absDiff;
                    maxRelErr = relErr;
                    maxI = i;
                    maxJ = j;
                }
            }
        }
        assertTrue(maxRelErr < 1e-6,
            String.format(
                "egrm mismatch at (%d, %d): expected=%.10f actual=%.10f rel_err=%.3e",
                maxI, maxJ, expectedMatrix[maxI][maxJ], 0.0, maxRelErr));
        assertTrue(maxAbsDiff < 1e-9,
            String.format("max abs diff %.3e too large (expected < 1e-9)", maxAbsDiff));

        // Opt-in: dump procedure output for Paper 2 Fig 1d driver.
        maybeDumpTsv(actual, "branch_grm_unconditional.tsv");
    }

    /**
     * Optional H1-schema TSV dump of a procedure matrix. Triggered
     * only when the system property {@code graphpop.bench.dump.dir}
     * is set (typically by the Paper 2 Fig 1d/1e driver
     * {@code paper/paper2_kinship_arg/benchmarks/run_fig1de.py}).
     *
     * <p>Schema: {@code sample_a, sample_b, kinship}, upper triangle
     * including diagonal. Sample IDs are the 0-based haplotype node
     * IDs serialised as strings, matching the H4 wrapper output and
     * the {@code egrm_expected_*.json} fixture sample_ids.</p>
     */
    private static void maybeDumpTsv(double[][] matrix, String filename) {
        String dir = System.getProperty("graphpop.bench.dump.dir");
        if (dir == null || dir.isBlank()) return;
        Path out = Paths.get(dir, filename);
        try {
            if (out.getParent() != null) {
                Files.createDirectories(out.getParent());
            }
            try (BufferedWriter w = Files.newBufferedWriter(out)) {
                w.write("sample_a\tsample_b\tkinship\n");
                int n = matrix.length;
                for (int i = 0; i < n; i++) {
                    for (int j = i; j < n; j++) {
                        w.write(String.format(
                            Locale.ROOT, "%d\t%d\t%.10g%n",
                            i, j, matrix[i][j]));
                    }
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(
                "failed to dump branch GRM TSV to " + out, e);
        }
    }

    @Test
    void branchGrm_isSymmetric() {
        double[][] actual = runProcedureMatrix();
        for (int i = 0; i < nSamples; i++) {
            for (int j = i + 1; j < nSamples; j++) {
                assertEquals(actual[i][j], actual[j][i], 1e-12,
                    "asymmetric at (" + i + "," + j + ")");
            }
        }
    }

    @Test
    void branchGrm_rowSumsZeroAfterCentering() {
        double[][] actual = runProcedureMatrix();
        for (int i = 0; i < nSamples; i++) {
            double sum = 0;
            for (int j = 0; j < nSamples; j++) sum += actual[i][j];
            assertEquals(0.0, sum, 1e-9, "row " + i + " sum non-zero");
        }
    }

    @Test
    void branchGrm_emitsSelfPairsByDefault() {
        try (Session session = driver.session()) {
            long selfCount = session.run(
                "CALL graphpop.kinship.branch_grm($rid) "
              + "YIELD sample_a, sample_b "
              + "WHERE sample_a = sample_b RETURN count(*) AS c",
                Map.of("rid", RUN_ID)).single().get("c").asLong();
            assertEquals(nSamples, selfCount,
                "default include_self=true must emit n self-pairs");
        }
    }

    @Test
    void branchGrm_unknownRunReturnsEmpty() {
        try (Session session = driver.session()) {
            long c = session.run(
                "CALL graphpop.kinship.branch_grm('does_not_exist') "
              + "YIELD sample_a RETURN count(*) AS c"
            ).single().get("c").asLong();
            assertEquals(0L, c);
        }
    }

    @Test
    void branchGrm_regionRestrictionSubsetEqualsFullWhenWidened() {
        // The fixture is sequence_length=50000. Calling with start=0, end=50000
        // must agree with the default (full-region) output. The full-region
        // case is already validated against egrm; this exercises explicit
        // bounds.
        try (Session session = driver.session()) {
            List<Record> full = session.run(
                "CALL graphpop.kinship.branch_grm($rid) "
              + "YIELD sample_a, sample_b, phi RETURN sample_a, sample_b, phi",
                Map.of("rid", RUN_ID)).list();
            List<Record> bounded = session.run(
                "CALL graphpop.kinship.branch_grm($rid, {start: 0, end: 50000}) "
              + "YIELD sample_a, sample_b, phi RETURN sample_a, sample_b, phi",
                Map.of("rid", RUN_ID)).list();
            assertEquals(full.size(), bounded.size());

            // Build maps and compare element-wise.
            Map<String, Double> fm = new HashMap<>(), bm = new HashMap<>();
            for (Record r : full) {
                fm.put(r.get("sample_a").asString() + "|" + r.get("sample_b").asString(),
                       r.get("phi").asDouble());
            }
            for (Record r : bounded) {
                bm.put(r.get("sample_a").asString() + "|" + r.get("sample_b").asString(),
                       r.get("phi").asDouble());
            }
            for (var entry : fm.entrySet()) {
                assertEquals(entry.getValue(), bm.get(entry.getKey()), 1e-12,
                    "mismatch for " + entry.getKey());
            }
        }
    }

    @Test
    void branchGrm_emitsExpectedNumberOfPairs() {
        // n samples, include_self=true => n*(n+1)/2 unique pairs.
        try (Session session = driver.session()) {
            long c = session.run(
                "CALL graphpop.kinship.branch_grm($rid) YIELD sample_a "
              + "RETURN count(*) AS c", Map.of("rid", RUN_ID)
            ).single().get("c").asLong();
            assertEquals((long) nSamples * (nSamples + 1) / 2, c);
        }
    }

    // ---- Step 4 conditional predicates ---------------------------------

    @Test
    void timeWindow_matchesReference() throws Exception {
        double tLo = readJsonDouble("egrm_expected_time_window.json", "t_lo");
        double tHi = readJsonDouble("egrm_expected_time_window.json", "t_hi");
        double[][] expected = loadExpectedMatrix("egrm_expected_time_window.json");
        double[][] actual = runProcedureMatrix(Map.of(
            "time_window", List.of(tLo, tHi)));
        assertMatricesAgree(expected, actual, 1e-6, 1e-9, "time_window");
    }

    @Test
    void pathway_matchesReference() throws Exception {
        double[][] expected = loadExpectedMatrix("egrm_expected_pathway_half.json");
        double[][] actual = runProcedureMatrix(Map.of(
            "restrict_to_pathway", "P_test"));
        assertMatricesAgree(expected, actual, 1e-6, 1e-9, "restrict_to_pathway");

        // Opt-in: dump procedure output for Paper 2 Fig 1e driver.
        maybeDumpTsv(actual, "branch_grm_pathway_half.tsv");
    }

    @Test
    void mutationFilter_matchesReference() throws Exception {
        double[][] expected = loadExpectedMatrix(
            "egrm_expected_consequence_missense.json");
        double[][] actual = runProcedureMatrix(Map.of(
            "mutation_filter", "missense_variant"));
        assertMatricesAgree(expected, actual, 1e-6, 1e-9, "mutation_filter");
    }

    @Test
    void composition_pathwayAndTimeWindow_matchesReference() throws Exception {
        double tLo = readJsonDouble(
            "egrm_expected_composed_pathway_time.json", "t_lo");
        double tHi = readJsonDouble(
            "egrm_expected_composed_pathway_time.json", "t_hi");
        double[][] expected = loadExpectedMatrix(
            "egrm_expected_composed_pathway_time.json");
        Map<String, Object> opts = new HashMap<>();
        opts.put("restrict_to_pathway", "P_test");
        opts.put("time_window", List.of(tLo, tHi));
        double[][] actual = runProcedureMatrix(opts);
        assertMatricesAgree(expected, actual, 1e-6, 1e-9, "composed");
    }

    @Test
    void pathway_unknownReturnsZeroMatrix() {
        double[][] m = runProcedureMatrix(Map.of(
            "restrict_to_pathway", "P_does_not_exist"));
        for (int i = 0; i < nSamples; i++)
            for (int j = 0; j < nSamples; j++)
                assertEquals(0.0, m[i][j], 1e-12,
                    "non-existent pathway must yield zero entry [" + i + "," + j + "]");
    }

    @Test
    void mutationFilter_unknownReturnsZeroMatrix() {
        double[][] m = runProcedureMatrix(Map.of(
            "mutation_filter", "frameshift_variant"));
        for (int i = 0; i < nSamples; i++)
            for (int j = 0; j < nSamples; j++)
                assertEquals(0.0, m[i][j], 1e-12,
                    "non-existent consequence must yield zero entry [" + i + "," + j + "]");
    }
}
