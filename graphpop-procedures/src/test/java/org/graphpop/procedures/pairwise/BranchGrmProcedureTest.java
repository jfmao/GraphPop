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
        }
    }

    @AfterAll
    static void tearDown() {
        if (driver != null) driver.close();
        if (embeddedNeo4j != null) embeddedNeo4j.close();
    }

    private static String sampleId(int tskitSampleNodeId) {
        return "hap_" + tskitSampleNodeId;
    }

    private static double[][] runProcedureMatrix() {
        try (Session session = driver.session()) {
            List<Record> rows = session.run(
                "CALL graphpop.kinship.branch_grm($rid) "
              + "YIELD sample_a, sample_b, phi RETURN sample_a, sample_b, phi",
                Map.of("rid", RUN_ID)).list();

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
}
