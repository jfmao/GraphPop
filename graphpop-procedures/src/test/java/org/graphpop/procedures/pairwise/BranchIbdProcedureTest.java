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
 * Cypher integration test for {@link BranchIbdProcedure}.
 *
 * <p><b>Headline correctness gate</b>: the procedure-emitted IBD
 * segment set on the 20-sample fixture must equal the
 * {@code tskit.TreeSequence.ibd_segments} reference set element-wise.</p>
 */
class BranchIbdProcedureTest {

    private static final String RUN_ID = "egrm_fixture_20samples";

    private static Neo4j embeddedNeo4j;
    private static Driver driver;
    private static Set<String> expectedSegmentKeys;

    @BeforeAll
    static void setUp() throws Exception {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(BranchIbdProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());

        ObjectMapper mapper = new ObjectMapper();
        JsonNode arg;
        try (InputStream is = BranchIbdProcedureTest.class.getResourceAsStream(
                "/egrm_fixture_20samples_arg.json")) {
            arg = mapper.readTree(is);
        }
        JsonNode ibd;
        try (InputStream is = BranchIbdProcedureTest.class.getResourceAsStream(
                "/egrm_fixture_20samples_ibd.json")) {
            ibd = mapper.readTree(is);
        }

        // Build the reference key set: "sa,sb,start,end,mrca".
        expectedSegmentKeys = new HashSet<>();
        for (JsonNode s : ibd.get("segments")) {
            expectedSegmentKeys.add(segmentKey(
                    s.get("sample_a").asInt(),
                    s.get("sample_b").asInt(),
                    s.get("start").asLong(),
                    s.get("end").asLong(),
                    s.get("mrca_node_id").asInt()));
        }

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

    private static String segmentKey(int sa, int sb, long start, long end, int mrca) {
        return sa + "," + sb + "," + start + "," + end + "," + mrca;
    }

    private static int parseHap(String sid) {
        int colon = sid.indexOf(':');
        String head = (colon >= 0) ? sid.substring(0, colon) : sid;
        return Integer.parseInt(head.substring("hap_".length()));
    }

    @Test
    void fromArg_matchesTskitReference() {
        try (Session session = driver.session()) {
            List<Record> rows = session.run(
                "CALL graphpop.ibd.from_arg($rid) "
              + "YIELD sample_a, sample_b, start, end, mrca_node_id "
              + "RETURN sample_a, sample_b, start, end, mrca_node_id",
                Map.of("rid", RUN_ID)).list();

            Set<String> actual = new HashSet<>();
            for (Record r : rows) {
                actual.add(segmentKey(
                        parseHap(r.get("sample_a").asString()),
                        parseHap(r.get("sample_b").asString()),
                        r.get("start").asLong(),
                        r.get("end").asLong(),
                        (int) r.get("mrca_node_id").asLong()));
            }
            assertEquals(expectedSegmentKeys, actual,
                "Java IBD segment set must equal tskit reference exactly. "
              + "actual size=" + actual.size()
              + ", expected size=" + expectedSegmentKeys.size());
        }
    }

    @Test
    void fromArg_writesIbdEdges() {
        try (Session session = driver.session()) {
            // Run again to verify idempotent re-write.
            session.run("CALL graphpop.ibd.from_arg($rid) YIELD sample_a "
                      + "RETURN count(*) AS c", Map.of("rid", RUN_ID))
                    .single();
            long edgeCount = session.run(
                "MATCH ()-[r:IBD_SEGMENT {runId: $rid, source: 'arg_derived'}]->() "
              + "RETURN count(r) AS c",
                Map.of("rid", RUN_ID)).single().get("c").asLong();
            assertEquals(expectedSegmentKeys.size(), edgeCount,
                "Should write one :IBD_SEGMENT per emitted segment");
        }
    }

    @Test
    void fromArg_idempotentOnReRun() {
        try (Session session = driver.session()) {
            session.run("CALL graphpop.ibd.from_arg($rid) YIELD sample_a "
                      + "RETURN count(*) AS c", Map.of("rid", RUN_ID)).single();
            session.run("CALL graphpop.ibd.from_arg($rid) YIELD sample_a "
                      + "RETURN count(*) AS c", Map.of("rid", RUN_ID)).single();
            long edgeCount = session.run(
                "MATCH ()-[r:IBD_SEGMENT {runId: $rid, source: 'arg_derived'}]->() "
              + "RETURN count(r) AS c",
                Map.of("rid", RUN_ID)).single().get("c").asLong();
            assertEquals(expectedSegmentKeys.size(), edgeCount,
                "Idempotent re-run must not duplicate edges");
        }
    }

    @Test
    void fromArg_tmrcaCapDropsLongerSegments() {
        try (Session session = driver.session()) {
            // The fixture's median internal time was ~0.327 (see the
            // egrm_expected_time_window.json). Cap at 0.2 should drop a
            // strict subset.
            List<Record> all = session.run(
                "CALL graphpop.ibd.from_arg($rid) YIELD tmrca RETURN tmrca",
                Map.of("rid", RUN_ID)).list();
            List<Record> capped = session.run(
                "CALL graphpop.ibd.from_arg($rid, {max_tmrca: 0.2}) "
              + "YIELD tmrca RETURN tmrca",
                Map.of("rid", RUN_ID)).list();
            assertTrue(capped.size() < all.size(),
                "TMRCA cap must drop some segments");
            for (Record r : capped) {
                assertTrue(r.get("tmrca").asDouble() <= 0.2 + 1e-12,
                    "all kept segments must have tmrca <= 0.2");
            }
        }
    }

    @Test
    void fromArg_minLengthFilter() {
        try (Session session = driver.session()) {
            List<Record> rows = session.run(
                "CALL graphpop.ibd.from_arg($rid, {min_length_bp: 1000}) "
              + "YIELD length_bp RETURN length_bp",
                Map.of("rid", RUN_ID)).list();
            for (Record r : rows) {
                assertTrue(r.get("length_bp").asLong() >= 1000,
                    "min_length_bp must be enforced");
            }
        }
    }

    @Test
    void fromArg_unknownRunReturnsEmpty() {
        try (Session session = driver.session()) {
            long c = session.run(
                "CALL graphpop.ibd.from_arg('does_not_exist') "
              + "YIELD sample_a RETURN count(*) AS c"
            ).single().get("c").asLong();
            assertEquals(0L, c);
        }
    }
}
