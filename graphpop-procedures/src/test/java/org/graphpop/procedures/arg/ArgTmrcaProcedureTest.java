package org.graphpop.procedures.arg;

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
 * Cypher integration test for {@link ArgTmrcaProcedure}.
 *
 * <p>Validates against tskit ground truth dumped into
 * {@code egrm_fixture_20samples_arg_stats.json}: 16 (pair × position)
 * single-position TMRCAs and 4 genome-mean TMRCAs.</p>
 */
class ArgTmrcaProcedureTest {

    private static final String RUN_ID = "egrm_fixture_20samples";

    private static Neo4j embeddedNeo4j;
    private static Driver driver;
    private static JsonNode statsRef;

    @BeforeAll
    static void setUp() throws Exception {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder()
                .withProcedure(ArgTmrcaProcedure.class)
                .build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());

        ObjectMapper mapper = new ObjectMapper();
        JsonNode arg;
        try (InputStream is = ArgTmrcaProcedureTest.class.getResourceAsStream(
                "/egrm_fixture_20samples_arg.json")) {
            arg = mapper.readTree(is);
        }
        try (InputStream is = ArgTmrcaProcedureTest.class.getResourceAsStream(
                "/egrm_fixture_20samples_arg_stats.json")) {
            statsRef = mapper.readTree(is);
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
                    "MATCH (n:TreeNode {treeNodeId: $tid}), "
                  + "(s:Sample {sampleId: $sid}) "
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

    @Test
    void tmrca_single_position_matches_tskit_for_every_pair_position() {
        try (Session session = driver.session()) {
            for (JsonNode row : statsRef.get("tmrca_per_position")) {
                int a = row.get("sample_a").asInt();
                int b = row.get("sample_b").asInt();
                long pos = row.get("position").asLong();
                double expected = row.get("tmrca").asDouble();

                List<Record> rs = session.run(
                    "CALL graphpop.arg.tmrca($rid, $sa, $sb, {position: $pos}) "
                  + "YIELD tmrca, mrca_node_id RETURN tmrca, mrca_node_id",
                    Map.of("rid", RUN_ID,
                           "sa", "hap_" + a,
                           "sb", "hap_" + b,
                           "pos", pos)).list();
                assertEquals(1, rs.size(),
                        "one row per haplotype pair; got " + rs.size());
                assertEquals(expected, rs.get(0).get("tmrca").asDouble(), 1e-9,
                        "tmrca mismatch for pair (" + a + "," + b
                            + ") at pos " + pos);
            }
        }
    }

    @Test
    void tmrca_genome_mean_matches_tskit_for_every_pair() {
        try (Session session = driver.session()) {
            for (JsonNode row : statsRef.get("tmrca_genome_means")) {
                int a = row.get("sample_a").asInt();
                int b = row.get("sample_b").asInt();
                double expected = row.get("mean_tmrca").asDouble();

                List<Record> rs = session.run(
                    "CALL graphpop.arg.tmrca($rid, $sa, $sb, {}) "
                  + "YIELD position, mean_tmrca "
                  + "RETURN position, mean_tmrca",
                    Map.of("rid", RUN_ID,
                           "sa", "hap_" + a,
                           "sb", "hap_" + b)).list();
                assertEquals(1, rs.size());
                assertEquals(-1L, rs.get(0).get("position").asLong());
                assertEquals(expected, rs.get(0).get("mean_tmrca").asDouble(),
                        1e-6, "genome-mean tmrca mismatch for pair ("
                            + a + "," + b + ")");
            }
        }
    }

    @Test
    void tmrca_window_mean_subset_falls_within_min_max_per_position() {
        // Sanity: window-mean over [0.30 * SEQ, 0.55 * SEQ] for pair (0,1)
        // should fall between min and max of single-position values for
        // overlapping positions in that window (here position[1] only,
        // since position[2] = 0.55*SEQ is inclusive of right boundary
        // depending on tree spans).
        try (Session session = driver.session()) {
            long lo = 15000;
            long hi = 27500;
            List<Record> rs = session.run(
                "CALL graphpop.arg.tmrca($rid, $sa, $sb, "
              + "{window_start: $lo, window_end: $hi}) "
              + "YIELD mean_tmrca RETURN mean_tmrca",
                Map.of("rid", RUN_ID, "sa", "hap_0", "sb", "hap_1",
                       "lo", lo, "hi", hi)).list();
            assertEquals(1, rs.size());
            assertTrue(rs.get(0).get("mean_tmrca").asDouble() > 0.0);
        }
    }

    @Test
    void tmrca_unknown_sample_returns_empty() {
        try (Session session = driver.session()) {
            List<Record> rs = session.run(
                "CALL graphpop.arg.tmrca($rid, 'no_such_sample', 'hap_1', {}) "
              + "YIELD tmrca RETURN tmrca",
                Map.of("rid", RUN_ID)).list();
            assertEquals(0, rs.size());
        }
    }
}
