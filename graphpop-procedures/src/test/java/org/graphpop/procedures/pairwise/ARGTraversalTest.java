package org.graphpop.procedures.pairwise;

import org.junit.jupiter.api.AfterAll;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;
import org.neo4j.driver.Driver;
import org.neo4j.driver.GraphDatabase;
import org.neo4j.driver.Session;
import org.neo4j.dbms.api.DatabaseManagementService;
import org.neo4j.graphdb.GraphDatabaseService;
import org.neo4j.graphdb.Transaction;
import org.neo4j.harness.Neo4j;
import org.neo4j.harness.Neo4jBuilders;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Integration test for {@link ARGTraversal} against a hand-built
 * {@code :TreeNode}/{@code :PARENT_OF} fixture.
 *
 * <pre>
 * Toy tree (run "r0", interval [0, 100]):
 *
 *   tskit ids:    4
 *                / \
 *              3    \
 *              / \   \
 *             0   1   2     (samples; time 0)
 *
 *   parent_of edges (parent -> child, [start, end]):
 *     4 -> 3  [0, 100]
 *     3 -> 0  [0, 100]
 *     3 -> 1  [0, 100]
 *     4 -> 2  [0, 100]
 *
 *   times: 0,1,2 -> 0;  3 -> 0.5;  4 -> 1.0
 * </pre>
 */
class ARGTraversalTest {

    private static Neo4j embeddedNeo4j;
    private static Driver driver;
    private static GraphDatabaseService gds;

    @BeforeAll
    static void setUp() {
        embeddedNeo4j = Neo4jBuilders.newInProcessBuilder().build();
        driver = GraphDatabase.driver(embeddedNeo4j.boltURI());
        gds = embeddedNeo4j.defaultDatabaseService();

        try (Session session = driver.session()) {
            session.run("CREATE (:ARGRun {runId: 'r0', source: 'test'})");
            for (int n = 0; n < 5; n++) {
                double t = (n < 3) ? 0.0 : (n == 3 ? 0.5 : 1.0);
                boolean isSample = n < 3;
                long flags = isSample ? 1L : 0L;  // tskit NODE_IS_SAMPLE = 1
                session.run(
                    "CREATE (:TreeNode {treeNodeId: $tid, runId: 'r0', "
                  + "nodeId: $nid, time: $t, is_sample: $is, flags: $f})",
                    Map.of("tid", "r0:" + n, "nid", n,
                           "t", t, "is", isSample, "f", flags));
            }
            int[][] edges = {{4, 3}, {3, 0}, {3, 1}, {4, 2}};
            for (int[] e : edges) {
                session.run(
                    "MATCH (p:TreeNode {treeNodeId: $pid}), "
                  + "(c:TreeNode {treeNodeId: $cid}) "
                  + "CREATE (p)-[:PARENT_OF {runId: 'r0', start: 0, end: 100}]->(c)",
                    Map.of("pid", "r0:" + e[0], "cid", "r0:" + e[1]));
            }
        }
    }

    @AfterAll
    static void tearDown() {
        if (driver != null) driver.close();
        if (embeddedNeo4j != null) embeddedNeo4j.close();
    }

    @Test
    void load_recoversNodesAndEdges() {
        try (Transaction tx = gds.beginTx()) {
            ARG arg = ARGTraversal.load(tx, "r0", 0, 100);

            assertEquals(5, arg.nNodes);
            assertEquals(4, arg.nEdges);
            assertEquals(3, arg.nSamples());
            assertEquals("r0", arg.runId);

            // Nodes are returned ORDER BY nodeId, so packed index == tskit id.
            for (int i = 0; i < 5; i++) {
                assertEquals(i, arg.tskitNodeId[i]);
            }
            assertArrayEquals(
                new double[]{0.0, 0.0, 0.0, 0.5, 1.0}, arg.time, 1e-12);
            assertArrayEquals(
                new boolean[]{true, true, true, false, false}, arg.isSample);

            // Sample nodes in order.
            assertArrayEquals(new int[]{0, 1, 2}, arg.sampleNodes);
        }
    }

    @Test
    void load_edgeIntervalsAndPackedIds() {
        try (Transaction tx = gds.beginTx()) {
            ARG arg = ARGTraversal.load(tx, "r0", 0, 100);

            // Verify the four expected (parent_packed, child_packed) pairs are present.
            Set<String> expected = new HashSet<>(Arrays.asList(
                "4->3", "3->0", "3->1", "4->2"));
            Set<String> actual = new HashSet<>();
            for (int e = 0; e < arg.nEdges; e++) {
                actual.add(arg.edgeParent[e] + "->" + arg.edgeChild[e]);
                assertEquals(0L, arg.edgeStart[e]);
                assertEquals(100L, arg.edgeEnd[e]);
            }
            assertEquals(expected, actual);
        }
    }

    @Test
    void load_emptyRunReturnsEmptyArg() {
        try (Transaction tx = gds.beginTx()) {
            ARG arg = ARGTraversal.load(tx, "nonexistent_run", 0, 100);
            assertEquals(0, arg.nNodes);
            assertEquals(0, arg.nEdges);
            assertEquals(0, arg.nSamples());
        }
    }

    @Test
    void load_regionRestrictionDropsNonOverlappingEdges() {
        // The fixture's only run has edges spanning [0, 100]. A region
        // [200, 300] must drop them all.
        try (Transaction tx = gds.beginTx()) {
            ARG arg = ARGTraversal.load(tx, "r0", 200, 300);
            assertEquals(5, arg.nNodes, "all nodes still loaded");
            assertEquals(0, arg.nEdges, "no edges overlap [200,300]");
        }
    }

    @Test
    void breakpoints_areSortedAndUnique() {
        try (Transaction tx = gds.beginTx()) {
            ARG arg = ARGTraversal.load(tx, "r0", 0, 100);
            long[] bp = arg.breakpointsClipped(0L, 100L);
            // Only two distinct breakpoints: 0 and 100.
            assertArrayEquals(new long[]{0L, 100L}, bp);
        }
    }

    @Test
    void tskitIdToIndexMap_isInverse() {
        try (Transaction tx = gds.beginTx()) {
            ARG arg = ARGTraversal.load(tx, "r0", 0, 100);
            Map<Integer, Integer> m = arg.tskitIdToIndexMap();
            for (int i = 0; i < arg.nNodes; i++) {
                assertEquals((Integer) i, m.get(arg.tskitNodeId[i]));
            }
        }
    }
}
