package org.graphpop.procedures.arg;


import org.graphpop.procedures.pairwise.ARG;
import org.graphpop.procedures.pairwise.ARGTraversal;
import org.neo4j.graphdb.Result;
import org.neo4j.graphdb.Transaction;
import org.neo4j.procedure.Context;
import org.neo4j.procedure.Description;
import org.neo4j.procedure.Mode;
import org.neo4j.procedure.Name;
import org.neo4j.procedure.Procedure;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

/**
 * Allele age from an ARG (M6 part D).
 *
 * <pre>
 * CALL graphpop.arg.allele_age($runId, $variantId)
 *   YIELD variant_id, child_node_id, parent_node_id,
 *         child_time, parent_time, midpoint_time,
 *         n_carriers, runId
 * </pre>
 *
 * <p>Reads the {@code :MUTATED_ON {runId, parent_node_id}} edge for
 * the variant, looks up both endpoint times, and returns the bracket
 * plus a midpoint and the carrier count (descendants of the child
 * TreeNode in the marginal tree containing the mutation site).</p>
 */
public class ArgAlleleAgeProcedure {

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.arg.allele_age", mode = Mode.READ)
    @Description("Allele age from an ARG: returns the time bracket of the "
            + ":MUTATED_ON edge for a variant plus midpoint + carrier count.")
    public Stream<AlleleAgeResult> alleleAge(
            @Name("runId") String runId,
            @Name("variantId") String variantId
    ) {
        Result r = tx.execute(
                "MATCH (v:Variant {variantId: $vid})"
              + "-[m:MUTATED_ON {runId: $rid}]->(c:TreeNode) "
              + "RETURN c.nodeId AS child_id, c.time AS child_time, "
              + "m.parent_node_id AS parent_id, v.position AS position",
                Map.of("vid", variantId, "rid", runId));
        if (!r.hasNext()) {
            r.close();
            return Stream.empty();
        }

        Map<String, Object> row = r.next();
        r.close();

        long childTskit = ((Number) row.get("child_id")).longValue();
        double childTime = ((Number) row.get("child_time")).doubleValue();
        long parentTskit = ((Number) row.get("parent_id")).longValue();
        long position = ((Number) row.get("position")).longValue();

        Result pr = tx.execute(
                "MATCH (n:TreeNode {runId: $rid, nodeId: $pid}) "
              + "RETURN n.time AS t",
                Map.of("rid", runId, "pid", parentTskit));
        if (!pr.hasNext()) {
            pr.close();
            return Stream.empty();
        }
        double parentTime = ((Number) pr.next().get("t")).doubleValue();
        pr.close();

        // Carrier count: descendants of `child` in the marginal tree at
        // `position`. Reuse the in-memory ARG to make the descendant DFS
        // cheap.
        long nCarriers = 0L;
        ARG arg = ARGTraversal.load(tx, runId, 0L, Long.MAX_VALUE);
        if (arg.nNodes > 0) {
            Map<Integer, Integer> idMap = arg.tskitIdToIndexMap();
            Integer childPacked = idMap.get((int) childTskit);
            if (childPacked != null) {
                int[] parent = ArgTraversalUtils.marginalTree(arg, position);
                nCarriers = countSampleDescendants(arg, parent, childPacked);
            }
        }

        AlleleAgeResult res = new AlleleAgeResult(
                variantId, childTskit, parentTskit,
                childTime, parentTime, 0.5 * (childTime + parentTime),
                nCarriers, runId);
        return List.of(res).stream();
    }

    /**
     * Count sample-flagged descendants of {@code root} in the marginal
     * tree defined by {@code parent[]}.
     */
    private static long countSampleDescendants(ARG arg, int[] parent, int root) {
        long count = 0L;
        // Walk from every sample upward, count those whose path passes
        // through root.
        for (int leaf : arg.sampleNodes) {
            int cur = leaf;
            while (cur != -1) {
                if (cur == root) { count++; break; }
                cur = parent[cur];
            }
        }
        return count;
    }
}
