package org.graphpop.procedures.pairwise;

import org.neo4j.graphdb.Result;
import org.neo4j.graphdb.Transaction;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Cypher-driven loader for the {@link ARG} in-memory graph.
 *
 * <p>Loads every {@code :TreeNode} and {@code :PARENT_OF} edge for a given
 * {@code runId} into packed parallel arrays in one read transaction.
 * Mirrors the {@code HaplotypeMatrix} pattern (single Cypher pass, then
 * compute in-memory) used elsewhere in graphpop-procedures.</p>
 */
public final class ARGTraversal {

    private ARGTraversal() {}

    /**
     * Load every {@code :TreeNode} for the run, plus every {@code :PARENT_OF}
     * edge whose interval overlaps {@code [start, end]}.
     *
     * @param tx       active read transaction
     * @param runId    target run-id
     * @param start    region start (bp, inclusive)
     * @param end      region end (bp, exclusive); pass {@link Long#MAX_VALUE} for "no upper bound"
     * @return packed {@link ARG}; nodes are in insertion order from Cypher
     */
    public static ARG load(Transaction tx, String runId, long start, long end) {
        // 1. Nodes
        List<Integer> tskitIds = new ArrayList<>();
        List<Double> times = new ArrayList<>();
        List<Boolean> isSampleFlags = new ArrayList<>();
        List<Long> nodeFlags = new ArrayList<>();
        Map<Integer, Integer> tskitToPacked = new HashMap<>();

        Result r = tx.execute(
                "MATCH (n:TreeNode {runId: $runId}) "
              + "RETURN n.nodeId AS nodeId, n.time AS time, "
              + "n.is_sample AS is_sample, n.flags AS flags "
              + "ORDER BY n.nodeId",
                Map.of("runId", runId));
        try {
            while (r.hasNext()) {
                Map<String, Object> row = r.next();
                int packedIdx = tskitIds.size();
                int tskitId = ((Number) row.get("nodeId")).intValue();
                tskitIds.add(tskitId);
                times.add(((Number) row.get("time")).doubleValue());
                isSampleFlags.add((Boolean) row.get("is_sample"));
                nodeFlags.add(((Number) row.get("flags")).longValue());
                tskitToPacked.put(tskitId, packedIdx);
            }
        } finally {
            r.close();
        }

        int nNodes = tskitIds.size();
        if (nNodes == 0) {
            return new ARG(runId, 0, 0, new double[0], new boolean[0],
                    new long[0], new int[0],
                    new int[0], new int[0], new long[0], new long[0],
                    new int[0], 0L);
        }

        double[] timeArr = new double[nNodes];
        boolean[] isSampleArr = new boolean[nNodes];
        long[] flagsArr = new long[nNodes];
        int[] tskitArr = new int[nNodes];
        for (int i = 0; i < nNodes; i++) {
            timeArr[i] = times.get(i);
            isSampleArr[i] = isSampleFlags.get(i);
            flagsArr[i] = nodeFlags.get(i);
            tskitArr[i] = tskitIds.get(i);
        }

        // 2. Edges (filter by region overlap)
        List<int[]> edges = new ArrayList<>();
        List<long[]> intervals = new ArrayList<>();
        long maxEnd = 0L;

        Result er = tx.execute(
                "MATCH (p:TreeNode {runId: $runId})-[rel:PARENT_OF]->(c:TreeNode) "
              + "WHERE rel.runId = $runId "
              + "  AND rel.start < $end AND rel.end > $start "
              + "RETURN p.nodeId AS parent, c.nodeId AS child, "
              + "rel.start AS rel_start, rel.end AS rel_end",
                Map.of("runId", runId, "start", start, "end", end));
        try {
            while (er.hasNext()) {
                Map<String, Object> row = er.next();
                int parentTskit = ((Number) row.get("parent")).intValue();
                int childTskit = ((Number) row.get("child")).intValue();
                long s = ((Number) row.get("rel_start")).longValue();
                long t = ((Number) row.get("rel_end")).longValue();
                Integer parentPacked = tskitToPacked.get(parentTskit);
                Integer childPacked = tskitToPacked.get(childTskit);
                if (parentPacked == null || childPacked == null) continue;
                edges.add(new int[]{parentPacked, childPacked});
                intervals.add(new long[]{s, t});
                if (t > maxEnd) maxEnd = t;
            }
        } finally {
            er.close();
        }

        int nEdges = edges.size();
        int[] edgeParent = new int[nEdges];
        int[] edgeChild = new int[nEdges];
        long[] edgeStart = new long[nEdges];
        long[] edgeEnd = new long[nEdges];
        for (int e = 0; e < nEdges; e++) {
            edgeParent[e] = edges.get(e)[0];
            edgeChild[e] = edges.get(e)[1];
            edgeStart[e] = intervals.get(e)[0];
            edgeEnd[e] = intervals.get(e)[1];
        }

        // 3. Sample nodes (sample-flagged, in tskit order)
        List<Integer> samplesPacked = new ArrayList<>();
        for (int i = 0; i < nNodes; i++) {
            if (isSampleArr[i]) samplesPacked.add(i);
        }
        int[] sampleNodes = new int[samplesPacked.size()];
        for (int i = 0; i < samplesPacked.size(); i++) sampleNodes[i] = samplesPacked.get(i);

        long sequenceLength = (end == Long.MAX_VALUE) ? maxEnd : end;

        return new ARG(runId, nNodes, nEdges,
                timeArr, isSampleArr, flagsArr, tskitArr,
                edgeParent, edgeChild, edgeStart, edgeEnd,
                sampleNodes, sequenceLength);
    }
}
