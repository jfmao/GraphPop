package org.graphpop.procedures.pairwise;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

/**
 * Ancestry-decomposed branch GRM (Tang &amp; Chiang 2025
 * <i>Genetics</i>, "On ARGs, pedigrees, and genetic relatedness
 * matrices"). Buckets each branch's eGRM contribution by the
 * ancestry painting on the branch's child node, normalises every
 * per-ancestry sub-matrix by the <em>shared</em> unconditional
 * {@code total_mu}, and double-centers each independently.
 *
 * <p>For deterministic single-label painting (one
 * {@code :HAS_ANCESTRY} edge per node with
 * {@code posterior_prob = 1.0}), this is the natural partition
 * {@code Σ_a B_ij^a == B_ij^uncond} (modulo unpainted nodes). For
 * probabilistic painting (multiple edges per node with probabilities
 * summing to ≤ 1), the contribution is split fractionally.</p>
 *
 * <p>Same {@code n &le; 64} sample cap as {@link BranchGrmComputer}.</p>
 */
public final class BranchGrmByAncestryComputer {

    private BranchGrmByAncestryComputer() {}

    /**
     * Per-ancestry weighted edge: ancestry label and weight in [0, 1].
     */
    public static final class AncestryAssignment {
        public final String ancestry;
        public final double prob;

        public AncestryAssignment(String ancestry, double prob) {
            this.ancestry = ancestry;
            this.prob = prob;
        }
    }

    public static final class Result {
        public final int n;
        /** Insertion-ordered ancestry-label set (deterministic output ordering). */
        public final List<String> ancestries;
        public final Map<String, double[][]> matrices;
        public final double totalMu;

        Result(int n, List<String> ancestries,
                Map<String, double[][]> matrices, double totalMu) {
            this.n = n;
            this.ancestries = ancestries;
            this.matrices = matrices;
            this.totalMu = totalMu;
        }

        public double bij(String ancestry, int i, int j) {
            return matrices.get(ancestry)[i][j];
        }
    }

    /**
     * Compute the ancestry-decomposed eGRM.
     *
     * @param arg              the loaded ARG
     * @param regionStart      bp start (inclusive)
     * @param regionEnd        bp end (exclusive); {@link Long#MAX_VALUE} for full
     * @param weightFn         step-4 conditional weight; pass
     *                         {@link BranchWeightFn#UNIT} for unconditional
     * @param paintingByTskitId map from tskit node id -&gt; list of
     *                          (ancestry, prob) assignments. Nodes
     *                          missing from the map are treated as
     *                          unpainted (skipped from any ancestry's
     *                          numerator but still counted in
     *                          {@code total_mu}).
     */
    public static Result compute(ARG arg, long regionStart, long regionEnd,
                                  BranchWeightFn weightFn,
                                  Map<Integer, List<AncestryAssignment>> paintingByTskitId) {
        final int n = arg.nSamples();
        if (n > 64) {
            throw new IllegalArgumentException(
                "BranchGrmByAncestryComputer (v1) requires n_samples <= 64; got " + n);
        }
        if (n == 0) {
            return new Result(0, List.of(), Map.of(), 0.0);
        }

        // Collect ancestry label set in deterministic (sorted) order.
        TreeSet<String> labelSet = new TreeSet<>();
        for (List<AncestryAssignment> list : paintingByTskitId.values()) {
            for (AncestryAssignment a : list) labelSet.add(a.ancestry);
        }
        if (labelSet.isEmpty()) {
            // Empty painting -> empty result. The procedure layer will
            // skip emission entirely.
            return new Result(n, List.of(), Map.of(), 0.0);
        }
        List<String> ancestries = List.copyOf(labelSet);
        Map<String, double[][]> mats = new LinkedHashMap<>();
        for (String a : ancestries) mats.put(a, new double[n][n]);

        final int nNodes = arg.nNodes;
        final int[] packedToSampleBit = new int[nNodes];
        Arrays.fill(packedToSampleBit, -1);
        for (int k = 0; k < n; k++) {
            packedToSampleBit[arg.sampleNodes[k]] = k;
        }

        final long lo = Math.max(regionStart, 0L);
        final long hi = (regionEnd == Long.MAX_VALUE) ? arg.sequenceLength
                : Math.min(regionEnd, arg.sequenceLength);
        if (hi <= lo) {
            return new Result(n, ancestries, mats, 0.0);
        }

        final long[] breakpoints = arg.breakpointsClipped(lo, hi);
        if (breakpoints.length < 2) {
            return new Result(n, ancestries, mats, 0.0);
        }

        double totalMu = 0.0;
        final int[] parent = new int[nNodes];
        final long[] descMask = new long[nNodes];

        for (int k = 0; k < breakpoints.length - 1; k++) {
            final long b0 = breakpoints[k];
            final long b1 = breakpoints[k + 1];
            final long intervalLen = b1 - b0;
            if (intervalLen <= 0) continue;

            Arrays.fill(parent, -1);
            for (int e = 0; e < arg.nEdges; e++) {
                if (arg.edgeStart[e] <= b0 && arg.edgeEnd[e] > b0) {
                    parent[arg.edgeChild[e]] = arg.edgeParent[e];
                }
            }

            Arrays.fill(descMask, 0L);
            for (int sBit = 0; sBit < n; sBit++) {
                int cur = arg.sampleNodes[sBit];
                long bit = 1L << sBit;
                while (cur != -1) {
                    descMask[cur] |= bit;
                    cur = parent[cur];
                }
            }

            for (int c = 0; c < nNodes; c++) {
                final int p = parent[c];
                if (p == -1) continue;
                final long mask = descMask[c];
                final int nDesc = Long.bitCount(mask);
                if (nDesc == 0 || nDesc == n) continue;

                final double branchLen = arg.time[p] - arg.time[c];
                if (branchLen <= 0.0) continue;

                final double w = weightFn.weight(
                        arg.tskitNodeId[p], arg.tskitNodeId[c],
                        arg.time[p], arg.time[c], b0, b1);
                if (w <= 0.0) continue;

                final double mu = intervalLen * branchLen * w * 1e-8;
                totalMu += mu;

                final double pFreq = (double) nDesc / (double) n;
                final double base = mu / (pFreq * (1.0 - pFreq));

                int childTskitId = arg.tskitNodeId[c];
                List<AncestryAssignment> assignments =
                        paintingByTskitId.getOrDefault(
                                childTskitId, Collections.emptyList());
                if (assignments.isEmpty()) continue;

                long m = mask;
                int[] desc = new int[nDesc];
                int di = 0;
                while (m != 0L) {
                    desc[di++] = Long.numberOfTrailingZeros(m);
                    m &= m - 1L;
                }

                for (AncestryAssignment a : assignments) {
                    if (a.prob <= 0.0) continue;
                    double weight = base * a.prob;
                    double[][] mat = mats.get(a.ancestry);
                    for (int i = 0; i < nDesc; i++) {
                        double[] r = mat[desc[i]];
                        for (int j = 0; j < nDesc; j++) {
                            r[desc[j]] += weight;
                        }
                    }
                }
            }
        }

        if (totalMu == 0.0) {
            return new Result(n, ancestries, mats, 0.0);
        }

        // Normalize by shared total_mu, then double-center each ancestry.
        for (String a : ancestries) {
            double[][] m = mats.get(a);
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    m[i][j] /= totalMu;
            // Subtract column means.
            double[] colMeans = new double[n];
            for (int j = 0; j < n; j++) {
                double s = 0.0;
                for (int i = 0; i < n; i++) s += m[i][j];
                colMeans[j] = s / n;
            }
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    m[i][j] -= colMeans[j];
            // Subtract row means of the result.
            double[] rowMeans = new double[n];
            for (int i = 0; i < n; i++) {
                double s = 0.0;
                for (int j = 0; j < n; j++) s += m[i][j];
                rowMeans[i] = s / n;
            }
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    m[i][j] -= rowMeans[i];
        }

        return new Result(n, ancestries, mats, totalMu);
    }

    /**
     * Load the {@code :HAS_ANCESTRY} painting for a run into a map
     * keyed on tskit node id. Optional {@code painter} filter when a
     * run has multiple paintings.
     */
    public static Map<Integer, List<AncestryAssignment>> loadPainting(
            org.neo4j.graphdb.Transaction tx, String runId, String painter) {
        Map<Integer, List<AncestryAssignment>> out = new HashMap<>();
        String cypher = (painter == null)
                ? "MATCH (n:TreeNode {runId: $runId})-[r:HAS_ANCESTRY {runId: $runId}]->(p:Population) "
                + "RETURN n.nodeId AS nodeId, p.populationId AS pop, r.posterior_prob AS prob"
                : "MATCH (n:TreeNode {runId: $runId})"
                + "-[r:HAS_ANCESTRY {runId: $runId, painter: $painter}]->(p:Population) "
                + "RETURN n.nodeId AS nodeId, p.populationId AS pop, r.posterior_prob AS prob";
        Map<String, Object> params = (painter == null)
                ? Map.of("runId", runId)
                : Map.of("runId", runId, "painter", painter);
        org.neo4j.graphdb.Result r = tx.execute(cypher, params);
        try {
            while (r.hasNext()) {
                Map<String, Object> row = r.next();
                int nodeId = ((Number) row.get("nodeId")).intValue();
                String pop = (String) row.get("pop");
                double prob = ((Number) row.get("prob")).doubleValue();
                out.computeIfAbsent(nodeId, k -> new java.util.ArrayList<>())
                        .add(new AncestryAssignment(pop, prob));
            }
        } finally {
            r.close();
        }
        return out;
    }
}
