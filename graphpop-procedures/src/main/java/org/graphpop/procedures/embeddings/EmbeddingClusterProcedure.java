package org.graphpop.procedures.embeddings;

import org.graphpop.procedures.embeddings.EmbeddingKnnProcedure.EmbeddingTable;
import org.neo4j.graphdb.Transaction;
import org.neo4j.procedure.Context;
import org.neo4j.procedure.Description;
import org.neo4j.procedure.Mode;
import org.neo4j.procedure.Name;
import org.neo4j.procedure.Procedure;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Stream;

/**
 * Per-sample clustering over {@code :Sample.embedding} vectors (M10).
 *
 * <pre>
 * CALL graphpop.embedding.cluster('kmeans', $k,
 *         {seed: 42, max_iter: 100})
 *   YIELD sample_id, cluster_id, distance_to_centroid,
 *         n_clusters, method
 * </pre>
 *
 * <p>v1 supports {@code method = 'kmeans'} (Lloyd's algorithm with
 * k-means++ initialisation; deterministic given a seed). HDBSCAN
 * is the natural follow-up; defer until needed.</p>
 */
public class EmbeddingClusterProcedure {

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.embedding.cluster", mode = Mode.READ)
    @Description("Cluster samples by their :Sample.embedding vectors. "
            + "v1: method='kmeans', Lloyd's algorithm with k-means++ "
            + "initialisation.")
    public Stream<EmbeddingClusterResult> cluster(
            @Name("method") String method,
            @Name("k") long k,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        if (options == null) options = new HashMap<>();
        if (!"kmeans".equalsIgnoreCase(method)) {
            throw new IllegalArgumentException(
                    "Only method='kmeans' is supported in v1; got: " + method);
        }
        if (k <= 0) {
            throw new IllegalArgumentException("k must be >= 1");
        }
        long seed = ((Number) options.getOrDefault("seed", 42L)).longValue();
        int maxIter = ((Number) options.getOrDefault("max_iter", 100L)).intValue();

        EmbeddingTable table = EmbeddingTable.load(tx);
        if (table.size() == 0) return Stream.empty();
        if (k > table.size()) {
            throw new IllegalArgumentException(
                    "k (" + k + ") exceeds the number of embedded samples ("
                  + table.size() + ")");
        }
        List<String> ids = table.sampleIds();
        List<double[]> vecs = table.vectors();

        int n = ids.size();
        int dim = vecs.get(0).length;
        int kInt = (int) k;

        // k-means++ initialisation.
        Random rng = new Random(seed);
        int[] centroidIdx = new int[kInt];
        centroidIdx[0] = rng.nextInt(n);
        double[] minDistSq = new double[n];
        for (int i = 0; i < n; i++) {
            minDistSq[i] = squaredEuclidean(vecs.get(i), vecs.get(centroidIdx[0]));
        }
        for (int c = 1; c < kInt; c++) {
            double total = 0.0;
            for (double d : minDistSq) total += d;
            if (total <= 0.0) {
                centroidIdx[c] = rng.nextInt(n);
            } else {
                double r = rng.nextDouble() * total;
                int pick = n - 1;
                double cum = 0.0;
                for (int i = 0; i < n; i++) {
                    cum += minDistSq[i];
                    if (cum >= r) { pick = i; break; }
                }
                centroidIdx[c] = pick;
            }
            double[] cv = vecs.get(centroidIdx[c]);
            for (int i = 0; i < n; i++) {
                double dsq = squaredEuclidean(vecs.get(i), cv);
                if (dsq < minDistSq[i]) minDistSq[i] = dsq;
            }
        }

        double[][] centroids = new double[kInt][dim];
        for (int c = 0; c < kInt; c++) {
            System.arraycopy(vecs.get(centroidIdx[c]), 0,
                              centroids[c], 0, dim);
        }

        int[] assign = new int[n];
        for (int iter = 0; iter < maxIter; iter++) {
            boolean changed = false;
            for (int i = 0; i < n; i++) {
                double[] x = vecs.get(i);
                int best = 0;
                double bestDsq = squaredEuclidean(x, centroids[0]);
                for (int c = 1; c < kInt; c++) {
                    double dsq = squaredEuclidean(x, centroids[c]);
                    if (dsq < bestDsq) {
                        bestDsq = dsq;
                        best = c;
                    }
                }
                if (best != assign[i]) {
                    assign[i] = best;
                    changed = true;
                }
            }
            // Recompute centroids.
            int[] counts = new int[kInt];
            double[][] sums = new double[kInt][dim];
            for (int i = 0; i < n; i++) {
                int c = assign[i];
                counts[c]++;
                double[] x = vecs.get(i);
                for (int d = 0; d < dim; d++) sums[c][d] += x[d];
            }
            for (int c = 0; c < kInt; c++) {
                if (counts[c] == 0) {
                    // Empty cluster: leave centroid in place. (Lloyd's
                    // can yield empty clusters; defer fancier handling.)
                    continue;
                }
                for (int d = 0; d < dim; d++) {
                    centroids[c][d] = sums[c][d] / counts[c];
                }
            }
            if (!changed) break;
        }

        List<EmbeddingClusterResult> rows = new ArrayList<>(n);
        for (int i = 0; i < n; i++) {
            double dist = Math.sqrt(squaredEuclidean(
                    vecs.get(i), centroids[assign[i]]));
            rows.add(new EmbeddingClusterResult(
                    ids.get(i), assign[i], dist, k, "kmeans"));
        }
        return rows.stream();
    }

    static double squaredEuclidean(double[] a, double[] b) {
        int n = Math.min(a.length, b.length);
        double s = 0.0;
        for (int i = 0; i < n; i++) {
            double d = a[i] - b[i];
            s += d * d;
        }
        return s;
    }
}
