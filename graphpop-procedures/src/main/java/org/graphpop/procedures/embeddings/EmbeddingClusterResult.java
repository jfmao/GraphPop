package org.graphpop.procedures.embeddings;

/**
 * Per-sample row from {@code graphpop.embedding.cluster}.
 */
public class EmbeddingClusterResult {

    public String sample_id;
    public long cluster_id;
    public double distance_to_centroid;
    public long n_clusters;
    public String method;

    public EmbeddingClusterResult() {}

    public EmbeddingClusterResult(String sampleId, long clusterId,
                                   double distanceToCentroid,
                                   long nClusters, String method) {
        this.sample_id = sampleId;
        this.cluster_id = clusterId;
        this.distance_to_centroid = distanceToCentroid;
        this.n_clusters = nClusters;
        this.method = method;
    }
}
