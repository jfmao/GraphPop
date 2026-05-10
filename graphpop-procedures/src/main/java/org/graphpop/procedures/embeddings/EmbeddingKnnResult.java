package org.graphpop.procedures.embeddings;

/**
 * Per-neighbour row from {@code graphpop.embedding.knn}.
 */
public class EmbeddingKnnResult {

    public String query_sample_id;
    public String neighbor_sample_id;
    public double cosine_similarity;
    public long rank;

    public EmbeddingKnnResult() {}

    public EmbeddingKnnResult(String querySampleId, String neighborSampleId,
                               double cosineSimilarity, long rank) {
        this.query_sample_id = querySampleId;
        this.neighbor_sample_id = neighborSampleId;
        this.cosine_similarity = cosineSimilarity;
        this.rank = rank;
    }
}
