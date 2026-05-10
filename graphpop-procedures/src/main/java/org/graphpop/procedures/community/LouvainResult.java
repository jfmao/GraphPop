package org.graphpop.procedures.community;

/**
 * Per-sample row from {@code graphpop.community.louvain}.
 */
public class LouvainResult {

    public String sample_id;
    public long community_id;
    public double modularity;
    public long n_communities;
    public String source;

    public LouvainResult() {}

    public LouvainResult(String sampleId, long communityId,
                         double modularity, long nCommunities,
                         String source) {
        this.sample_id = sampleId;
        this.community_id = communityId;
        this.modularity = modularity;
        this.n_communities = nCommunities;
        this.source = source;
    }
}
