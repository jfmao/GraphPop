package org.graphpop.procedures.pairwise;

/**
 * Result row from {@code graphpop.kinship.branch_grm_apply}: one row
 * per (sample, vector-column) pair carrying the value of {@code G · v}
 * at that sample's index.
 */
public class KinshipVectorResult {

    public String sample_id;
    public long col;
    public double value;
    public long n_branches;
    public String method;

    public KinshipVectorResult() {}

    public KinshipVectorResult(String sampleId, long col, double value,
                                long nBranches, String method) {
        this.sample_id = sampleId;
        this.col = col;
        this.value = value;
        this.n_branches = nBranches;
        this.method = method;
    }
}
