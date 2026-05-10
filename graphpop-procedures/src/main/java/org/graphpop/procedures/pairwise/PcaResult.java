package org.graphpop.procedures.pairwise;

/**
 * Result row for {@code graphpop.kinship.branch_grm_pca}: one row
 * per (sample, principal-component) pair.
 */
public class PcaResult {

    public String sample_id;
    public long pc;            // 0-indexed
    public double value;
    public double eigenvalue;
    public String method;

    public PcaResult() {}

    public PcaResult(String sampleId, long pc, double value,
                     double eigenvalue, String method) {
        this.sample_id = sampleId;
        this.pc = pc;
        this.value = value;
        this.eigenvalue = eigenvalue;
        this.method = method;
    }
}
