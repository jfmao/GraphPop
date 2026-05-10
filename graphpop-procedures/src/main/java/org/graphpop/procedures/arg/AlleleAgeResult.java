package org.graphpop.procedures.arg;

/**
 * Result row for {@code graphpop.arg.allele_age}.
 */
public class AlleleAgeResult {

    public String variant_id;
    public long child_node_id;
    public long parent_node_id;
    public double child_time;
    public double parent_time;
    public double midpoint_time;
    public long n_carriers;
    public String runId;

    public AlleleAgeResult() {}

    public AlleleAgeResult(String variantId,
                            long childNodeId, long parentNodeId,
                            double childTime, double parentTime,
                            double midpointTime, long nCarriers,
                            String runId) {
        this.variant_id = variantId;
        this.child_node_id = childNodeId;
        this.parent_node_id = parentNodeId;
        this.child_time = childTime;
        this.parent_time = parentTime;
        this.midpoint_time = midpointTime;
        this.n_carriers = nCarriers;
        this.runId = runId;
    }
}
