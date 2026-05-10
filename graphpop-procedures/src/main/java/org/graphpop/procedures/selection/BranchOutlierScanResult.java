package org.graphpop.procedures.selection;

/**
 * Per-window row from {@code graphpop.selection.branch_outlier_scan}.
 */
public class BranchOutlierScanResult {

    public long start;
    public long end;
    public double total_branch_length;
    public double mean_total;
    public double sd_total;
    public double z_score;
    public long n_samples;
    public String runId;

    public BranchOutlierScanResult() {}

    public BranchOutlierScanResult(long start, long end,
                                    double total, double meanTotal,
                                    double sdTotal, double z,
                                    long nSamples, String runId) {
        this.start = start;
        this.end = end;
        this.total_branch_length = total;
        this.mean_total = meanTotal;
        this.sd_total = sdTotal;
        this.z_score = z;
        this.n_samples = nSamples;
        this.runId = runId;
    }
}
