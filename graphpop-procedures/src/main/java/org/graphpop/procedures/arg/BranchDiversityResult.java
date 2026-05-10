package org.graphpop.procedures.arg;

/**
 * Result row for {@code graphpop.arg.branch_diversity}.
 */
public class BranchDiversityResult {

    public long start;
    public long end;
    public double branch_pi;
    public long n_samples;
    public String mode;
    public String runId;

    public BranchDiversityResult() {}

    public BranchDiversityResult(long start, long end, double branchPi,
                                  long nSamples, String mode, String runId) {
        this.start = start;
        this.end = end;
        this.branch_pi = branchPi;
        this.n_samples = nSamples;
        this.mode = mode;
        this.runId = runId;
    }
}
