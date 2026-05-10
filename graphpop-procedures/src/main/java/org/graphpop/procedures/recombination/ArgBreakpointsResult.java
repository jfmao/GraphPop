package org.graphpop.procedures.recombination;

/**
 * Per-window row from {@code graphpop.recombination.arg_breakpoints}.
 */
public class ArgBreakpointsResult {

    public long start;
    public long end;
    public long n_breakpoints;
    public long n_marginal_trees;
    public double total_branch_length;
    public double rho_per_bp;
    public String runId;

    public ArgBreakpointsResult() {}

    public ArgBreakpointsResult(long start, long end,
                                 long nBreakpoints, long nMarginalTrees,
                                 double totalBranchLength,
                                 double rhoPerBp, String runId) {
        this.start = start;
        this.end = end;
        this.n_breakpoints = nBreakpoints;
        this.n_marginal_trees = nMarginalTrees;
        this.total_branch_length = totalBranchLength;
        this.rho_per_bp = rhoPerBp;
        this.runId = runId;
    }
}
