package org.graphpop.procedures.pairwise;

/**
 * Result row from {@code graphpop.kinship.branch_grm_posterior}.
 *
 * <p>{@code b_ij_sd} is the standard error of the posterior mean
 * (Welford SE = SD/√n_runs); multiply by 1.96 for a Wald 95 % CI.
 * Returns {@code NaN} when the posterior set has only one run.</p>
 */
public class PosteriorKinshipResult {

    public String sample_a;
    public String sample_b;
    public double b_ij_mean;
    public double b_ij_sd;
    public long n_runs;
    public String method;

    public PosteriorKinshipResult() {}

    public PosteriorKinshipResult(String a, String b,
                                   double mean, double sd,
                                   long nRuns, String method) {
        this.sample_a = a;
        this.sample_b = b;
        this.b_ij_mean = mean;
        this.b_ij_sd = sd;
        this.n_runs = nRuns;
        this.method = method;
    }
}
