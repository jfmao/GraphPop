package org.graphpop.procedures.recombination;

/**
 * Per-window row from {@code graphpop.recombination.ldhat_mcmc}.
 */
public class LdhatMcmcResult {

    public long start;
    public long end;
    public long n_variant_pairs;
    public double rho_posterior_mean;
    public double rho_lower_2_5;
    public double rho_upper_97_5;
    public long n_iter;
    public long n_accepted;
    public long n_samples;
    public String runId;

    public LdhatMcmcResult() {}

    public LdhatMcmcResult(long start, long end, long nVariantPairs,
                           double rhoMean, double rhoLower, double rhoUpper,
                           long nIter, long nAccepted, long nSamples,
                           String runId) {
        this.start = start;
        this.end = end;
        this.n_variant_pairs = nVariantPairs;
        this.rho_posterior_mean = rhoMean;
        this.rho_lower_2_5 = rhoLower;
        this.rho_upper_97_5 = rhoUpper;
        this.n_iter = nIter;
        this.n_accepted = nAccepted;
        this.n_samples = nSamples;
        this.runId = runId;
    }
}
