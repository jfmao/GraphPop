package org.graphpop.procedures.recombination;

/**
 * Per-window row from {@code graphpop.recombination.ld_decay}.
 */
public class LdDecayResult {

    public long start;
    public long end;
    public long n_variant_pairs;
    public double mean_r2;
    public double mean_pair_distance;
    public double rho_per_bp;
    public long n_samples;
    public String method;
    public String runId;

    public LdDecayResult() {}

    public LdDecayResult(long start, long end,
                         long nVariantPairs, double meanR2,
                         double meanPairDistance, double rhoPerBp,
                         long nSamples, String method, String runId) {
        this.start = start;
        this.end = end;
        this.n_variant_pairs = nVariantPairs;
        this.mean_r2 = meanR2;
        this.mean_pair_distance = meanPairDistance;
        this.rho_per_bp = rhoPerBp;
        this.n_samples = nSamples;
        this.method = method;
        this.runId = runId;
    }
}
