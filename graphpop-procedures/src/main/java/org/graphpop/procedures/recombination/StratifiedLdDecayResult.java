package org.graphpop.procedures.recombination;

/**
 * Per-(window, stratum) row from
 * {@code graphpop.recombination.stratified_ld_decay}.
 */
public class StratifiedLdDecayResult {

    public long start;
    public long end;
    public String stratum;
    public long n_variant_pairs;
    public double mean_r2;
    public double mean_pair_distance;
    public double rho_per_bp;
    public long n_samples;
    public String stratify_by;
    public String runId;

    public StratifiedLdDecayResult() {}

    public StratifiedLdDecayResult(long start, long end, String stratum,
                                    long nVariantPairs, double meanR2,
                                    double meanPairDistance,
                                    double rhoPerBp, long nSamples,
                                    String stratifyBy, String runId) {
        this.start = start;
        this.end = end;
        this.stratum = stratum;
        this.n_variant_pairs = nVariantPairs;
        this.mean_r2 = meanR2;
        this.mean_pair_distance = meanPairDistance;
        this.rho_per_bp = rhoPerBp;
        this.n_samples = nSamples;
        this.stratify_by = stratifyBy;
        this.runId = runId;
    }
}
