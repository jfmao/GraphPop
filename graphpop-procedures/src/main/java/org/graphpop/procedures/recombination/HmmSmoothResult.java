package org.graphpop.procedures.recombination;

/**
 * Per-window row from {@code graphpop.recombination.hmm_smooth}.
 */
public class HmmSmoothResult {

    public long start;
    public long end;
    public double rho_per_bp;
    public double rho_smoothed;
    public long hmm_state;
    public long n_variant_pairs;
    public String runId;

    public HmmSmoothResult() {}

    public HmmSmoothResult(long start, long end,
                            double rhoPerBp, double rhoSmoothed,
                            long hmmState, long nVariantPairs,
                            String runId) {
        this.start = start;
        this.end = end;
        this.rho_per_bp = rhoPerBp;
        this.rho_smoothed = rhoSmoothed;
        this.hmm_state = hmmState;
        this.n_variant_pairs = nVariantPairs;
        this.runId = runId;
    }
}
