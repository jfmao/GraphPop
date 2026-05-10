package org.graphpop.procedures.selection;

/**
 * Per-variant row from {@code graphpop.selection.allele_age_scan}.
 */
public class AlleleAgeScanResult {

    public String variant_id;
    public double freq;
    public long n_carriers;
    public long n_samples;
    public double age_midpoint;
    public double log_age;
    public double bin_mean_log_age;
    public double bin_sd_log_age;
    public long bin_index;
    public long bin_n;
    public double z_score;
    public String runId;

    public AlleleAgeScanResult() {}

    public AlleleAgeScanResult(String variantId, double freq,
                                long nCarriers, long nSamples,
                                double ageMid, double logAge,
                                double binMean, double binSd,
                                long binIndex, long binN,
                                double zScore, String runId) {
        this.variant_id = variantId;
        this.freq = freq;
        this.n_carriers = nCarriers;
        this.n_samples = nSamples;
        this.age_midpoint = ageMid;
        this.log_age = logAge;
        this.bin_mean_log_age = binMean;
        this.bin_sd_log_age = binSd;
        this.bin_index = binIndex;
        this.bin_n = binN;
        this.z_score = zScore;
        this.runId = runId;
    }
}
