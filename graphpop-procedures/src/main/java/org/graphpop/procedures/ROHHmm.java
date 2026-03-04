package org.graphpop.procedures;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * Hidden Markov Model for Runs of Homozygosity detection.
 *
 * <p>Implements the 2-state Viterbi HMM approach of Narasimhan et al. (2016,
 * <a href="https://doi.org/10.1093/bioinformatics/btw044">bcftools roh</a>).
 * Two hidden states:</p>
 * <ul>
 *   <li><b>AZ</b> (autozygous) — both chromosomes identical-by-descent</li>
 *   <li><b>HW</b> (Hardy–Weinberg) — normal outbred state</li>
 * </ul>
 *
 * <p>Emission probabilities are weighted by population allele frequency:
 * a rare-allele homozygote is strong evidence for AZ (P_AZ = p vs P_HW = p²).
 * Transitions are scaled by physical distance between variants.</p>
 *
 * <p>All computation in log-space to avoid underflow over long chromosomes.</p>
 */
final class ROHHmm {

    private static final int STATE_HW = 0;
    private static final int STATE_AZ = 1;

    /** Expected ROH length in bp. Controls transition AZ→HW rate. Default: 500 kb. */
    private final double azLength;

    /** Expected fraction of genome in ROH. Controls transition HW→AZ rate. Default: 1%. */
    private final double azFraction;

    /** P(het | AZ) — genotyping error rate. Default: 1e-5. */
    private final double hetErrorRate;

    /**
     * Minimum population AF to include a variant. Monomorphic and very rare
     * variants provide near-zero information for ROH detection. bcftools roh
     * defaults to AF > 0 and suggests AF ≥ 0.05 for WGS data. Default: 0.0
     * (skip only AF=0 monomorphic sites).
     */
    private final double minAf;

    /** Minimum AF/max AF clamp to avoid log(0). */
    private static final double AF_FLOOR = 0.001;
    private static final double AF_CEIL = 0.999;

    ROHHmm(double azLength, double azFraction, double hetErrorRate, double minAf) {
        this.azLength = azLength;
        this.azFraction = azFraction;
        this.hetErrorRate = hetErrorRate;
        this.minAf = minAf;
    }

    /**
     * Create from procedure options map.
     */
    static ROHHmm fromOptions(Map<String, Object> options) {
        double azLen = getDouble(options, "az_length", 500_000.0);
        double azFrac = getDouble(options, "az_fraction", 0.01);
        double hetErr = getDouble(options, "het_error_rate", 1e-5);
        double maf = getDouble(options, "hmm_min_af", 0.0);
        return new ROHHmm(azLen, azFrac, hetErr, maf);
    }

    /**
     * Run Viterbi decoding on one sample.
     *
     * @param positions variant positions (sorted, length N)
     * @param afs       population allele frequencies (length N)
     * @param genotypes 0=hom-ref, 1=het, 2=hom-alt (length N)
     * @param minLength minimum ROH length to report (bp)
     * @param minSnps   minimum SNPs in ROH to report
     * @return list of ROH segments as [startPos, endPos] arrays
     */
    List<long[]> viterbi(long[] positions, double[] afs, int[] genotypes,
                         long minLength, int minSnps) {
        int rawN = positions.length;
        if (rawN == 0) return List.of();

        // Pre-filter to informative variants (AF >= minAf and AF <= 1-minAf).
        // Monomorphic variants provide essentially zero signal and dilute the
        // transition model by making inter-variant distances artificially small.
        // Use folded AF: min(af, 1-af) to catch AF near both 0 and 1.
        int[] idx;
        if (minAf > 0) {
            int count = 0;
            for (int i = 0; i < rawN; i++) {
                double faf = Math.min(afs[i], 1.0 - afs[i]);
                if (faf >= minAf) count++;
            }
            idx = new int[count];
            int k = 0;
            for (int i = 0; i < rawN; i++) {
                double faf = Math.min(afs[i], 1.0 - afs[i]);
                if (faf >= minAf) idx[k++] = i;
            }
        } else {
            // Skip only truly monomorphic (AF == 0 or AF == 1)
            int count = 0;
            for (int i = 0; i < rawN; i++) {
                if (afs[i] > 0 && afs[i] < 1.0) count++;
            }
            idx = new int[count];
            int k = 0;
            for (int i = 0; i < rawN; i++) {
                if (afs[i] > 0 && afs[i] < 1.0) idx[k++] = i;
            }
        }

        int n = idx.length;
        if (n == 0) return List.of();

        // Initial state probabilities (log-space)
        double logPiHW = Math.log(1.0 - azFraction);
        double logPiAZ = Math.log(azFraction);

        // Viterbi arrays
        double[] prevLogProb = new double[2]; // [HW, AZ]
        int[][] backpointer = new int[n][2];

        // Emission at first variant
        double p = clampAF(afs[idx[0]]);
        prevLogProb[STATE_HW] = logPiHW + logEmissionHW(genotypes[idx[0]], p);
        prevLogProb[STATE_AZ] = logPiAZ + logEmissionAZ(genotypes[idx[0]], p);

        // Forward pass
        for (int i = 1; i < n; i++) {
            double d = positions[idx[i]] - positions[idx[i - 1]];
            p = clampAF(afs[idx[i]]);

            // Transition probabilities (log-space)
            double pAZtoHW = 1.0 - Math.exp(-d / azLength);
            double pHWtoAZ = 1.0 - Math.exp(-d / ((1.0 / azFraction - 1.0) * azLength));

            // Clamp transition probabilities away from 0 and 1
            pAZtoHW = clampProb(pAZtoHW);
            pHWtoAZ = clampProb(pHWtoAZ);

            double logAZtoHW = Math.log(pAZtoHW);
            double logAZtoAZ = Math.log(1.0 - pAZtoHW);
            double logHWtoAZ = Math.log(pHWtoAZ);
            double logHWtoHW = Math.log(1.0 - pHWtoAZ);

            double logEmHW = logEmissionHW(genotypes[idx[i]], p);
            double logEmAZ = logEmissionAZ(genotypes[idx[i]], p);

            double[] currLogProb = new double[2];

            // Best path to HW state
            double fromHW = prevLogProb[STATE_HW] + logHWtoHW;
            double fromAZ = prevLogProb[STATE_AZ] + logAZtoHW;
            if (fromHW >= fromAZ) {
                currLogProb[STATE_HW] = fromHW + logEmHW;
                backpointer[i][STATE_HW] = STATE_HW;
            } else {
                currLogProb[STATE_HW] = fromAZ + logEmHW;
                backpointer[i][STATE_HW] = STATE_AZ;
            }

            // Best path to AZ state
            fromHW = prevLogProb[STATE_HW] + logHWtoAZ;
            fromAZ = prevLogProb[STATE_AZ] + logAZtoAZ;
            if (fromAZ >= fromHW) {
                currLogProb[STATE_AZ] = fromAZ + logEmAZ;
                backpointer[i][STATE_AZ] = STATE_AZ;
            } else {
                currLogProb[STATE_AZ] = fromHW + logEmAZ;
                backpointer[i][STATE_AZ] = STATE_HW;
            }

            prevLogProb = currLogProb;
        }

        // Traceback
        int[] path = new int[n];
        path[n - 1] = prevLogProb[STATE_AZ] >= prevLogProb[STATE_HW] ? STATE_AZ : STATE_HW;
        for (int i = n - 2; i >= 0; i--) {
            path[i] = backpointer[i + 1][path[i + 1]];
        }

        // Merge consecutive AZ states into segments (using idx mapping)
        List<long[]> segments = new ArrayList<>();
        int segStartI = -1;
        for (int i = 0; i < n; i++) {
            if (path[i] == STATE_AZ) {
                if (segStartI < 0) segStartI = i;
            } else {
                if (segStartI >= 0) {
                    addSegmentIfPassing(segments, positions, idx,
                            segStartI, i - 1, minLength, minSnps);
                    segStartI = -1;
                }
            }
        }
        if (segStartI >= 0) {
            addSegmentIfPassing(segments, positions, idx,
                    segStartI, n - 1, minLength, minSnps);
        }

        return segments;
    }

    private void addSegmentIfPassing(List<long[]> segments, long[] positions,
                                     int[] idx, int startI, int endI,
                                     long minLength, int minSnps) {
        int nSnps = endI - startI + 1;
        long startPos = positions[idx[startI]];
        long endPos = positions[idx[endI]];
        long length = endPos - startPos;
        if (length >= minLength && nSnps >= minSnps) {
            segments.add(new long[]{startPos, endPos});
        }
    }

    /**
     * Log emission probability for HW state.
     * P(hom-ref | HW) = (1-p)^2
     * P(het     | HW) = 2*p*(1-p)
     * P(hom-alt | HW) = p^2
     */
    private double logEmissionHW(int genotype, double p) {
        switch (genotype) {
            case 0: return 2.0 * Math.log(1.0 - p);
            case 1: return Math.log(2.0 * p * (1.0 - p));
            case 2: return 2.0 * Math.log(p);
            default: return 0.0; // missing data — neutral
        }
    }

    /**
     * Log emission probability for AZ state.
     * P(hom-ref | AZ) = (1-p)
     * P(het     | AZ) = hetErrorRate
     * P(hom-alt | AZ) = p
     */
    private double logEmissionAZ(int genotype, double p) {
        switch (genotype) {
            case 0: return Math.log(1.0 - p);
            case 1: return Math.log(hetErrorRate);
            case 2: return Math.log(p);
            default: return 0.0;
        }
    }

    private static double clampAF(double af) {
        return Math.max(AF_FLOOR, Math.min(AF_CEIL, af));
    }

    private static double clampProb(double p) {
        return Math.max(1e-10, Math.min(1.0 - 1e-10, p));
    }

    private static double getDouble(Map<String, Object> options, String key, double defaultValue) {
        if (options == null) return defaultValue;
        Object val = options.get(key);
        if (val == null) return defaultValue;
        return ((Number) val).doubleValue();
    }
}
