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
 *   <li><b>HW</b> (Hardy-Weinberg) — normal outbred state</li>
 * </ul>
 *
 * <p>Transition probabilities use the exact 2-state continuous-time Markov chain
 * (CTMC) formula, matching bcftools roh's matrix power approach. For two
 * consecutive variants separated by physical distance d (bp):</p>
 * <pre>
 *   decay = (1 - hwToAz - azToHw)^d
 *   P(HW-&gt;AZ | d) = hwToAz/(hwToAz+azToHw) * (1 - decay)
 *   P(AZ-&gt;HW | d) = azToHw/(hwToAz+azToHw) * (1 - decay)
 * </pre>
 *
 * <p>Emission probabilities are weighted by population allele frequency:
 * a rare-allele homozygote is strong evidence for AZ (P_AZ = p vs P_HW = p^2).
 * All computation in log-space to avoid underflow over long chromosomes.</p>
 */
final class ROHHmm {

    private static final int STATE_HW = 0;
    private static final int STATE_AZ = 1;

    /** Per-bp transition rate from HW to AZ state. bcftools default: 6.7e-8. */
    private final double hwToAz;

    /** Per-bp transition rate from AZ to HW state. bcftools default: 5e-9. */
    private final double azToHw;

    /**
     * P(het | AZ) — genotyping error rate. Default: 0.001 (matches bcftools
     * -G30, where PL=30 → error probability = 10^(-30/10) = 0.001).
     * Lower values (e.g. 1e-5) make ROH calls more conservative (closer to
     * PLINK); higher values are more liberal (closer to bcftools).
     */
    private final double hetErrorRate;

    /**
     * Minimum population AF to include a variant. Monomorphic and very rare
     * variants provide near-zero information for ROH detection. Default: 0.0
     * (skip only AF=0 monomorphic sites).
     */
    private final double minAf;

    /** Minimum AF/max AF clamp to avoid log(0). */
    private static final double AF_FLOOR = 0.001;
    private static final double AF_CEIL = 0.999;

    /** bcftools default transition rates (per bp). */
    static final double DEFAULT_HW_TO_AZ = 6.7e-8;
    static final double DEFAULT_AZ_TO_HW = 5.0e-9;

    ROHHmm(double hwToAz, double azToHw, double hetErrorRate, double minAf) {
        this.hwToAz = hwToAz;
        this.azToHw = azToHw;
        this.hetErrorRate = hetErrorRate;
        this.minAf = minAf;
    }

    /**
     * Create from procedure options map.
     *
     * <p>Accepts two parameterization styles:</p>
     * <ul>
     *   <li><b>Per-bp rates</b> (preferred, matches bcftools):
     *       {@code hw_to_az} and {@code az_to_hw}</li>
     *   <li><b>Length/fraction</b> (legacy):
     *       {@code az_length} (expected ROH length in bp) and
     *       {@code az_fraction} (expected genome fraction in ROH).
     *       Converted to per-bp rates: azToHw = 1/az_length,
     *       hwToAz = az_fraction / ((1-az_fraction) * az_length).</li>
     * </ul>
     *
     * <p>If neither is specified, bcftools defaults are used:
     * hw_to_az=6.7e-8, az_to_hw=5e-9.</p>
     */
    static ROHHmm fromOptions(Map<String, Object> options) {
        double hetErr = getDouble(options, "het_error_rate", 1e-3);
        double maf = getDouble(options, "hmm_min_af", 0.0);

        // Check for per-bp rate parameterization first
        Double hwToAzOpt = getDoubleOrNull(options, "hw_to_az");
        Double azToHwOpt = getDoubleOrNull(options, "az_to_hw");

        if (hwToAzOpt != null || azToHwOpt != null) {
            // Per-bp rates specified directly
            double hwAz = hwToAzOpt != null ? hwToAzOpt : DEFAULT_HW_TO_AZ;
            double azHw = azToHwOpt != null ? azToHwOpt : DEFAULT_AZ_TO_HW;
            return new ROHHmm(hwAz, azHw, hetErr, maf);
        }

        // Check for legacy az_length/az_fraction parameterization
        Double azLenOpt = getDoubleOrNull(options, "az_length");
        Double azFracOpt = getDoubleOrNull(options, "az_fraction");

        if (azLenOpt != null || azFracOpt != null) {
            double azLen = azLenOpt != null ? azLenOpt : 500_000.0;
            double azFrac = azFracOpt != null ? azFracOpt : 0.01;
            // Convert to per-bp rates
            double azHw = 1.0 / azLen;
            double hwAz = azFrac / ((1.0 - azFrac) * azLen);
            return new ROHHmm(hwAz, azHw, hetErr, maf);
        }

        // Default: bcftools rates
        return new ROHHmm(DEFAULT_HW_TO_AZ, DEFAULT_AZ_TO_HW, hetErr, maf);
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

        // Initial state probabilities from stationary distribution (log-space)
        double sumRate = hwToAz + azToHw;
        double piAZ = hwToAz / sumRate;  // stationary probability of AZ
        double piHW = azToHw / sumRate;   // stationary probability of HW
        double logPiHW = Math.log(piHW);
        double logPiAZ = Math.log(piAZ);

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

            // Transition probabilities — exact 2-state CTMC (matrix power T^d).
            // For per-bp base matrix T = [[1-hwToAz, hwToAz], [azToHw, 1-azToHw]]:
            //   eigenvalues: 1 and (1-hwToAz-azToHw)
            //   T^d element [i][j] computed via eigendecomposition
            double decay = Math.pow(1.0 - sumRate, d);

            double pHWtoAZ = piAZ * (1.0 - decay);
            double pAZtoHW = piHW * (1.0 - decay);

            // Clamp transition probabilities away from 0 and 1
            pHWtoAZ = clampProb(pHWtoAZ);
            pAZtoHW = clampProb(pAZtoHW);

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

    private static Double getDoubleOrNull(Map<String, Object> options, String key) {
        if (options == null) return null;
        Object val = options.get(key);
        if (val == null) return null;
        return ((Number) val).doubleValue();
    }
}
