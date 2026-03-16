package org.graphpop.procedures;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Garud et al. (2015) haplotype homozygosity statistics.
 *
 * <p>Computes H1, H12, and H2/H1 from haplotype frequencies in a genomic window.
 * These statistics detect hard and soft selective sweeps:</p>
 * <ul>
 *   <li><b>H1</b> — expected haplotype homozygosity = Σ(fi²). High when few
 *       haplotypes dominate (sweep or low diversity).</li>
 *   <li><b>H12</b> — (f1+f2)² + Σ(fi² for i≥3). Elevated by both hard sweeps
 *       (one dominant haplotype) and soft sweeps (two common haplotypes).</li>
 *   <li><b>H2/H1</b> — (H1 - f1²)/H1. Distinguishes hard from soft sweeps:
 *       near 0 for hard sweeps, elevated for soft sweeps.</li>
 * </ul>
 */
final class GarudH {

    private GarudH() {}

    /**
     * Result container for Garud's H statistics.
     */
    static class Result {
        final double h1;
        final double h12;
        final double h2_h1;
        final double hap_diversity;
        final int n_haplotypes;

        Result(double h1, double h12, double h2_h1, double hap_diversity, int n_haplotypes) {
            this.h1 = h1;
            this.h12 = h12;
            this.h2_h1 = h2_h1;
            this.hap_diversity = hap_diversity;
            this.n_haplotypes = n_haplotypes;
        }
    }

    /**
     * Compute Garud's H statistics from a haplotype matrix window.
     *
     * <p>Each sample contributes two haplotypes. Haplotypes are identified by
     * their byte pattern across all variants in the window.</p>
     *
     * @param haplotypes haplotype matrix: [variantIdx][hapIdx], values 0/1
     * @param varStart   first variant index (inclusive)
     * @param varEnd     last variant index (inclusive)
     * @param nHaplotypes total number of haplotypes (2 × nSamples)
     * @return computed H statistics, or null if no haplotypes
     */
    static Result compute(byte[][] haplotypes, int varStart, int varEnd, int nHaplotypes) {
        if (nHaplotypes == 0) return null;

        int nVars = varEnd - varStart + 1;

        // Build haplotype strings and count frequencies
        // Use a Map<HaplotypeKey, Integer> for efficient hashing
        Map<HaplotypeKey, int[]> counts = new HashMap<>();

        for (int h = 0; h < nHaplotypes; h++) {
            byte[] hapPattern = new byte[nVars];
            for (int v = 0; v < nVars; v++) {
                byte[] row = haplotypes[varStart + v];
                hapPattern[v] = (byte)((row[h >> 3] >> (h & 7)) & 1);
            }
            HaplotypeKey key = new HaplotypeKey(hapPattern);
            counts.computeIfAbsent(key, k -> new int[]{0})[0]++;
        }

        int nDistinct = counts.size();

        // Sort frequencies in descending order
        int[] freqCounts = new int[nDistinct];
        int idx = 0;
        for (int[] c : counts.values()) {
            freqCounts[idx++] = c[0];
        }
        Arrays.sort(freqCounts);
        // Reverse to descending
        for (int i = 0; i < nDistinct / 2; i++) {
            int tmp = freqCounts[i];
            freqCounts[i] = freqCounts[nDistinct - 1 - i];
            freqCounts[nDistinct - 1 - i] = tmp;
        }

        // Compute frequencies
        double n = nHaplotypes;
        double[] freq = new double[nDistinct];
        for (int i = 0; i < nDistinct; i++) {
            freq[i] = freqCounts[i] / n;
        }

        // H1 = Σ fi²
        double h1 = 0.0;
        for (double f : freq) h1 += f * f;

        // H12 = (f1 + f2)² + Σ(fi² for i≥3)
        double h12;
        if (nDistinct >= 2) {
            double combined = freq[0] + freq[1];
            h12 = combined * combined;
            for (int i = 2; i < nDistinct; i++) {
                h12 += freq[i] * freq[i];
            }
        } else {
            h12 = h1; // only one haplotype
        }

        // H2/H1 = (H1 - f1²) / H1
        double h2_h1 = h1 > 0 ? (h1 - freq[0] * freq[0]) / h1 : 0.0;

        // Haplotype diversity: Hd = (n/(n-1)) × (1 - H1)
        double hap_diversity = n > 1 ? (n / (n - 1.0)) * (1.0 - h1) : 0.0;

        return new Result(h1, h12, h2_h1, hap_diversity, nDistinct);
    }

    /**
     * Wrapper for byte array to use as HashMap key with proper equals/hashCode.
     */
    private static final class HaplotypeKey {
        private final byte[] data;
        private final int hash;

        HaplotypeKey(byte[] data) {
            this.data = data;
            this.hash = Arrays.hashCode(data);
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof HaplotypeKey)) return false;
            return Arrays.equals(data, ((HaplotypeKey) o).data);
        }

        @Override
        public int hashCode() {
            return hash;
        }
    }
}
