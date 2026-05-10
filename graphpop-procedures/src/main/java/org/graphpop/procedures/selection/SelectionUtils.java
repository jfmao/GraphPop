package org.graphpop.procedures.selection;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Shared helpers for the M8 selection-scan procedures.
 *
 * <p>Provides:</p>
 *
 * <ul>
 *   <li>Welford running mean / variance over a stream of doubles.</li>
 *   <li>Logit-spaced frequency binning (N bins evenly spaced in
 *       {@code logit(f)} between {@code logit(min_f)} and
 *       {@code logit(max_f)}).</li>
 *   <li>Per-bin (mean, std) → per-variant z-score helper.</li>
 * </ul>
 */
public final class SelectionUtils {

    private SelectionUtils() {}

    /**
     * Welford online accumulator for mean and population standard
     * deviation over a stream of doubles. Numerically stable for
     * large streams.
     */
    public static final class Welford {
        private long n = 0L;
        private double mean = 0.0;
        private double m2 = 0.0;

        public void add(double x) {
            n++;
            double delta = x - mean;
            mean += delta / n;
            double delta2 = x - mean;
            m2 += delta * delta2;
        }

        public long n() { return n; }
        public double mean() { return mean; }

        /** Sample standard deviation (n-1 denominator). NaN when n < 2. */
        public double sd() {
            if (n < 2) return Double.NaN;
            return Math.sqrt(m2 / (n - 1));
        }
    }

    /**
     * Build {@code nBins} logit-spaced edges between {@code minFreq}
     * and {@code maxFreq} (both strictly in (0, 1)). Returns
     * {@code nBins + 1} edges; the bin {@code i} is
     * {@code [edges[i], edges[i + 1])}. Outside-range frequencies
     * fall in the first or last bin.
     */
    public static double[] logitFreqEdges(int nBins, double minFreq, double maxFreq) {
        if (nBins < 1) {
            throw new IllegalArgumentException("nBins must be >= 1");
        }
        if (!(minFreq > 0.0 && minFreq < 1.0)) {
            throw new IllegalArgumentException("minFreq must be in (0, 1)");
        }
        if (!(maxFreq > minFreq && maxFreq < 1.0)) {
            throw new IllegalArgumentException(
                    "maxFreq must be in (minFreq, 1)");
        }
        double lo = logit(minFreq);
        double hi = logit(maxFreq);
        double[] edges = new double[nBins + 1];
        for (int i = 0; i <= nBins; i++) {
            double t = i / (double) nBins;
            edges[i] = invLogit(lo + t * (hi - lo));
        }
        // Force exact endpoints to guard against tiny float drift.
        edges[0] = minFreq;
        edges[nBins] = maxFreq;
        return edges;
    }

    /**
     * Index of the bin containing {@code freq}. Frequencies below
     * {@code edges[0]} fall in bin 0; frequencies at or above
     * {@code edges[nBins]} fall in bin {@code nBins - 1}.
     */
    public static int binIndex(double[] edges, double freq) {
        int n = edges.length - 1;
        if (freq < edges[0]) return 0;
        if (freq >= edges[n]) return n - 1;
        // Binary-search; edges sorted ascending.
        int lo = 0, hi = n;
        while (lo < hi - 1) {
            int mid = (lo + hi) >>> 1;
            if (edges[mid] <= freq) lo = mid;
            else hi = mid;
        }
        return lo;
    }

    /**
     * z-score of {@code x} given the bin's running statistics. Returns
     * 0 when the bin has fewer than 2 samples or zero variance, since
     * the z-score is undefined / would inflate to ±∞.
     */
    public static double zScore(double x, Welford bin) {
        double sd = bin.sd();
        if (Double.isNaN(sd) || sd <= 0.0) return 0.0;
        return (x - bin.mean()) / sd;
    }

    /**
     * Sliding-window edges over {@code [0, sequenceLength)} given
     * {@code windowSize} and {@code step}. The last window is
     * truncated to fit the sequence end; windows shorter than
     * {@code windowSize / 2} are dropped.
     */
    public static long[][] slidingWindows(long sequenceLength,
                                           long windowSize,
                                           long step) {
        if (windowSize <= 0 || step <= 0 || sequenceLength <= 0) {
            return new long[0][];
        }
        List<long[]> out = new ArrayList<>();
        long start = 0L;
        while (start < sequenceLength) {
            long end = Math.min(start + windowSize, sequenceLength);
            if (end - start >= windowSize / 2) {
                out.add(new long[]{start, end});
            }
            start += step;
        }
        return out.toArray(new long[0][]);
    }

    private static double logit(double p) { return Math.log(p / (1.0 - p)); }

    private static double invLogit(double x) {
        if (x >= 0) {
            double e = Math.exp(-x);
            return 1.0 / (1.0 + e);
        } else {
            double e = Math.exp(x);
            return e / (1.0 + e);
        }
    }
}
