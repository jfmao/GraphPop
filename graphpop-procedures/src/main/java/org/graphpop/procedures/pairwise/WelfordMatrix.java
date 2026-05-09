package org.graphpop.procedures.pairwise;

import java.util.Arrays;

/**
 * Element-wise numerically-stable running mean and variance for an
 * {@code n × n} matrix (Welford 1962). Used by
 * {@link PosteriorBranchGrmProcedure} to aggregate {@code branch_grm}
 * results over an ARG posterior.
 *
 * <p>For each cell {@code (i, j)} the accumulator tracks {@code n_seen},
 * {@code mean}, and {@code M2}. After {@code k} updates:</p>
 * <ul>
 *   <li>{@link #mean()} returns the sample mean.</li>
 *   <li>{@link #stderr()} returns the standard error of the mean,
 *       {@code SE = sqrt(M2 / (k * (k - 1)))}, or {@code NaN} when
 *       {@code k < 2}.</li>
 * </ul>
 *
 * <p>Symmetric inputs produce symmetric outputs by construction
 * (per-cell streaming).</p>
 */
public final class WelfordMatrix {

    private final int n;
    private int nSeen;
    private final double[][] mean;
    private final double[][] m2;

    public WelfordMatrix(int n) {
        this.n = n;
        this.mean = new double[n][n];
        this.m2 = new double[n][n];
    }

    /** Number of updates absorbed so far. */
    public int nSeen() { return nSeen; }
    public int size() { return n; }

    /** Absorb one full {@code n × n} matrix into the running statistics. */
    public void update(double[][] x) {
        if (x.length != n) {
            throw new IllegalArgumentException(
                "row count mismatch: expected " + n + ", got " + x.length);
        }
        nSeen++;
        double inv = 1.0 / nSeen;
        for (int i = 0; i < n; i++) {
            if (x[i].length != n) {
                throw new IllegalArgumentException(
                    "row " + i + " column count mismatch");
            }
            double[] mRow = mean[i];
            double[] m2Row = m2[i];
            double[] xRow = x[i];
            for (int j = 0; j < n; j++) {
                double delta = xRow[j] - mRow[j];
                mRow[j] += delta * inv;
                double delta2 = xRow[j] - mRow[j];
                m2Row[j] += delta * delta2;
            }
        }
    }

    /** Defensive copy of the running mean matrix. */
    public double[][] mean() {
        double[][] out = new double[n][];
        for (int i = 0; i < n; i++) out[i] = Arrays.copyOf(mean[i], n);
        return out;
    }

    /**
     * Standard error of the mean per cell:
     * {@code SE = sqrt(M2 / (n * (n - 1)))}. Returns {@code NaN} per
     * cell when {@code n_seen < 2}.
     */
    public double[][] stderr() {
        double[][] se = new double[n][n];
        if (nSeen < 2) {
            for (int i = 0; i < n; i++) Arrays.fill(se[i], Double.NaN);
            return se;
        }
        double denom = (double) nSeen * (nSeen - 1);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double v = m2[i][j] / denom;
                se[i][j] = (v < 0.0) ? 0.0 : Math.sqrt(v);
            }
        }
        return se;
    }

    /** Sample standard deviation (ddof = 1) per cell; {@code NaN} when n_seen < 2. */
    public double[][] sampleStd() {
        double[][] sd = new double[n][n];
        if (nSeen < 2) {
            for (int i = 0; i < n; i++) Arrays.fill(sd[i], Double.NaN);
            return sd;
        }
        double denom = (double) (nSeen - 1);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double v = m2[i][j] / denom;
                sd[i][j] = (v < 0.0) ? 0.0 : Math.sqrt(v);
            }
        }
        return sd;
    }
}
