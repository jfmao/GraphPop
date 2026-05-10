package org.graphpop.procedures.pairwise;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

import java.util.Random;

/**
 * Lanczos iteration with full re-orthogonalisation on the branch
 * GRM. Computes the top-K eigenpairs without materialising the full
 * matrix; uses {@link BranchGrmMatVec#apply} as the only access to
 * {@code G}. Memory: {@code O(K · n + n_iter²)}.
 *
 * <p>The k×k tridiagonal subproblem is solved by Apache Commons
 * Math 3 {@link EigenDecomposition} (general symmetric).</p>
 *
 * <p>Reference: Demmel, <em>Applied Numerical Linear Algebra</em>,
 * Algorithm 7.1 (Lanczos with re-orthogonalisation). For symmetric
 * positive-semidefinite GRMs this gives geometric convergence on the
 * top eigenvalues.</p>
 */
public final class BranchGrmLanczos {

    private BranchGrmLanczos() {}

    public static final class Result {
        /** Top-K eigenvalues, sorted in descending order. */
        public final double[] eigenvalues;
        /** Top-K eigenvectors (each length n_samples), aligned with eigenvalues. */
        public final double[][] eigenvectors;
        /** Number of Lanczos iterations actually run. */
        public final int nIter;

        public Result(double[] eigenvalues, double[][] eigenvectors, int nIter) {
            this.eigenvalues = eigenvalues;
            this.eigenvectors = eigenvectors;
            this.nIter = nIter;
        }
    }

    /**
     * Run Lanczos and return the top-{@code k} eigenpairs.
     *
     * @param arg          loaded ARG
     * @param regionStart  bp inclusive
     * @param regionEnd    bp exclusive ({@link Long#MAX_VALUE} for full)
     * @param weightFn     per-branch weight (for conditional G-vector products)
     * @param k            number of top eigenpairs to return
     * @param nIter        Lanczos iteration count (default 3·k, capped at n-1).
     *                     Pass &le; 0 to use the default.
     * @param seed         RNG seed for the starting vector (deterministic output)
     * @param tol          convergence tolerance on β_k (default 1e-12)
     */
    public static Result runTopK(ARG arg, long regionStart, long regionEnd,
                                  BranchWeightFn weightFn,
                                  int k, int nIter, long seed, double tol) {
        final int n = arg.nSamples();
        if (k < 1) throw new IllegalArgumentException("k must be ≥ 1; got " + k);
        if (k >= n) throw new IllegalArgumentException(
                "k (" + k + ") must be < n_samples (" + n + ")");
        final int maxIter = (nIter > 0) ? Math.min(nIter, n - 1)
                                         : Math.min(3 * k, n - 1);
        if (maxIter < k) {
            throw new IllegalArgumentException(
                "n_iter (" + maxIter + ") must be ≥ k (" + k + ")");
        }

        // q_0: random Gaussian, centred, normalised.
        Random rng = new Random(seed);
        double[] q0 = new double[n];
        double sum = 0;
        for (int i = 0; i < n; i++) {
            q0[i] = rng.nextGaussian();
            sum += q0[i];
        }
        double mean = sum / n;
        double sq = 0;
        for (int i = 0; i < n; i++) {
            q0[i] -= mean;
            sq += q0[i] * q0[i];
        }
        double norm = Math.sqrt(sq);
        if (norm < 1e-15) {
            // Degenerate (extremely unlikely for n ≥ 2). Fall back to a
            // simple alternating ±1 sequence centred to zero.
            for (int i = 0; i < n; i++) q0[i] = (i % 2 == 0) ? 1.0 : -1.0;
            norm = Math.sqrt(n);
        }
        for (int i = 0; i < n; i++) q0[i] /= norm;

        double[][] Q = new double[maxIter + 1][];
        Q[0] = q0;
        double[] alpha = new double[maxIter];
        double[] beta = new double[maxIter + 1];   // beta[0] unused

        int kIter = 0;
        for (int iter = 0; iter < maxIter; iter++) {
            double[] w = BranchGrmMatVec.apply(
                    arg, regionStart, regionEnd, weightFn, Q[iter]);
            // alpha[iter] = q[iter] · w
            double a = 0;
            for (int i = 0; i < n; i++) a += Q[iter][i] * w[i];
            alpha[iter] = a;
            // w -= alpha · q[iter]
            for (int i = 0; i < n; i++) w[i] -= a * Q[iter][i];
            // w -= beta[iter] · q[iter-1]   (for iter > 0)
            if (iter > 0) {
                double b = beta[iter];
                for (int i = 0; i < n; i++) w[i] -= b * Q[iter - 1][i];
            }
            // Full re-orthogonalisation against previous Q columns.
            for (int j = 0; j <= iter; j++) {
                double dot = 0;
                for (int i = 0; i < n; i++) dot += w[i] * Q[j][i];
                for (int i = 0; i < n; i++) w[i] -= dot * Q[j][i];
            }
            double bnext = 0;
            for (int i = 0; i < n; i++) bnext += w[i] * w[i];
            bnext = Math.sqrt(bnext);
            beta[iter + 1] = bnext;
            kIter = iter + 1;
            if (bnext < tol) break;
            // q[iter+1] = w / bnext
            double[] qNext = new double[n];
            double inv = 1.0 / bnext;
            for (int i = 0; i < n; i++) qNext[i] = w[i] * inv;
            Q[iter + 1] = qNext;
        }

        // Build the kIter × kIter tridiagonal matrix and diagonalise.
        double[][] T = new double[kIter][kIter];
        for (int i = 0; i < kIter; i++) {
            T[i][i] = alpha[i];
            if (i + 1 < kIter) {
                T[i][i + 1] = beta[i + 1];
                T[i + 1][i] = beta[i + 1];
            }
        }
        RealMatrix Tm = new Array2DRowRealMatrix(T, false);
        EigenDecomposition ed = new EigenDecomposition(Tm);
        double[] evals = ed.getRealEigenvalues();   // not necessarily sorted

        // Sort indices by eigenvalue descending.
        Integer[] order = new Integer[kIter];
        for (int i = 0; i < kIter; i++) order[i] = i;
        java.util.Arrays.sort(order, (x, y) -> Double.compare(evals[y], evals[x]));

        int kReturn = Math.min(k, kIter);
        double[] topEvals = new double[kReturn];
        double[][] topVecs = new double[kReturn][n];

        for (int idx = 0; idx < kReturn; idx++) {
            int j = order[idx];
            topEvals[idx] = evals[j];
            double[] tv = ed.getEigenvector(j).toArray();   // length kIter
            // Project to sample space: pc = Q · tv
            double[] pc = new double[n];
            for (int t = 0; t < kIter; t++) {
                double w = tv[t];
                if (w == 0) continue;
                double[] qt = Q[t];
                for (int i = 0; i < n; i++) pc[i] += w * qt[i];
            }
            // Canonicalise sign: largest-magnitude entry positive.
            int argMax = 0;
            double maxAbs = Math.abs(pc[0]);
            for (int i = 1; i < n; i++) {
                if (Math.abs(pc[i]) > maxAbs) {
                    maxAbs = Math.abs(pc[i]);
                    argMax = i;
                }
            }
            if (pc[argMax] < 0) {
                for (int i = 0; i < n; i++) pc[i] = -pc[i];
            }
            topVecs[idx] = pc;
        }
        return new Result(topEvals, topVecs, kIter);
    }

    /** Convenience overload using the unit weight function. */
    public static Result runTopK(ARG arg, long regionStart, long regionEnd,
                                  int k, int nIter, long seed) {
        return runTopK(arg, regionStart, regionEnd,
                BranchWeightFn.UNIT, k, nIter, seed, 1e-12);
    }
}
