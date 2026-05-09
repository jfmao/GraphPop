package org.graphpop.procedures.pairwise;

import java.util.Arrays;
import java.util.BitSet;

/**
 * Branch GRM (eGRM) computation on an in-memory {@link ARG}, matching
 * the node-centric algorithm from {@code egrm.varGRM(ts)}
 * (Fan, Mancuso & Chiang 2022 <i>AJHG</i>).
 *
 * <p>For each marginal tree interval {@code [b_k, b_{k+1}]} and each
 * non-root, non-trivial node {@code c} with descendant-sample count
 * {@code n_c ∈ (0, N)}:</p>
 *
 * <pre>
 *   mu_c     = interval_length * branch_length(c) * 1e-8
 *   p_c      = n_c / N
 *   w_c      = 1 / (p_c * (1 - p_c))
 *   egrm[i,j] += mu_c * w_c   for every (i, j) in descendants(c) × descendants(c)
 *   total_mu += mu_c
 * </pre>
 *
 * <p>After all intervals, the matrix is divided by {@code total_mu} and
 * double-centered (subtract column means, then subtract row means of the
 * result), matching {@code egrm.varGRM}'s output exactly.</p>
 *
 * <p>Soft cap: the full matrix is materialised as {@code double[n][n]},
 * which scales as {@code O(n²)} memory. Above the runtime threshold
 * (default {@value #DEFAULT_FULL_MATRIX_CAP}) the kernel throws and
 * recommends {@code branch_grm_apply} (Algorithm V matrix–vector form;
 * M4.1 step 6) for biobank-scale cohorts.</p>
 */
public final class BranchGrmComputer {

    /**
     * Soft cap on {@code n_samples} for the full-matrix path.
     * 5 000 haplotypes ⇒ 25 M cells × 8 bytes ≈ 200 MB. Set higher
     * via {@link #compute(ARG, long, long, BranchWeightFn, int)} if
     * the host has more heap.
     */
    public static final int DEFAULT_FULL_MATRIX_CAP = 5_000;

    private BranchGrmComputer() {}

    public static final class Result {
        public final int n;
        public final double[][] matrix;
        public final double totalMu;

        Result(int n, double[][] matrix, double totalMu) {
            this.n = n;
            this.matrix = matrix;
            this.totalMu = totalMu;
        }

        public double bij(int i, int j) { return matrix[i][j]; }
    }

    /**
     * Compute the unconditional eGRM over the region {@code [regionStart,
     * regionEnd)} (bp). Pass {@code 0} and {@link Long#MAX_VALUE} for the
     * full ARG. Equivalent to
     * {@link #compute(ARG, long, long, BranchWeightFn)} with
     * {@link BranchWeightFn#UNIT}.
     */
    public static Result compute(ARG arg, long regionStart, long regionEnd) {
        return compute(arg, regionStart, regionEnd, BranchWeightFn.UNIT);
    }

    /**
     * Compute the conditional eGRM with a per-branch weight multiplier
     * {@code weightFn}. Used by step 4 conditional predicates
     * ({@code restrict_to_pathway}, {@code mutation_filter},
     * {@code time_window}). The unit weight reproduces the unconditional
     * matrix exactly.
     */
    public static Result compute(ARG arg, long regionStart, long regionEnd,
                                 BranchWeightFn weightFn) {
        return compute(arg, regionStart, regionEnd, weightFn, DEFAULT_FULL_MATRIX_CAP);
    }

    /**
     * Same as {@link #compute(ARG, long, long, BranchWeightFn)} but
     * with a configurable soft cap on sample count for the full
     * matrix path. Use {@link BranchGrmMatVec} for biobank scale
     * (matrix-free).
     */
    public static Result compute(ARG arg, long regionStart, long regionEnd,
                                 BranchWeightFn weightFn,
                                 int fullMatrixCap) {
        final int n = arg.nSamples();
        if (n > fullMatrixCap) {
            throw new IllegalArgumentException(
                "BranchGrmComputer requires n_samples <= " + fullMatrixCap
              + " for the full-matrix path; got " + n
              + ". Use graphpop.kinship.branch_grm_apply (Algorithm V) "
              + "for larger cohorts.");
        }
        if (n == 0) {
            return new Result(0, new double[0][0], 0.0);
        }

        // Map packed node index -> sample bit (or -1 for non-sample).
        final int nNodes = arg.nNodes;
        final int[] packedToSampleBit = new int[nNodes];
        Arrays.fill(packedToSampleBit, -1);
        for (int k = 0; k < n; k++) {
            packedToSampleBit[arg.sampleNodes[k]] = k;
        }

        final long lo = Math.max(regionStart, 0L);
        final long hi = (regionEnd == Long.MAX_VALUE) ? arg.sequenceLength
                                                      : Math.min(regionEnd, arg.sequenceLength);
        if (hi <= lo) {
            return new Result(n, new double[n][n], 0.0);
        }

        final long[] breakpoints = arg.breakpointsClipped(lo, hi);
        if (breakpoints.length < 2) {
            return new Result(n, new double[n][n], 0.0);
        }

        final double[][] egrm = new double[n][n];
        double totalMu = 0.0;

        // Reusable per-tree buffers.
        final int[] parent = new int[nNodes];
        final BitSet[] descMask = new BitSet[nNodes];
        for (int i = 0; i < nNodes; i++) descMask[i] = new BitSet(n);

        for (int k = 0; k < breakpoints.length - 1; k++) {
            final long b0 = breakpoints[k];
            final long b1 = breakpoints[k + 1];
            final long intervalLen = b1 - b0;
            if (intervalLen <= 0) continue;

            // Build the marginal-tree parent[] from edges active at pivot b0.
            Arrays.fill(parent, -1);
            for (int e = 0; e < arg.nEdges; e++) {
                if (arg.edgeStart[e] <= b0 && arg.edgeEnd[e] > b0) {
                    parent[arg.edgeChild[e]] = arg.edgeParent[e];
                }
            }

            // Compute descendant-sample-set per node by walking each sample
            // up to its (transitive) root, ORing its bit at every ancestor.
            for (int i = 0; i < nNodes; i++) descMask[i].clear();
            for (int sBit = 0; sBit < n; sBit++) {
                int cur = arg.sampleNodes[sBit];
                while (cur != -1) {
                    descMask[cur].set(sBit);
                    cur = parent[cur];
                }
            }

            // Per non-root non-trivial node, accumulate into egrm.
            for (int c = 0; c < nNodes; c++) {
                final int p = parent[c];
                if (p == -1) continue;  // root
                final BitSet mask = descMask[c];
                final int nDesc = mask.cardinality();
                if (nDesc == 0 || nDesc == n) continue;

                final double branchLen = arg.time[p] - arg.time[c];
                if (branchLen <= 0.0) continue;

                final double w = weightFn.weight(
                        arg.tskitNodeId[p], arg.tskitNodeId[c],
                        arg.time[p], arg.time[c], b0, b1);
                if (w <= 0.0) continue;

                final double mu = intervalLen * branchLen * w * 1e-8;
                final double pFreq = (double) nDesc / (double) n;
                final double weight = mu / (pFreq * (1.0 - pFreq));

                // Materialise descendant indices once.
                final int[] desc = new int[nDesc];
                int di = 0;
                for (int b = mask.nextSetBit(0); b >= 0; b = mask.nextSetBit(b + 1)) {
                    desc[di++] = b;
                }

                // Outer-product update on the (descendants x descendants) sub-matrix.
                for (int i = 0; i < nDesc; i++) {
                    double[] r = egrm[desc[i]];
                    for (int j = 0; j < nDesc; j++) {
                        r[desc[j]] += weight;
                    }
                }
                totalMu += mu;
            }
        }

        if (totalMu == 0.0) {
            return new Result(n, egrm, 0.0);
        }

        // Normalize by total_mu.
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                egrm[i][j] /= totalMu;
            }
        }

        // Double-center matching numpy:
        //   egrm -= egrm.mean(axis=0)               # subtract column means
        //   egrm -= egrm.mean(axis=1, keepdims=True) # subtract row means OF THE RESULT
        final double[] colMeans = new double[n];
        for (int j = 0; j < n; j++) {
            double s = 0.0;
            for (int i = 0; i < n; i++) s += egrm[i][j];
            colMeans[j] = s / n;
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                egrm[i][j] -= colMeans[j];
            }
        }
        final double[] rowMeans = new double[n];
        for (int i = 0; i < n; i++) {
            double s = 0.0;
            for (int j = 0; j < n; j++) s += egrm[i][j];
            rowMeans[i] = s / n;
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                egrm[i][j] -= rowMeans[i];
            }
        }

        return new Result(n, egrm, totalMu);
    }
}
