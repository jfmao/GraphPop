package org.graphpop.procedures.pairwise;

import java.util.Arrays;
import java.util.BitSet;

/**
 * Algorithm V (Tang &amp; Chiang 2025): matrix-vector form of
 * {@code branch_grm} that computes {@code G · v} without
 * materialising {@code G}, in {@code O(T · N log N)} per
 * multiplication.
 *
 * <p>Per branch with descendants {@code D ⊂ samples}, the contribution
 * to {@code G v} is</p>
 * <pre>
 *   (G v)_i  +=  α · 𝟙{i ∈ D} · Σ_{j ∈ D} v_j
 * </pre>
 * <p>with {@code α = mu / (p · (1 − p))} as in the full-matrix
 * kernel. Per-branch cost shrinks from {@code O(|D|²)} to
 * {@code O(|D|)}; biobank-scale (n &gt; 5 000 haplotypes) cohorts
 * become tractable without ever building the n × n matrix.</p>
 *
 * <p>Centring matches the full-matrix path
 * (cf. {@link BranchGrmComputer}): the doubly-centred matrix
 * {@code B = (M / total_mu) − colMeans − rowMeans} satisfies
 * {@code B v = (M / total_mu) v − mean(M v) · 1} for any zero-mean
 * vector. We pre-centre the input vector by subtracting
 * {@code mean(v)}, accumulate {@code (M / total_mu) v}, then
 * subtract {@code mean(result) · 1} as the final step.</p>
 */
public final class BranchGrmMatVec {

    private BranchGrmMatVec() {}

    /**
     * Compute {@code G · v} where {@code G} is the (un-conditioned)
     * branch GRM. {@code v} length must equal {@code arg.nSamples()}.
     */
    public static double[] apply(ARG arg, long regionStart, long regionEnd,
                                  double[] v) {
        return apply(arg, regionStart, regionEnd, BranchWeightFn.UNIT, v);
    }

    /**
     * Compute {@code G · v} with a per-branch weight function (matches
     * step-4 conditional predicates).
     */
    public static double[] apply(ARG arg, long regionStart, long regionEnd,
                                  BranchWeightFn weightFn, double[] v) {
        final int n = arg.nSamples();
        if (v.length != n) {
            throw new IllegalArgumentException(
                "vector length " + v.length + " != n_samples " + n);
        }
        if (n == 0) return new double[0];

        // Pre-centre v: subtract mean(v) so that Σ v_centred = 0.
        double meanV = 0.0;
        for (double x : v) meanV += x;
        meanV /= n;
        final double[] vc = new double[n];
        for (int i = 0; i < n; i++) vc[i] = v[i] - meanV;

        final long lo = Math.max(regionStart, 0L);
        final long hi = (regionEnd == Long.MAX_VALUE) ? arg.sequenceLength
                : Math.min(regionEnd, arg.sequenceLength);
        if (hi <= lo) return new double[n];

        final long[] breakpoints = arg.breakpointsClipped(lo, hi);
        if (breakpoints.length < 2) return new double[n];

        final int nNodes = arg.nNodes;
        final int[] parent = new int[nNodes];
        final BitSet[] descMask = new BitSet[nNodes];
        for (int i = 0; i < nNodes; i++) descMask[i] = new BitSet(n);

        final double[] gv = new double[n];
        double totalMu = 0.0;

        for (int k = 0; k < breakpoints.length - 1; k++) {
            final long b0 = breakpoints[k];
            final long b1 = breakpoints[k + 1];
            final long intervalLen = b1 - b0;
            if (intervalLen <= 0) continue;

            Arrays.fill(parent, -1);
            for (int e = 0; e < arg.nEdges; e++) {
                if (arg.edgeStart[e] <= b0 && arg.edgeEnd[e] > b0) {
                    parent[arg.edgeChild[e]] = arg.edgeParent[e];
                }
            }
            for (int i = 0; i < nNodes; i++) descMask[i].clear();
            for (int sBit = 0; sBit < n; sBit++) {
                int cur = arg.sampleNodes[sBit];
                while (cur != -1) {
                    descMask[cur].set(sBit);
                    cur = parent[cur];
                }
            }

            for (int c = 0; c < nNodes; c++) {
                final int p = parent[c];
                if (p == -1) continue;
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
                totalMu += mu;

                final double pFreq = (double) nDesc / (double) n;
                final double alpha = mu / (pFreq * (1.0 - pFreq));

                // s = Σ_{j ∈ D} v_centred[j]   -- O(|D|)
                double s = 0.0;
                for (int b = mask.nextSetBit(0); b >= 0; b = mask.nextSetBit(b + 1)) {
                    s += vc[b];
                }
                final double increment = alpha * s;
                if (increment == 0.0) continue;
                // gv[i] += alpha · s · 𝟙{i ∈ D}   -- O(|D|)
                for (int b = mask.nextSetBit(0); b >= 0; b = mask.nextSetBit(b + 1)) {
                    gv[b] += increment;
                }
            }
        }

        if (totalMu == 0.0) return gv;
        // Normalise by total_mu and centre on the result side.
        double meanGv = 0.0;
        for (int i = 0; i < n; i++) {
            gv[i] /= totalMu;
            meanGv += gv[i];
        }
        meanGv /= n;
        for (int i = 0; i < n; i++) gv[i] -= meanGv;
        return gv;
    }

    /**
     * Apply {@code G} to a block of {@code k} vectors. Each input
     * column becomes one output column. Used for stacked Lanczos /
     * randomised SVD / multi-trait HE-regression.
     */
    public static double[][] applyBlock(ARG arg, long regionStart, long regionEnd,
                                         BranchWeightFn weightFn,
                                         double[][] vectors) {
        if (vectors.length == 0) return new double[0][];
        int k = vectors.length;
        double[][] out = new double[k][];
        for (int i = 0; i < k; i++) {
            out[i] = apply(arg, regionStart, regionEnd, weightFn, vectors[i]);
        }
        return out;
    }
}
