package org.graphpop.procedures.arg;

import org.graphpop.procedures.pairwise.ARG;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

/**
 * Pure-Java per-bin coalescence-rate aggregator. Used by both
 * {@link ArgCoalescenceRateProcedure} (M6) and the demographic
 * inference procedure (M7), which inverts the rate to Ne(t).
 *
 * <p>Implementation mirrors the Python reference at
 * {@code build_egrm_fixture.py:coalescence_rate_reference}
 * bit-for-bit so the existing fixture validates both consumers.</p>
 */
public final class ArgCoalescenceRateComputer {

    private ArgCoalescenceRateComputer() {}

    /**
     * Per-bin coalescent stats over the focal sample set.
     */
    public static final class BinResult {
        public final double timeLo;
        public final double timeHi;
        public final double events;
        public final double pairTime;
        public final double rate;

        BinResult(double timeLo, double timeHi,
                  double events, double pairTime, double rate) {
            this.timeLo = timeLo;
            this.timeHi = timeHi;
            this.events = events;
            this.pairTime = pairTime;
            this.rate = rate;
        }
    }

    /**
     * Compute per-bin coalescent rate for a focal sample set on the
     * given ARG. The bin edges are the row boundaries of the output;
     * {@code bins} must be strictly increasing.
     */
    public static List<BinResult> compute(ARG arg, int[] focalNodes, double[] bins) {
        int nBins = bins.length - 1;
        double[] events = new double[nBins];
        double[] pairTime = new double[nBins];
        double seqLen = arg.sequenceLength;

        long[] bp = arg.breakpointsClipped(0L, arg.sequenceLength);
        for (int k = 0; k < bp.length - 1; k++) {
            long b0 = bp[k];
            long b1 = bp[k + 1];
            if (b1 <= b0) continue;
            double s = (b1 - b0) / seqLen;

            int[] parent = ArgTraversalUtils.marginalTree(arg, b0);
            int[] desc = ArgTraversalUtils.descendantCounts(arg, parent, focalNodes);

            // Build children lists for this marginal tree.
            Map<Integer, List<Integer>> children = new HashMap<>();
            for (int i = 0; i < arg.nNodes; i++) {
                if (parent[i] != -1) {
                    children.computeIfAbsent(parent[i], x -> new ArrayList<>())
                            .add(i);
                }
            }

            // Coalescence events: at each internal node, the focal-pair
            // events here = C(total, 2) - Σ C(c_i, 2).
            for (int u = 0; u < arg.nNodes; u++) {
                List<Integer> cs = children.get(u);
                if (cs == null || cs.size() < 2) continue;
                int total = 0;
                int sumChooseTwo = 0;
                for (int c : cs) {
                    int ck = desc[c];
                    total += ck;
                    sumChooseTwo += ck * (ck - 1) / 2;
                }
                if (total < 2) continue;
                int cross = (total * (total - 1) / 2) - sumChooseTwo;
                if (cross == 0) continue;
                int bi = binIndex(bins, arg.time[u]);
                if (bi >= 0) events[bi] += cross * s;
            }

            // Lineage-pair time: piecewise-constant k_t over node-time grid.
            TreeSet<Double> timeSet = new TreeSet<>();
            for (int i = 0; i < arg.nNodes; i++) {
                if (parent[i] != -1 || children.containsKey(i)) {
                    timeSet.add(arg.time[i]);
                }
            }
            Double[] timeArr = timeSet.toArray(new Double[0]);
            for (int t = 0; t < timeArr.length - 1; t++) {
                double tLo = timeArr[t];
                double tHi = timeArr[t + 1];
                int kT = 0;
                for (int u = 0; u < arg.nNodes; u++) {
                    int p = parent[u];
                    if (p == -1) continue;
                    double ct = arg.time[u];
                    double pt = arg.time[p];
                    if (ct <= tLo && pt >= tHi && desc[u] > 0) {
                        kT++;
                    }
                }
                double pairs = kT * (kT - 1) / 2.0;
                if (pairs == 0) continue;
                for (int bi = 0; bi < nBins; bi++) {
                    double lo = Math.max(bins[bi], tLo);
                    double hi = Math.min(bins[bi + 1], tHi);
                    if (hi > lo) pairTime[bi] += pairs * (hi - lo) * s;
                }
            }
        }

        List<BinResult> out = new ArrayList<>(nBins);
        for (int bi = 0; bi < nBins; bi++) {
            double rate = (pairTime[bi] > 0) ? events[bi] / pairTime[bi] : 0.0;
            out.add(new BinResult(bins[bi], bins[bi + 1],
                    events[bi], pairTime[bi], rate));
        }
        return out;
    }

    /**
     * Validate that {@code bins} is non-empty, strictly monotonically
     * increasing, and has at least 2 edges. Throws
     * {@link IllegalArgumentException} on violation.
     */
    public static void validateBins(double[] bins) {
        if (bins.length < 2) {
            throw new IllegalArgumentException(
                    "time_bins must have at least 2 edges");
        }
        for (int i = 1; i < bins.length; i++) {
            if (bins[i] <= bins[i - 1]) {
                throw new IllegalArgumentException(
                        "time_bins must be strictly increasing");
            }
        }
    }

    private static int binIndex(double[] bins, double t) {
        for (int i = 0; i < bins.length - 1; i++) {
            if (bins[i] <= t && t < bins[i + 1]) return i;
        }
        return -1;
    }
}
