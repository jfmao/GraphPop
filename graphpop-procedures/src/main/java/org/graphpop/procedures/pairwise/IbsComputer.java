package org.graphpop.procedures.pairwise;

import org.graphpop.procedures.PackedGenotypeReader;

/**
 * Identity-by-state pairwise statistic on packed 2-bit genotype
 * arrays.
 *
 * <p>For samples i, j and variant v:</p>
 * <pre>
 *   ibs(i, j, v) = (2 − |gt_i(v) − gt_j(v)|) / 2     ∈ {0, 0.5, 1}
 * </pre>
 * <p>with {@code gt ∈ {0, 1, 2}} (alt-allele dosage). Variants where
 * either sample is missing are excluded. The pair-level statistic is
 * the mean over called variants:</p>
 * <pre>
 *   IBS_ij = (2·ibs2 + ibs1) / (2 · n_used)
 * </pre>
 * <p>where {@code ibs0} = opposite homozygotes, {@code ibs1} = one
 * allele match (het pair, or het vs homozygote), {@code ibs2} =
 * identical genotypes.</p>
 *
 * <p>Self-pair ({@code i == i}) is always 1. Same packed-array
 * structure as {@link KingRobustComputer}.</p>
 */
public final class IbsComputer {

    private IbsComputer() {}

    public static final class PairAccumulator {
        public int ibs0;   // opposite homozygotes
        public int ibs1;   // one allele match
        public int ibs2;   // identical genotypes
        public int nUsed;  // variants where neither is missing

        /** Update counters for one variant given i and j's 2-bit genotypes. */
        public void accumulate(int gtI, int gtJ) {
            if (gtI == PackedGenotypeReader.GT_MISSING
                    || gtJ == PackedGenotypeReader.GT_MISSING) {
                return;
            }
            nUsed++;
            int diff = Math.abs(gtI - gtJ);
            if (diff == 0) {
                ibs2++;
            } else if (diff == 1) {
                ibs1++;
            } else {
                ibs0++;
            }
        }

        /** Mean IBS in {@code [0, 1]}. NaN when {@code n_used == 0}. */
        public double ibs() {
            if (nUsed == 0) return Double.NaN;
            return (2.0 * ibs2 + ibs1) / (2.0 * nUsed);
        }
    }

    public static final class Result {
        private final int n;
        private final double[] ibs;
        private final int[] nUsed;
        private final int[] ibs0;
        private final int[] ibs2;

        Result(int n) {
            this.n = n;
            int sz = n * (n + 1) / 2;
            this.ibs = new double[sz];
            this.nUsed = new int[sz];
            this.ibs0 = new int[sz];
            this.ibs2 = new int[sz];
        }

        private int idx(int i, int j) {
            int a = Math.min(i, j);
            int b = Math.max(i, j);
            return b * (b + 1) / 2 + a;
        }

        void put(int i, int j, double ibsValue, int nUsedValue,
                 int ibs0Value, int ibs2Value) {
            int k = idx(i, j);
            ibs[k] = ibsValue;
            nUsed[k] = nUsedValue;
            ibs0[k] = ibs0Value;
            ibs2[k] = ibs2Value;
        }

        public double ibs(int i, int j) { return ibs[idx(i, j)]; }
        public int nUsed(int i, int j) { return nUsed[idx(i, j)]; }
        public int ibs0(int i, int j) { return ibs0[idx(i, j)]; }
        public int ibs2(int i, int j) { return ibs2[idx(i, j)]; }
        public int nSamples() { return n; }
    }

    /**
     * Compute IBS for every pair (including self) over packed-genotype
     * variants.
     */
    public static Result computeAllPairs(byte[][] gtPackedPerVariant,
                                          int nSamples,
                                          int[] packedIndices) {
        if (packedIndices.length != nSamples) {
            throw new IllegalArgumentException(
                "packedIndices length " + packedIndices.length
              + " != nSamples " + nSamples);
        }
        int nPairs = nSamples * (nSamples + 1) / 2;
        PairAccumulator[] acc = new PairAccumulator[nPairs];
        for (int p = 0; p < nPairs; p++) acc[p] = new PairAccumulator();

        int[] gtBuf = new int[nSamples];
        for (byte[] gtPacked : gtPackedPerVariant) {
            for (int s = 0; s < nSamples; s++) {
                gtBuf[s] = PackedGenotypeReader.genotype(gtPacked, packedIndices[s]);
            }
            int p = 0;
            for (int b = 0; b < nSamples; b++) {
                for (int a = 0; a <= b; a++) {
                    acc[p++].accumulate(gtBuf[a], gtBuf[b]);
                }
            }
        }

        Result result = new Result(nSamples);
        int p = 0;
        for (int b = 0; b < nSamples; b++) {
            for (int a = 0; a <= b; a++) {
                result.put(a, b, acc[p].ibs(), acc[p].nUsed,
                           acc[p].ibs0, acc[p].ibs2);
                p++;
            }
        }
        return result;
    }
}
