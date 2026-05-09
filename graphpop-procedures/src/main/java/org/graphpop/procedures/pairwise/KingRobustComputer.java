package org.graphpop.procedures.pairwise;

import org.graphpop.procedures.PackedGenotypeReader;

/**
 * KING-robust pairwise kinship from packed 2-bit genotype arrays.
 *
 * <p>Manichaikul et al. 2010, eq. 9:
 * {@code phi_ij = (HetHet - 2 * IBS0) / (2 * min(N_Aa(i), N_Aa(j)))}
 * where {@code HetHet} is the count of variants where both samples are
 * heterozygous, {@code IBS0} is the count where one is hom-ref and the other
 * hom-alt, and {@code N_Aa(s)} is the heterozygous-site count of sample s.
 * Variants where either sample is missing are excluded from all counters.</p>
 *
 * <p>This computer operates directly on {@code byte[]} packed genotype arrays
 * (2 bits per sample) as stored on Variant nodes. It does not require Neo4j;
 * the {@code KingRobustProcedure} facade is responsible for streaming
 * variants out of the graph and feeding them in.</p>
 */
public final class KingRobustComputer {

    private KingRobustComputer() {}

    /**
     * Per-pair accumulator. Public mutable counters keep the loop body
     * branch-free for SIMD-friendly extension later.
     */
    public static final class PairAccumulator {
        public int hetHet;
        public int ibs0;
        public int nAaI;
        public int nAaJ;
        public int nUsed;

        /** Update counters from one variant's two-bit genotypes for samples i and j. */
        public void accumulate(int gtI, int gtJ) {
            if (gtI == PackedGenotypeReader.GT_MISSING
                    || gtJ == PackedGenotypeReader.GT_MISSING) {
                return;
            }
            nUsed++;
            boolean iHet = (gtI == PackedGenotypeReader.GT_HET);
            boolean jHet = (gtJ == PackedGenotypeReader.GT_HET);
            if (iHet) nAaI++;
            if (jHet) nAaJ++;
            if (iHet && jHet) hetHet++;
            // IBS0: opposite homozygotes
            boolean ibs0Pair =
                    (gtI == PackedGenotypeReader.GT_HOM_REF && gtJ == PackedGenotypeReader.GT_HOM_ALT)
                 || (gtI == PackedGenotypeReader.GT_HOM_ALT && gtJ == PackedGenotypeReader.GT_HOM_REF);
            if (ibs0Pair) ibs0++;
        }

        /** KING-robust kinship coefficient. Returns NaN when {@code min(N_Aa)} is zero. */
        public double phi() {
            int minHet = Math.min(nAaI, nAaJ);
            if (minHet == 0) return Double.NaN;
            return (hetHet - 2.0 * ibs0) / (2.0 * minHet);
        }
    }

    /**
     * Symmetric upper-triangular result. {@code phi[i][j]} is stored only for
     * {@code i &lt;= j}; queries are symmetric. Diagonal is always {@code 0.5}
     * by the KING formula (HetHet=N_Aa, IBS0=0).
     */
    public static final class Result {
        private final int n;
        private final double[] phi;
        private final int[] nUsed;

        Result(int n) {
            this.n = n;
            int sz = n * (n + 1) / 2;
            this.phi = new double[sz];
            this.nUsed = new int[sz];
        }

        private int idx(int i, int j) {
            int a = Math.min(i, j);
            int b = Math.max(i, j);
            return b * (b + 1) / 2 + a;
        }

        void put(int i, int j, double phiValue, int nUsedValue) {
            int k = idx(i, j);
            phi[k] = phiValue;
            nUsed[k] = nUsedValue;
        }

        public double phi(int i, int j) { return phi[idx(i, j)]; }
        public int nUsed(int i, int j) { return nUsed[idx(i, j)]; }
        public int nSamples() { return n; }
    }

    /**
     * Compute KING-robust kinship for every pair (including self) over the
     * given packed-genotype variants.
     *
     * @param gtPackedPerVariant one byte[] per variant in the analysis region
     * @param nSamples number of samples in the analysis
     * @param packedIndices map from analysis-position to packed_index in each
     *                      variant's gt_packed; length must equal {@code nSamples}
     * @return symmetric kinship matrix as a {@link Result}
     */
    public static Result computeAllPairs(byte[][] gtPackedPerVariant,
                                         int nSamples,
                                         int[] packedIndices) {
        if (packedIndices.length != nSamples) {
            throw new IllegalArgumentException(
                    "packedIndices length " + packedIndices.length
                  + " != nSamples " + nSamples);
        }
        // Allocate accumulators for upper-triangle pairs (including diagonal).
        int nPairs = nSamples * (nSamples + 1) / 2;
        PairAccumulator[] acc = new PairAccumulator[nPairs];
        for (int p = 0; p < nPairs; p++) acc[p] = new PairAccumulator();

        int[] gtBuf = new int[nSamples];
        for (byte[] gtPacked : gtPackedPerVariant) {
            // Decode all samples once per variant.
            for (int s = 0; s < nSamples; s++) {
                gtBuf[s] = PackedGenotypeReader.genotype(gtPacked, packedIndices[s]);
            }
            // Update each pair accumulator.
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
                result.put(a, b, acc[p].phi(), acc[p].nUsed);
                p++;
            }
        }
        return result;
    }
}
