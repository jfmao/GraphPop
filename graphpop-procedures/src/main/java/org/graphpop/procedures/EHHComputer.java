package org.graphpop.procedures;

import java.util.*;

/**
 * Computes Extended Haplotype Homozygosity (EHH) decay from a focal variant
 * using a dense in-memory {@link HaplotypeMatrix}.
 *
 * <p>EHH measures the probability that two randomly chosen haplotypes carrying
 * a given allele are identical by state from the focal site to a given distance.
 * This implementation walks position-sorted arrays — no database access during
 * computation, enabling parallel execution.</p>
 *
 * <p>Shared by {@link IHSProcedure} and {@link XPEHHProcedure}.</p>
 */
final class EHHComputer {

    private static final double DEFAULT_MIN_EHH = 0.05;

    private EHHComputer() {}

    /**
     * Compute integrated EHH (iHH) by walking both directions from a focal variant
     * in the in-memory haplotype matrix.
     *
     * @param matrix         dense haplotype matrix
     * @param focalIdx       index of the focal variant in the matrix
     * @param carrierHapIdxs haplotype indices (into matrix columns) that carry the focal allele
     * @param minEHH         stop when EHH drops below this threshold
     * @return integrated EHH (area under curve, in bp units)
     */
    static double computeIHH(HaplotypeMatrix matrix, int focalIdx,
                             int[] carrierHapIdxs, double minEHH) {
        if (carrierHapIdxs.length < 2) return 0.0;

        long focalPos = matrix.positions[focalIdx];

        double ihhUp = walkDirection(matrix, focalIdx, carrierHapIdxs,
                focalPos, minEHH, -1);
        double ihhDown = walkDirection(matrix, focalIdx, carrierHapIdxs,
                focalPos, minEHH, +1);

        return ihhUp + ihhDown;
    }

    /**
     * Overload with default minEHH of 0.05.
     */
    static double computeIHH(HaplotypeMatrix matrix, int focalIdx,
                             int[] carrierHapIdxs) {
        return computeIHH(matrix, focalIdx, carrierHapIdxs, DEFAULT_MIN_EHH);
    }

    /**
     * Walk in one direction from the focal variant, tracking haplotype group splits.
     *
     * @param step -1 for upstream (lower positions), +1 for downstream
     */
    private static double walkDirection(HaplotypeMatrix matrix, int focalIdx,
                                        int[] carrierHapIdxs,
                                        long focalPos, double minEHH, int step) {
        int nHaps = carrierHapIdxs.length;
        if (nHaps < 2) return 0.0;

        double nPairs = (double) nHaps * (nHaps - 1) / 2.0;

        // Each haplotype starts in group 0 (all identical)
        int[] groupOf = new int[nHaps];
        int nextGroupId = 1;

        double prevEHH = 1.0;
        long prevPos = focalPos;
        double ihh = 0.0;

        // Walk the position-sorted array
        int idx = focalIdx + step;
        while (idx >= 0 && idx < matrix.nVariants) {
            long pos = matrix.positions[idx];
            byte[] hapRow = matrix.haplotypes[idx];

            // Split groups by allele at this site
            // Use a compact approach: for each group, track sub-groups by allele
            // Since alleles are 0 or 1, we only need two sub-lists per group
            Map<Integer, int[]> groupCounts = new HashMap<>(); // groupId → [count0, count1]
            int[] alleles = new int[nHaps]; // allele for each carrier haplotype

            for (int h = 0; h < nHaps; h++) {
                int allele = hapRow[carrierHapIdxs[h]];
                alleles[h] = allele;
                int grp = groupOf[h];
                int[] counts = groupCounts.computeIfAbsent(grp, k -> new int[2]);
                counts[allele]++;
            }

            // Only split groups that are heterogeneous at this site
            for (Map.Entry<Integer, int[]> entry : groupCounts.entrySet()) {
                int grp = entry.getKey();
                int[] counts = entry.getValue();
                if (counts[0] > 0 && counts[1] > 0) {
                    // This group needs splitting — assign new group ID to allele=1
                    int newGrp = nextGroupId++;
                    for (int h = 0; h < nHaps; h++) {
                        if (groupOf[h] == grp && alleles[h] == 1) {
                            groupOf[h] = newGrp;
                        }
                    }
                }
            }

            // Compute EHH = sum(nCr(k_i, 2)) / nCr(n, 2)
            // Count group sizes
            Map<Integer, Integer> groupSizes = new HashMap<>();
            for (int h = 0; h < nHaps; h++) {
                groupSizes.merge(groupOf[h], 1, Integer::sum);
            }

            double sumPairs = 0.0;
            for (int k : groupSizes.values()) {
                if (k >= 2) {
                    sumPairs += (double) k * (k - 1) / 2.0;
                }
            }
            double ehh = sumPairs / nPairs;

            // Integrate using trapezoidal rule
            long dist = Math.abs(pos - prevPos);
            ihh += (prevEHH + ehh) / 2.0 * dist;

            if (ehh < minEHH) break;

            prevEHH = ehh;
            prevPos = pos;
            idx += step;
        }

        return ihh;
    }
}
