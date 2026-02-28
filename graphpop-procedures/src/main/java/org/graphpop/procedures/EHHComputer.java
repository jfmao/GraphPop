package org.graphpop.procedures;

/**
 * Computes Extended Haplotype Homozygosity (EHH) decay from a focal variant
 * using a dense in-memory {@link HaplotypeMatrix}.
 *
 * <p>EHH measures the probability that two randomly chosen haplotypes carrying
 * a given allele are identical by state from the focal site to a given distance.
 * This implementation walks position-sorted arrays — no database access during
 * computation, enabling parallel execution.</p>
 *
 * <p>Optimized hot loop: flat arrays instead of HashMaps, incremental sumPairs
 * tracking via delta on group splits. O(nGroupsSplit) per step instead of
 * O(totalGroups).</p>
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
     * Walk in one direction from the focal variant, tracking haplotype group splits
     * using flat arrays and incremental sumPairs.
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

        // Flat arrays for group tracking — reused across steps
        // Max possible groups = nHaps (each haplotype in its own group)
        int[] groupSize = new int[nHaps];
        groupSize[0] = nHaps;  // all start in group 0

        // Incremental sumPairs: track sum of C(k_i, 2) across all groups
        // Initially one group of size nHaps → sumPairs = C(nHaps, 2)
        double sumPairs = nPairs;

        // Temporary arrays for per-step group splitting
        int[] count0 = new int[nHaps];  // count of allele=0 per group this step
        int[] count1 = new int[nHaps];  // count of allele=1 per group this step
        // Track which groups were seen this step (for reset)
        int[] seenGroups = new int[nHaps];
        int nSeen;

        double prevEHH = 1.0;
        long prevPos = focalPos;
        double ihh = 0.0;

        // Walk the position-sorted array
        int idx = focalIdx + step;
        while (idx >= 0 && idx < matrix.nVariants) {
            long pos = matrix.positions[idx];
            byte[] hapRow = matrix.haplotypes[idx];

            // Count alleles per group at this site
            nSeen = 0;
            for (int h = 0; h < nHaps; h++) {
                int grp = groupOf[h];
                // Track this group if first time seeing it
                if (count0[grp] == 0 && count1[grp] == 0) {
                    seenGroups[nSeen++] = grp;
                }
                if (hapRow[carrierHapIdxs[h]] == 1) {
                    count1[grp]++;
                } else {
                    count0[grp]++;
                }
            }

            // Split heterogeneous groups and update sumPairs incrementally
            for (int s = 0; s < nSeen; s++) {
                int grp = seenGroups[s];
                int c0 = count0[grp];
                int c1 = count1[grp];

                if (c0 > 0 && c1 > 0) {
                    // Group is heterogeneous — split it
                    int oldSize = groupSize[grp];

                    // Remove old group's pair contribution
                    sumPairs -= (double) oldSize * (oldSize - 1) / 2.0;

                    // Existing group keeps allele=0 haplotypes
                    groupSize[grp] = c0;
                    sumPairs += (double) c0 * (c0 - 1) / 2.0;

                    // New group gets allele=1 haplotypes
                    int newGrp = nextGroupId++;
                    groupSize[newGrp] = c1;
                    sumPairs += (double) c1 * (c1 - 1) / 2.0;

                    // Reassign haplotypes with allele=1 to new group
                    for (int h = 0; h < nHaps; h++) {
                        if (groupOf[h] == grp && hapRow[carrierHapIdxs[h]] == 1) {
                            groupOf[h] = newGrp;
                        }
                    }
                }

                // Reset counts for this group
                count0[grp] = 0;
                count1[grp] = 0;
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
