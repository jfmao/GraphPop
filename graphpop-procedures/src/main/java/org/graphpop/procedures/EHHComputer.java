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
    /** Maximum physical gap (bp) before walk is aborted. Matches selscan/scikit-allel default. */
    private static final long DEFAULT_MAX_GAP = 200_000;
    /** Gap scaling threshold (bp). Gaps larger than this are down-weighted. Matches selscan/scikit-allel. */
    private static final long DEFAULT_GAP_SCALE = 20_000;

    private EHHComputer() {}

    /**
     * Compute integrated EHH (iHH) by walking both directions from a focal variant
     * in the in-memory haplotype matrix.
     *
     * <p>Matches selscan/scikit-allel behavior:
     * <ul>
     *   <li>EHH truncated to 0 when below minEHH before final trapezoid step</li>
     *   <li>Walk aborted (returns -1) if any gap exceeds maxGap (200kb)</li>
     *   <li>Gaps larger than gapScale (20kb) are down-weighted by gapScale/gap</li>
     * </ul>
     * </p>
     *
     * @param matrix            dense haplotype matrix
     * @param focalIdx          index of the focal variant in the matrix
     * @param carrierHapIdxs    haplotype indices (into matrix columns) that carry the focal allele
     * @param minEHH            stop when EHH drops below this threshold
     * @param skipMonomorphic   if true, skip variants monomorphic in this population
     *                          (correct for iHS where input is pre-filtered; wrong for
     *                          XP-EHH where both pops must walk the same variant set)
     * @return integrated EHH (area under curve, in bp units), or -1 if walk hit a gap > maxGap
     */
    static double computeIHH(HaplotypeMatrix matrix, int focalIdx,
                             int[] carrierHapIdxs, double minEHH,
                             boolean skipMonomorphic) {
        if (carrierHapIdxs.length < 2) return 0.0;

        long focalPos = matrix.positions[focalIdx];

        double ihhUp = walkDirection(matrix, focalIdx, carrierHapIdxs,
                focalPos, minEHH, -1, skipMonomorphic);
        if (ihhUp < 0) return -1;  // max_gap hit — mark locus invalid

        double ihhDown = walkDirection(matrix, focalIdx, carrierHapIdxs,
                focalPos, minEHH, +1, skipMonomorphic);
        if (ihhDown < 0) return -1;

        return ihhUp + ihhDown;
    }

    /**
     * Overload that skips monomorphic variants (default for iHS).
     */
    static double computeIHH(HaplotypeMatrix matrix, int focalIdx,
                             int[] carrierHapIdxs, double minEHH) {
        return computeIHH(matrix, focalIdx, carrierHapIdxs, minEHH, true);
    }

    /**
     * Overload with default minEHH of 0.05 and skipMonomorphic=true.
     */
    static double computeIHH(HaplotypeMatrix matrix, int focalIdx,
                             int[] carrierHapIdxs) {
        return computeIHH(matrix, focalIdx, carrierHapIdxs, DEFAULT_MIN_EHH, true);
    }

    /**
     * Compute Segregating sites by Length (SL) — the number of SNP positions traversed
     * before EHH decays below minEHH. Used by nSL.
     *
     * <p>Unlike iHH which integrates EHH over physical distance (bp), SL counts the
     * number of variant positions traversed. This makes nSL independent of recombination
     * rate variation across the genome.</p>
     *
     * @param matrix         dense haplotype matrix
     * @param focalIdx       index of the focal variant in the matrix
     * @param carrierHapIdxs haplotype indices that carry the focal allele
     * @param minEHH         stop when EHH drops below this threshold
     * @return total number of SNP positions traversed (upstream + downstream)
     */
    static int computeSL(HaplotypeMatrix matrix, int focalIdx,
                         int[] carrierHapIdxs, double minEHH) {
        if (carrierHapIdxs.length < 2) return 0;

        int slUp = walkDirectionSL(matrix, focalIdx, carrierHapIdxs, minEHH, -1);
        int slDown = walkDirectionSL(matrix, focalIdx, carrierHapIdxs, minEHH, +1);

        return slUp + slDown;
    }

    /**
     * Overload with default minEHH of 0.05.
     */
    static int computeSL(HaplotypeMatrix matrix, int focalIdx,
                         int[] carrierHapIdxs) {
        return computeSL(matrix, focalIdx, carrierHapIdxs, DEFAULT_MIN_EHH);
    }

    /**
     * Walk in one direction from the focal variant, tracking haplotype group splits
     * using flat arrays and incremental sumPairs.
     *
     * <p>Matches selscan/scikit-allel behavior:
     * <ul>
     *   <li>max_gap: if physical gap > 200kb, return -1 (abort entire locus)</li>
     *   <li>gap_scale: if gap > 20kb, scale integration distance by gapScale/gap</li>
     *   <li>EHH truncation: when EHH drops below minEHH, set to 0 before final trapezoid</li>
     * </ul>
     * </p>
     *
     * @param step              -1 for upstream (lower positions), +1 for downstream
     * @param skipMonomorphic   skip variants where AC==0 or AC==AN in this matrix
     * @return integrated EHH, or -1 if a gap > maxGap was encountered
     */
    private static double walkDirection(HaplotypeMatrix matrix, int focalIdx,
                                        int[] carrierHapIdxs,
                                        long focalPos, double minEHH, int step,
                                        boolean skipMonomorphic) {
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
            // Skip variants monomorphic in this population (iHS mode).
            // For iHS, scikit-allel/selscan pre-filter to polymorphic variants,
            // so monomorphic sites shouldn't contribute distance to iHH.
            // For XP-EHH, both populations must walk the SAME variant set
            // (shared polymorphic), so monomorphic-in-one-pop sites must NOT
            // be skipped — they still contribute distance to the trapezoid.
            if (skipMonomorphic) {
                int acHere = matrix.acs[idx];
                if (acHere == 0 || acHere == matrix.ans[idx]) {
                    idx += step;
                    continue;
                }
            }

            long pos = matrix.positions[idx];
            long physGap = Math.abs(pos - prevPos);

            // max_gap: abort entire locus if gap exceeds threshold
            // (matches selscan skipLocus behavior)
            if (physGap > DEFAULT_MAX_GAP) {
                return -1;
            }

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

            // EHH truncation: when EHH drops below threshold, set to 0 before
            // the final trapezoid step. Matches selscan/scikit-allel behavior.
            boolean belowThreshold = ehh < minEHH;
            if (belowThreshold) {
                ehh = 0.0;
            }

            // gap_scale: down-weight integration distance for large gaps
            // effective_dist = physGap * min(1, gapScale / physGap)
            double dist = physGap;
            if (physGap > DEFAULT_GAP_SCALE) {
                dist = DEFAULT_GAP_SCALE;
            }

            // Integrate using trapezoidal rule
            ihh += (prevEHH + ehh) / 2.0 * dist;

            if (belowThreshold) break;

            prevEHH = ehh;
            prevPos = pos;
            idx += step;
        }

        return ihh;
    }

    /**
     * Walk in one direction counting SNP positions traversed (for SL/nSL).
     * Same group-splitting logic as walkDirection, but counts steps instead of
     * integrating physical distance.
     */
    private static int walkDirectionSL(HaplotypeMatrix matrix, int focalIdx,
                                       int[] carrierHapIdxs,
                                       double minEHH, int step) {
        int nHaps = carrierHapIdxs.length;
        if (nHaps < 2) return 0;

        double nPairs = (double) nHaps * (nHaps - 1) / 2.0;

        int[] groupOf = new int[nHaps];
        int nextGroupId = 1;

        int[] groupSize = new int[nHaps];
        groupSize[0] = nHaps;

        double sumPairs = nPairs;

        int[] count0 = new int[nHaps];
        int[] count1 = new int[nHaps];
        int[] seenGroups = new int[nHaps];
        int nSeen;

        int snpCount = 0;

        int idx = focalIdx + step;
        while (idx >= 0 && idx < matrix.nVariants) {
            byte[] hapRow = matrix.haplotypes[idx];

            nSeen = 0;
            for (int h = 0; h < nHaps; h++) {
                int grp = groupOf[h];
                if (count0[grp] == 0 && count1[grp] == 0) {
                    seenGroups[nSeen++] = grp;
                }
                if (hapRow[carrierHapIdxs[h]] == 1) {
                    count1[grp]++;
                } else {
                    count0[grp]++;
                }
            }

            for (int s = 0; s < nSeen; s++) {
                int grp = seenGroups[s];
                int c0 = count0[grp];
                int c1 = count1[grp];

                if (c0 > 0 && c1 > 0) {
                    int oldSize = groupSize[grp];
                    sumPairs -= (double) oldSize * (oldSize - 1) / 2.0;

                    groupSize[grp] = c0;
                    sumPairs += (double) c0 * (c0 - 1) / 2.0;

                    int newGrp = nextGroupId++;
                    groupSize[newGrp] = c1;
                    sumPairs += (double) c1 * (c1 - 1) / 2.0;

                    for (int h = 0; h < nHaps; h++) {
                        if (groupOf[h] == grp && hapRow[carrierHapIdxs[h]] == 1) {
                            groupOf[h] = newGrp;
                        }
                    }
                }

                count0[grp] = 0;
                count1[grp] = 0;
            }

            double ehh = sumPairs / nPairs;

            snpCount++;

            if (ehh < minEHH) break;

            idx += step;
        }

        return snpCount;
    }
}
