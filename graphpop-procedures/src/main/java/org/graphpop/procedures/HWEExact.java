package org.graphpop.procedures;

/**
 * Hardy–Weinberg Equilibrium exact test following Wigginton et al. (2005).
 *
 * <p>Implements the mid-p exact test, which enumerates all possible heterozygote
 * counts for a fixed allele count and total sample count, computes the probability
 * of each configuration under HWE, and sums the probabilities of configurations
 * as or more extreme than the observed one.</p>
 *
 * <p>Reference: Wigginton JE, Cutler DJ, Abecasis GR. A note on exact tests of
 * Hardy-Weinberg equilibrium. Am J Hum Genet. 2005;76(5):887-93.</p>
 */
final class HWEExact {

    private HWEExact() {}

    /**
     * Compute the HWE exact test mid-p value.
     *
     * @param nHet      observed heterozygote count
     * @param nHomMinor observed count of the minor-allele homozygote
     * @param nTotal    total number of diploid individuals
     * @return mid-p value; values near 0 indicate departure from HWE
     */
    static double hweExactMidP(int nHet, int nHomMinor, int nTotal) {
        if (nTotal <= 0) return 1.0;

        // Derive allele counts
        int nHomMajor = nTotal - nHet - nHomMinor;
        if (nHomMajor < 0) return 1.0;

        int minorCount = 2 * nHomMinor + nHet;
        int majorCount = 2 * nHomMajor + nHet;
        int nAlleles = 2 * nTotal;

        if (minorCount > majorCount) {
            // Swap so minor is the less frequent allele
            int tmp = minorCount;
            minorCount = majorCount;
            majorCount = tmp;
        }

        if (minorCount == 0) return 1.0;

        // The possible het counts have the same parity as minorCount
        // and range from 0 (or 1) to minorCount
        // We compute probabilities relative to the mode using ratios

        // Find the expected het count under HWE (the mode of the distribution)
        // hetMax = minorCount if minorCount <= nTotal, else nAlleles - minorCount (bounded)
        int hetMax = minorCount;
        if (hetMax > nTotal) hetMax = nTotal;
        // het values must have same parity as minorCount
        // Adjust hetMax to have correct parity
        if ((hetMax & 1) != (minorCount & 1)) hetMax--;

        // Use log probabilities to avoid overflow
        // P(het_i) is proportional to:
        //   nTotal! / (het_i! * homMinor_i! * homMajor_i!) * (minorCount! * majorCount!) / nAlleles!
        // We compute log-probabilities relative to the observed count.

        // Actually, use the standard iterative approach: compute probs via recurrence
        // P(het - 2) / P(het) = het * (het - 1) / ((homMinor + 1) * (homMajor + 1) * 4)
        // where homMinor = (minorCount - het) / 2, homMajor = (majorCount - het) / 2

        // Compute probabilities for all valid het counts
        int nSteps = hetMax / 2 + 1;  // number of possible het values (step by 2)
        double[] probs = new double[nSteps];

        // Start at observed and normalize later
        // Index mapping: het = minorCount - 2*i when stepping from max to 0
        // Better: index het values as het = minorCount % 2, minorCount % 2 + 2, ..., hetMax
        // Let's enumerate from het=0 (or 1) stepping by 2

        int startHet = minorCount % 2;  // 0 if minor even, 1 if odd

        // Compute using ratios starting from the midpoint
        // Find index of observed het
        int obsIdx = (nHet - startHet) / 2;
        if (obsIdx < 0 || obsIdx >= nSteps) return 1.0;  // observed het impossible

        probs[obsIdx] = 1.0;  // reference probability

        // Walk upward from observed het
        for (int i = obsIdx + 1; i < nSteps; i++) {
            int het = startHet + 2 * i;
            int hetPrev = het - 2;
            int homMinorPrev = (minorCount - hetPrev) / 2;
            int homMajorPrev = (majorCount - hetPrev) / 2;
            // P(het) / P(het-2) = 4 * homMinor(het-2) * homMajor(het-2) / (het * (het - 1))
            double ratio = 4.0 * homMinorPrev * homMajorPrev / ((double) het * (het - 1));
            probs[i] = probs[i - 1] * ratio;
        }

        // Walk downward from observed het
        for (int i = obsIdx - 1; i >= 0; i--) {
            int het = startHet + 2 * i;
            int hetNext = het + 2;
            int homMinor = (minorCount - het) / 2;
            int homMajor = (majorCount - het) / 2;
            // P(het) / P(het+2) = hetNext * (hetNext - 1) / (4 * homMinor * homMajor)
            double ratio = (double) hetNext * (hetNext - 1) / (4.0 * homMinor * homMajor);
            probs[i] = probs[i + 1] * ratio;
        }

        // Normalize to sum to 1
        double total = 0.0;
        for (double p : probs) total += p;
        if (total == 0.0) return 1.0;

        double obsProb = probs[obsIdx] / total;

        // Mid-p: sum of probs strictly less than observed + 0.5 * observed
        double pval = 0.0;
        for (int i = 0; i < nSteps; i++) {
            double p = probs[i] / total;
            if (p < obsProb) {
                pval += p;
            } else if (i == obsIdx) {
                pval += 0.5 * p;
            }
        }

        return pval;
    }
}
