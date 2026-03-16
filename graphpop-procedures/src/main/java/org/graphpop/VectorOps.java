package org.graphpop;

import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.VectorSpecies;

/**
 * SIMD-accelerated numeric operations for population genomics.
 *
 * <p>Uses the Java Vector API ({@code jdk.incubator.vector}) to perform
 * vectorised arithmetic on allele-frequency arrays and genotype vectors.</p>
 */
public final class VectorOps {

    private static final VectorSpecies<Double> SPECIES = DoubleVector.SPECIES_PREFERRED;

    private VectorOps() {
        // utility class
    }

    /**
     * Compute the dot product of two double arrays using SIMD lanes.
     *
     * @param a first array
     * @param b second array (must be same length as {@code a})
     * @return dot product
     */
    public static double dotProduct(double[] a, double[] b) {
        int i = 0;
        int upperBound = SPECIES.loopBound(a.length);
        DoubleVector sum = DoubleVector.zero(SPECIES);

        for (; i < upperBound; i += SPECIES.length()) {
            DoubleVector va = DoubleVector.fromArray(SPECIES, a, i);
            DoubleVector vb = DoubleVector.fromArray(SPECIES, b, i);
            sum = va.fma(vb, sum);
        }

        double result = sum.reduceLanes(jdk.incubator.vector.VectorOperators.ADD);

        // Scalar tail
        for (; i < a.length; i++) {
            result += a[i] * b[i];
        }
        return result;
    }

    /**
     * Compute the sum of all elements in an array using SIMD lanes.
     *
     * @param a the array
     * @return sum of elements
     */
    public static double sum(double[] a) {
        int i = 0;
        int upperBound = SPECIES.loopBound(a.length);
        DoubleVector acc = DoubleVector.zero(SPECIES);

        for (; i < upperBound; i += SPECIES.length()) {
            DoubleVector va = DoubleVector.fromArray(SPECIES, a, i);
            acc = acc.add(va);
        }

        double result = acc.reduceLanes(jdk.incubator.vector.VectorOperators.ADD);

        for (; i < a.length; i++) {
            result += a[i];
        }
        return result;
    }

    /**
     * Compute the sum of squares of all elements using SIMD lanes.
     *
     * @param a the array
     * @return sum of squared elements
     */
    public static double sumOfSquares(double[] a) {
        return dotProduct(a, a);
    }

    /**
     * Compute the Euclidean distance between two vectors.
     *
     * @param a first vector
     * @param b second vector (same length as {@code a})
     * @return Euclidean distance
     */
    public static double euclideanDistance(double[] a, double[] b) {
        int i = 0;
        int upperBound = SPECIES.loopBound(a.length);
        DoubleVector acc = DoubleVector.zero(SPECIES);

        for (; i < upperBound; i += SPECIES.length()) {
            DoubleVector va = DoubleVector.fromArray(SPECIES, a, i);
            DoubleVector vb = DoubleVector.fromArray(SPECIES, b, i);
            DoubleVector diff = va.sub(vb);
            acc = diff.fma(diff, acc);
        }

        double result = acc.reduceLanes(jdk.incubator.vector.VectorOperators.ADD);

        for (; i < a.length; i++) {
            double d = a[i] - b[i];
            result += d * d;
        }
        return Math.sqrt(result);
    }

    /**
     * Compute cosine similarity between two vectors.
     *
     * @param a first vector
     * @param b second vector (same length as {@code a})
     * @return cosine similarity in [-1, 1]
     */
    public static double cosineSimilarity(double[] a, double[] b) {
        int i = 0;
        int upperBound = SPECIES.loopBound(a.length);
        DoubleVector dotAcc = DoubleVector.zero(SPECIES);
        DoubleVector normAAcc = DoubleVector.zero(SPECIES);
        DoubleVector normBAcc = DoubleVector.zero(SPECIES);

        for (; i < upperBound; i += SPECIES.length()) {
            DoubleVector va = DoubleVector.fromArray(SPECIES, a, i);
            DoubleVector vb = DoubleVector.fromArray(SPECIES, b, i);
            dotAcc = va.fma(vb, dotAcc);
            normAAcc = va.fma(va, normAAcc);
            normBAcc = vb.fma(vb, normBAcc);
        }

        double dot = dotAcc.reduceLanes(jdk.incubator.vector.VectorOperators.ADD);
        double normA = normAAcc.reduceLanes(jdk.incubator.vector.VectorOperators.ADD);
        double normB = normBAcc.reduceLanes(jdk.incubator.vector.VectorOperators.ADD);

        for (; i < a.length; i++) {
            dot += a[i] * b[i];
            normA += a[i] * a[i];
            normB += b[i] * b[i];
        }

        double denom = Math.sqrt(normA) * Math.sqrt(normB);
        return denom == 0.0 ? 0.0 : dot / denom;
    }

    /**
     * Element-wise subtraction: result[i] = a[i] - b[i].
     *
     * @param a first array
     * @param b second array (same length)
     * @return new array with differences
     */
    public static double[] subtract(double[] a, double[] b) {
        double[] result = new double[a.length];
        int i = 0;
        int upperBound = SPECIES.loopBound(a.length);

        for (; i < upperBound; i += SPECIES.length()) {
            DoubleVector va = DoubleVector.fromArray(SPECIES, a, i);
            DoubleVector vb = DoubleVector.fromArray(SPECIES, b, i);
            va.sub(vb).intoArray(result, i);
        }

        for (; i < a.length; i++) {
            result[i] = a[i] - b[i];
        }
        return result;
    }

    /**
     * Compute expected heterozygosity for each population: He = 2*p*(1-p).
     *
     * @param af allele frequency array (one per population)
     * @return array of He values
     */
    public static double[] expectedHeterozygosity(double[] af) {
        double[] he = new double[af.length];
        int i = 0;
        int upperBound = SPECIES.loopBound(af.length);
        DoubleVector two = DoubleVector.broadcast(SPECIES, 2.0);
        DoubleVector one = DoubleVector.broadcast(SPECIES, 1.0);

        for (; i < upperBound; i += SPECIES.length()) {
            DoubleVector p = DoubleVector.fromArray(SPECIES, af, i);
            // He = 2 * p * (1 - p)
            DoubleVector q = one.sub(p);
            two.mul(p).mul(q).intoArray(he, i);
        }

        for (; i < af.length; i++) {
            he[i] = 2.0 * af[i] * (1.0 - af[i]);
        }
        return he;
    }

    /**
     * Compute per-site nucleotide diversity contribution: 2*p*(1-p) * n/(n-1).
     * This applies the finite-sample correction factor.
     *
     * @param af allele frequency
     * @param an allele number (2 * number of diploid samples with data)
     * @return per-site pi contribution (0 if an &lt; 2)
     */
    public static double piPerSite(double af, int an) {
        if (an < 2) return 0.0;
        return 2.0 * af * (1.0 - af) * an / (an - 1.0);
    }

    /**
     * Compute Hudson's Fst for a single biallelic site between two populations.
     *
     * <p>Hudson's estimator: Fst = 1 - (Hw / Hb)
     * where Hw = average within-population heterozygosity,
     * Hb = between-population heterozygosity.</p>
     *
     * @param p1  allele frequency in population 1
     * @param n1  allele number in population 1
     * @param p2  allele frequency in population 2
     * @param n2  allele number in population 2
     * @return per-site Fst numerator and denominator as a double[2]
     */
    public static double[] hudsonFstComponents(double p1, int n1, double p2, int n2) {
        if (n1 < 2 || n2 < 2) return new double[]{0.0, 0.0};

        double hw1 = 2.0 * p1 * (1.0 - p1) * n1 / (n1 - 1.0);
        double hw2 = 2.0 * p2 * (1.0 - p2) * n2 / (n2 - 1.0);
        double hw = (hw1 + hw2) / 2.0;

        double hb = (p1 - p2) * (p1 - p2) - hw1 / (2.0 * n1) - hw2 / (2.0 * n2);

        double numerator = hb;
        double denominator = hb + hw;

        return new double[]{numerator, denominator};
    }

    /**
     * Compute Pearson's correlation coefficient squared (r²) between two
     * genotype dosage vectors.
     *
     * <p>r² = [cov(X,Y)]² / [var(X) · var(Y)]
     * where dosage is 0 (ref/ref), 1 (het), 2 (hom_alt).</p>
     *
     * @param geno1 dosage vector for variant 1 (length = n_samples)
     * @param geno2 dosage vector for variant 2 (same length)
     * @return r² in [0, 1], or 0 if either variant is monomorphic
     */
    public static double pearsonR2(double[] geno1, double[] geno2) {
        int n = geno1.length;
        if (n == 0) return 0.0;

        // Compute sums using SIMD
        double sum1 = sum(geno1);
        double sum2 = sum(geno2);
        double mean1 = sum1 / n;
        double mean2 = sum2 / n;

        // Compute centered dot product and variances in one SIMD pass
        int i = 0;
        int upperBound = SPECIES.loopBound(n);
        DoubleVector covAcc = DoubleVector.zero(SPECIES);
        DoubleVector var1Acc = DoubleVector.zero(SPECIES);
        DoubleVector var2Acc = DoubleVector.zero(SPECIES);
        DoubleVector vMean1 = DoubleVector.broadcast(SPECIES, mean1);
        DoubleVector vMean2 = DoubleVector.broadcast(SPECIES, mean2);

        for (; i < upperBound; i += SPECIES.length()) {
            DoubleVector v1 = DoubleVector.fromArray(SPECIES, geno1, i).sub(vMean1);
            DoubleVector v2 = DoubleVector.fromArray(SPECIES, geno2, i).sub(vMean2);
            covAcc = v1.fma(v2, covAcc);
            var1Acc = v1.fma(v1, var1Acc);
            var2Acc = v2.fma(v2, var2Acc);
        }

        double cov = covAcc.reduceLanes(jdk.incubator.vector.VectorOperators.ADD);
        double var1 = var1Acc.reduceLanes(jdk.incubator.vector.VectorOperators.ADD);
        double var2 = var2Acc.reduceLanes(jdk.incubator.vector.VectorOperators.ADD);

        // Scalar tail
        for (; i < n; i++) {
            double d1 = geno1[i] - mean1;
            double d2 = geno2[i] - mean2;
            cov += d1 * d2;
            var1 += d1 * d1;
            var2 += d2 * d2;
        }

        double denom = var1 * var2;
        if (denom == 0.0) return 0.0;

        double r = cov / Math.sqrt(denom);
        return r * r;
    }

    /**
     * Compute D' (normalized linkage disequilibrium coefficient) from
     * haplotype vectors.
     *
     * <p>Each haplotype vector contains 0 (ref) or 1 (alt) for each haplotype.
     * D' = D / D_max, where D = p_AB - p_A * p_B.</p>
     *
     * @param hap1 haplotype vector for variant 1 (length = n_haplotypes)
     * @param hap2 haplotype vector for variant 2 (same length)
     * @param n    number of haplotypes
     * @return |D'| in [0, 1], or 0 if either variant is monomorphic
     */
    public static double dPrime(int[] hap1, int[] hap2, int n) {
        if (n == 0) return 0.0;

        // Build 2x2 haplotype frequency table
        int n11 = 0; // both alt
        int n10 = 0; // alt at 1, ref at 2
        int n01 = 0; // ref at 1, alt at 2
        for (int i = 0; i < n; i++) {
            if (hap1[i] == 1) {
                if (hap2[i] == 1) n11++;
                else n10++;
            } else {
                if (hap2[i] == 1) n01++;
            }
        }
        double pA = (double) (n11 + n10) / n; // freq of alt at variant 1
        double pB = (double) (n11 + n01) / n; // freq of alt at variant 2
        double pAB = (double) n11 / n;

        // Check monomorphic
        if (pA == 0.0 || pA == 1.0 || pB == 0.0 || pB == 1.0) return 0.0;

        double D = pAB - pA * pB;

        // D_max depends on sign of D
        double Dmax;
        if (D >= 0) {
            Dmax = Math.min(pA * (1.0 - pB), (1.0 - pA) * pB);
        } else {
            Dmax = Math.min(pA * pB, (1.0 - pA) * (1.0 - pB));
        }

        if (Dmax == 0.0) return 0.0;
        return Math.abs(D) / Dmax;
    }

    /**
     * Compute D' from byte[] haplotype vectors directly (avoids int[] copy).
     *
     * <p>Identical algorithm to {@link #dPrime(int[], int[], int)} but uses
     * branchless counting on byte values (0 or 1).</p>
     *
     * @param hap1 haplotype vector for variant 1 (byte values 0 or 1)
     * @param hap2 haplotype vector for variant 2 (same length)
     * @param n    number of haplotypes
     * @return |D'| in [0, 1], or 0 if either variant is monomorphic
     */
    public static double dPrime(byte[] hap1, byte[] hap2, int n) {
        if (n == 0) return 0.0;

        int n11 = 0, n10 = 0, n01 = 0;
        for (int i = 0; i < n; i++) {
            int a = hap1[i];
            int b = hap2[i];
            n11 += a & b;
            n10 += a & (b ^ 1);
            n01 += (a ^ 1) & b;
        }

        double pA = (double) (n11 + n10) / n;
        double pB = (double) (n11 + n01) / n;
        double pAB = (double) n11 / n;

        if (pA == 0.0 || pA == 1.0 || pB == 0.0 || pB == 1.0) return 0.0;

        double D = pAB - pA * pB;

        double Dmax;
        if (D >= 0) {
            Dmax = Math.min(pA * (1.0 - pB), (1.0 - pA) * pB);
        } else {
            Dmax = Math.min(pA * pB, (1.0 - pA) * (1.0 - pB));
        }

        if (Dmax == 0.0) return 0.0;
        return Math.abs(D) / Dmax;
    }

    /**
     * Compute D' from bit-packed haplotype byte arrays.
     *
     * <p>Each byte stores 8 haplotypes (bit 0 = haplotype 0, bit 7 = haplotype 7).
     * Uses {@link Integer#bitCount} for efficient bulk counting.</p>
     *
     * @param packed1     bit-packed haplotype row for variant 1
     * @param packed2     bit-packed haplotype row for variant 2
     * @param nHaplotypes total number of haplotypes
     * @return |D'| in [0, 1], or 0 if either variant is monomorphic
     */
    public static double dPrimePacked(byte[] packed1, byte[] packed2, int nHaplotypes) {
        if (nHaplotypes == 0) return 0.0;

        int fullBytes = nHaplotypes >> 3;
        int remainBits = nHaplotypes & 7;

        int n11 = 0, nA = 0, nB = 0;
        for (int i = 0; i < fullBytes; i++) {
            int a = packed1[i] & 0xFF;
            int b = packed2[i] & 0xFF;
            n11 += Integer.bitCount(a & b);
            nA += Integer.bitCount(a);
            nB += Integer.bitCount(b);
        }
        if (remainBits > 0) {
            int mask = (1 << remainBits) - 1;
            int a = packed1[fullBytes] & mask;
            int b = packed2[fullBytes] & mask;
            n11 += Integer.bitCount(a & b);
            nA += Integer.bitCount(a);
            nB += Integer.bitCount(b);
        }

        double pA = (double) nA / nHaplotypes;
        double pB = (double) nB / nHaplotypes;
        double pAB = (double) n11 / nHaplotypes;

        if (pA == 0.0 || pA == 1.0 || pB == 0.0 || pB == 1.0) return 0.0;

        double D = pAB - pA * pB;

        double Dmax;
        if (D >= 0) {
            Dmax = Math.min(pA * (1.0 - pB), (1.0 - pA) * pB);
        } else {
            Dmax = Math.min(pA * pB, (1.0 - pA) * (1.0 - pB));
        }

        if (Dmax == 0.0) return 0.0;
        return Math.abs(D) / Dmax;
    }

    /**
     * Compute Dxy (net divergence) for a single biallelic site.
     *
     * @param p1  allele frequency in population 1
     * @param p2  allele frequency in population 2
     * @return per-site Dxy = p1*(1-p2) + p2*(1-p1)
     */
    public static double dxyPerSite(double p1, double p2) {
        return p1 * (1.0 - p2) + p2 * (1.0 - p1);
    }

    /**
     * Compute Weir &amp; Cockerham (1984) Fst variance components for a single
     * biallelic locus between two populations.
     *
     * <p>Returns {a, b, c} where:
     * <ul>
     *   <li>a = between-population variance component</li>
     *   <li>b = within-population, between-individual variance component</li>
     *   <li>c = within-individual (gametic) variance component</li>
     * </ul>
     * Fst = a / (a + b + c), computed as ratio of sums across loci.</p>
     *
     * <p>Uses diploid sample counts (n_i = number of diploid individuals = an_i / 2).
     * Formulas follow Weir &amp; Cockerham (1984) equations 2, 3, and 4 for r=2
     * populations, matching scikit-allel's {@code weir_cockerham_fst}.</p>
     *
     * @param ac1  allele count in population 1
     * @param an1  allele number in population 1 (2 × diploid count)
     * @param het1 heterozygote count in population 1
     * @param ac2  allele count in population 2
     * @param an2  allele number in population 2
     * @param het2 heterozygote count in population 2
     * @return double[3] = {a, b, c}; all zeros if either population has fewer than 2 alleles
     */
    public static double[] wcFstComponents(int ac1, int an1, int het1,
                                           int ac2, int an2, int het2) {
        if (an1 < 2 || an2 < 2) return new double[]{0.0, 0.0, 0.0};

        int r = 2; // number of populations
        double n1 = an1 / 2.0; // diploid sample count
        double n2 = an2 / 2.0;
        double nTotal = n1 + n2;
        double nBar = nTotal / r;

        // Guard: nBar <= 1 causes division by zero in 1/(nBar-1) below.
        // W&C Fst is undefined with fewer than 2 diploid individuals per pop on average.
        if (nBar <= 1.0) return new double[]{0.0, 0.0, 0.0};

        // n_c: sample size correction factor (Eq. in W&C 1984 section 3)
        double nC = nTotal - (n1 * n1 + n2 * n2) / nTotal;

        double p1 = (double) ac1 / an1;
        double p2 = (double) ac2 / an2;

        // Weighted average allele frequency
        double pBar = (n1 * p1 + n2 * p2) / nTotal;

        // Sample variance of allele frequencies
        double s2 = (n1 * (p1 - pBar) * (p1 - pBar)
                    + n2 * (p2 - pBar) * (p2 - pBar)) / ((r - 1) * nBar);

        // Average observed heterozygote frequency
        double hBar = (het1 / n1 + het2 / n2) / r;

        // Variance components (W&C 1984 Eqs. 2-4)
        double a = nBar / nC * (s2 - 1.0 / (nBar - 1.0)
                * (pBar * (1.0 - pBar) - (r - 1.0) / r * s2 - hBar / 4.0));

        double b = nBar / (nBar - 1.0)
                * (pBar * (1.0 - pBar) - (r - 1.0) / r * s2
                   - (2.0 * nBar - 1.0) / (4.0 * nBar) * hBar);

        double c = hBar / 2.0;

        return new double[]{a, b, c};
    }
}
