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
     * Compute Dxy (net divergence) for a single biallelic site.
     *
     * @param p1  allele frequency in population 1
     * @param p2  allele frequency in population 2
     * @return per-site Dxy = p1*(1-p2) + p2*(1-p1)
     */
    public static double dxyPerSite(double p1, double p2) {
        return p1 * (1.0 - p2) + p2 * (1.0 - p1);
    }
}
