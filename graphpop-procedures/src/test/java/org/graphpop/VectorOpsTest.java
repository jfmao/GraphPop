package org.graphpop;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.stream.Stream;

import static org.junit.jupiter.api.Assertions.*;

class VectorOpsTest {

    private static final double EPS = 1e-12;

    // -----------------------------------------------------------------------
    // dotProduct
    // -----------------------------------------------------------------------

    @Test
    void dotProductSimple() {
        double[] a = {1.0, 2.0, 3.0};
        double[] b = {4.0, 5.0, 6.0};
        assertEquals(32.0, VectorOps.dotProduct(a, b), EPS);
    }

    @Test
    void dotProductSingleElement() {
        assertEquals(6.0, VectorOps.dotProduct(new double[]{2.0}, new double[]{3.0}), EPS);
    }

    @Test
    void dotProductEmpty() {
        assertEquals(0.0, VectorOps.dotProduct(new double[]{}, new double[]{}), EPS);
    }

    @Test
    void dotProductLargeArray() {
        // Ensure SIMD loop + scalar tail both get exercised
        int n = 1023; // odd size to hit tail
        double[] a = new double[n];
        double[] b = new double[n];
        double expected = 0.0;
        for (int i = 0; i < n; i++) {
            a[i] = i + 1.0;
            b[i] = 1.0;
            expected += a[i];
        }
        assertEquals(expected, VectorOps.dotProduct(a, b), EPS);
    }

    // -----------------------------------------------------------------------
    // sum
    // -----------------------------------------------------------------------

    @Test
    void sumSimple() {
        assertEquals(10.0, VectorOps.sum(new double[]{1, 2, 3, 4}), EPS);
    }

    @Test
    void sumEmpty() {
        assertEquals(0.0, VectorOps.sum(new double[]{}), EPS);
    }

    @Test
    void sumLargeArray() {
        int n = 1000;
        double[] a = new double[n];
        for (int i = 0; i < n; i++) a[i] = 1.0;
        assertEquals(1000.0, VectorOps.sum(a), EPS);
    }

    // -----------------------------------------------------------------------
    // sumOfSquares
    // -----------------------------------------------------------------------

    @Test
    void sumOfSquaresSimple() {
        assertEquals(14.0, VectorOps.sumOfSquares(new double[]{1, 2, 3}), EPS);
    }

    // -----------------------------------------------------------------------
    // euclideanDistance
    // -----------------------------------------------------------------------

    @Test
    void euclideanDistanceSimple() {
        double[] a = {1.0, 0.0};
        double[] b = {0.0, 1.0};
        assertEquals(Math.sqrt(2.0), VectorOps.euclideanDistance(a, b), EPS);
    }

    @Test
    void euclideanDistanceSameVector() {
        double[] a = {3.0, 4.0, 5.0};
        assertEquals(0.0, VectorOps.euclideanDistance(a, a), EPS);
    }

    @Test
    void euclideanDistanceLargeArray() {
        int n = 500;
        double[] a = new double[n];
        double[] b = new double[n];
        for (int i = 0; i < n; i++) {
            a[i] = i;
            b[i] = i + 1;
        }
        // Each diff is 1.0, so distance = sqrt(n)
        assertEquals(Math.sqrt(n), VectorOps.euclideanDistance(a, b), EPS);
    }

    // -----------------------------------------------------------------------
    // cosineSimilarity
    // -----------------------------------------------------------------------

    @Test
    void cosineSimilarityParallel() {
        double[] a = {1.0, 2.0, 3.0};
        double[] b = {2.0, 4.0, 6.0};
        assertEquals(1.0, VectorOps.cosineSimilarity(a, b), EPS);
    }

    @Test
    void cosineSimilarityOrthogonal() {
        double[] a = {1.0, 0.0};
        double[] b = {0.0, 1.0};
        assertEquals(0.0, VectorOps.cosineSimilarity(a, b), EPS);
    }

    @Test
    void cosineSimilarityAntiparallel() {
        double[] a = {1.0, 2.0};
        double[] b = {-1.0, -2.0};
        assertEquals(-1.0, VectorOps.cosineSimilarity(a, b), EPS);
    }

    @Test
    void cosineSimilarityZeroVector() {
        double[] a = {0.0, 0.0};
        double[] b = {1.0, 2.0};
        assertEquals(0.0, VectorOps.cosineSimilarity(a, b), EPS);
    }

    // -----------------------------------------------------------------------
    // subtract
    // -----------------------------------------------------------------------

    @Test
    void subtractSimple() {
        double[] a = {5.0, 3.0, 1.0};
        double[] b = {1.0, 2.0, 3.0};
        double[] result = VectorOps.subtract(a, b);
        assertArrayEquals(new double[]{4.0, 1.0, -2.0}, result, EPS);
    }

    // -----------------------------------------------------------------------
    // expectedHeterozygosity
    // -----------------------------------------------------------------------

    @Test
    void expectedHeterozygosityFixed() {
        // Monomorphic: p=0 → He=0, p=1 → He=0
        double[] af = {0.0, 1.0};
        double[] he = VectorOps.expectedHeterozygosity(af);
        assertEquals(0.0, he[0], EPS);
        assertEquals(0.0, he[1], EPS);
    }

    @Test
    void expectedHeterozygosityMaximal() {
        // p=0.5 → He=0.5 (maximum)
        double[] af = {0.5};
        double[] he = VectorOps.expectedHeterozygosity(af);
        assertEquals(0.5, he[0], EPS);
    }

    @Test
    void expectedHeterozygosityGeneral() {
        // p=0.3 → He = 2*0.3*0.7 = 0.42
        double[] af = {0.3};
        double[] he = VectorOps.expectedHeterozygosity(af);
        assertEquals(0.42, he[0], EPS);
    }

    // -----------------------------------------------------------------------
    // piPerSite
    // -----------------------------------------------------------------------

    @Test
    void piPerSiteWithCorrection() {
        // p=0.5, an=100 → 2*0.5*0.5 * 100/99 = 0.50505...
        double pi = VectorOps.piPerSite(0.5, 100);
        assertEquals(2.0 * 0.5 * 0.5 * 100.0 / 99.0, pi, EPS);
    }

    @Test
    void piPerSiteTooFewAlleles() {
        assertEquals(0.0, VectorOps.piPerSite(0.5, 0), EPS);
        assertEquals(0.0, VectorOps.piPerSite(0.5, 1), EPS);
    }

    @Test
    void piPerSiteMonomorphic() {
        assertEquals(0.0, VectorOps.piPerSite(0.0, 100), EPS);
        assertEquals(0.0, VectorOps.piPerSite(1.0, 100), EPS);
    }

    // -----------------------------------------------------------------------
    // hudsonFstComponents
    // -----------------------------------------------------------------------

    @Test
    void hudsonFstIdenticalPops() {
        // Same frequency → Fst numerator should be near zero
        double[] fst = VectorOps.hudsonFstComponents(0.3, 100, 0.3, 100);
        // With same frequency, numerator ≈ -(hw1/(2*n1) + hw2/(2*n2)) which is slightly negative
        // But this is per-site; averaged over many sites it would be ~0
        assertTrue(Math.abs(fst[0]) < 0.01); // numerator very small
    }

    @Test
    void hudsonFstFixedDifference() {
        // p1=0, p2=1 → maximum differentiation
        double[] fst = VectorOps.hudsonFstComponents(0.0, 100, 1.0, 100);
        // Fst = numerator/denominator should be ~1
        assertTrue(fst[1] > 0); // denominator positive
        double fstVal = fst[0] / fst[1];
        assertTrue(fstVal > 0.95, "Fst for fixed difference should be ~1, got " + fstVal);
    }

    @Test
    void hudsonFstTooFewAlleles() {
        double[] fst = VectorOps.hudsonFstComponents(0.5, 1, 0.5, 100);
        assertEquals(0.0, fst[0], EPS);
        assertEquals(0.0, fst[1], EPS);
    }

    static Stream<Arguments> hudsonFstSymmetryArgs() {
        return Stream.of(
                Arguments.of(0.1, 100, 0.9, 200),
                Arguments.of(0.5, 50, 0.2, 150),
                Arguments.of(0.0, 100, 0.5, 100)
        );
    }

    @ParameterizedTest
    @MethodSource("hudsonFstSymmetryArgs")
    void hudsonFstSymmetric(double p1, int n1, double p2, int n2) {
        double[] ab = VectorOps.hudsonFstComponents(p1, n1, p2, n2);
        double[] ba = VectorOps.hudsonFstComponents(p2, n2, p1, n1);
        assertEquals(ab[0], ba[0], EPS, "Fst numerator should be symmetric");
        assertEquals(ab[1], ba[1], EPS, "Fst denominator should be symmetric");
    }

    // -----------------------------------------------------------------------
    // dxyPerSite
    // -----------------------------------------------------------------------

    @Test
    void dxyPerSiteSameFrequency() {
        // p1=p2 → Dxy = 2*p*(1-p), same as He
        double dxy = VectorOps.dxyPerSite(0.3, 0.3);
        assertEquals(2.0 * 0.3 * 0.7, dxy, EPS);
    }

    @Test
    void dxyPerSiteFixedDifference() {
        assertEquals(1.0, VectorOps.dxyPerSite(0.0, 1.0), EPS);
        assertEquals(1.0, VectorOps.dxyPerSite(1.0, 0.0), EPS);
    }

    @Test
    void dxyPerSiteSymmetric() {
        assertEquals(
                VectorOps.dxyPerSite(0.2, 0.8),
                VectorOps.dxyPerSite(0.8, 0.2),
                EPS
        );
    }

    // -----------------------------------------------------------------------
    // pearsonR2
    // -----------------------------------------------------------------------

    @Test
    void pearsonR2PerfectCorrelation() {
        // Identical dosage vectors → r² = 1.0
        double[] g = {0, 1, 2, 0, 1, 2, 1, 0};
        assertEquals(1.0, VectorOps.pearsonR2(g, g), EPS);
    }

    @Test
    void pearsonR2InverseCorrelation() {
        // g1 and g2 are perfectly inversely correlated: r = -1, r² = 1
        double[] g1 = {0, 0, 1, 1, 2, 2};
        double[] g2 = {2, 2, 1, 1, 0, 0};
        assertEquals(1.0, VectorOps.pearsonR2(g1, g2), EPS);
    }

    @Test
    void pearsonR2Uncorrelated() {
        // Construct orthogonal vectors after centering
        double[] g1 = {0, 2, 0, 2};
        double[] g2 = {0, 0, 2, 2};
        // mean1 = 1, mean2 = 1
        // centered: (-1,1,-1,1) and (-1,-1,1,1) → dot = 1-1-1+1 = 0 → r² = 0
        assertEquals(0.0, VectorOps.pearsonR2(g1, g2), EPS);
    }

    @Test
    void pearsonR2Monomorphic() {
        // One vector is constant → variance = 0 → r² = 0
        double[] g1 = {1, 1, 1, 1};
        double[] g2 = {0, 1, 2, 1};
        assertEquals(0.0, VectorOps.pearsonR2(g1, g2), EPS);
    }

    @Test
    void pearsonR2Empty() {
        assertEquals(0.0, VectorOps.pearsonR2(new double[]{}, new double[]{}), EPS);
    }

    @Test
    void pearsonR2SingleCarrier() {
        // Only one non-ref sample
        double[] g1 = {0, 0, 0, 1, 0};
        double[] g2 = {0, 0, 0, 1, 0};
        assertEquals(1.0, VectorOps.pearsonR2(g1, g2), EPS);
    }

    @Test
    void pearsonR2KnownValue() {
        // Manually computed: g1={0,1,2}, g2={0,0,2}
        // mean1=1, mean2=2/3
        // cov = (-1)(-2/3) + (0)(-2/3) + (1)(4/3) = 2/3 + 0 + 4/3 = 2
        // var1 = 1+0+1 = 2, var2 = 4/9+4/9+16/9 = 24/9 = 8/3
        // r = 2/sqrt(2*8/3) = 2/sqrt(16/3) = 2/(4/sqrt(3)) = sqrt(3)/2
        // r² = 3/4 = 0.75
        double[] g1 = {0, 1, 2};
        double[] g2 = {0, 0, 2};
        assertEquals(0.75, VectorOps.pearsonR2(g1, g2), 1e-10);
    }

    @Test
    void pearsonR2LargeArray() {
        // Ensure SIMD + scalar tail both work
        int n = 1023;
        double[] g1 = new double[n];
        double[] g2 = new double[n];
        for (int i = 0; i < n; i++) {
            g1[i] = i % 3; // 0,1,2,0,1,2,...
            g2[i] = i % 3;
        }
        assertEquals(1.0, VectorOps.pearsonR2(g1, g2), 1e-10);
    }

    // -----------------------------------------------------------------------
    // dPrime
    // -----------------------------------------------------------------------

    @Test
    void dPrimeCompleteLD() {
        // All haplotypes are either (0,0) or (1,1) → D' = 1.0
        int[] h1 = {0, 0, 1, 1, 0, 1};
        int[] h2 = {0, 0, 1, 1, 0, 1};
        assertEquals(1.0, VectorOps.dPrime(h1, h2, h1.length), EPS);
    }

    @Test
    void dPrimeNoLD() {
        // All four haplotype combinations in equal proportions → D = 0
        int[] h1 = {0, 0, 1, 1};
        int[] h2 = {0, 1, 0, 1};
        // pA=0.5, pB=0.5, pAB=0.25, D = 0.25-0.25 = 0
        assertEquals(0.0, VectorOps.dPrime(h1, h2, h1.length), EPS);
    }

    @Test
    void dPrimeMonomorphic() {
        int[] h1 = {0, 0, 0, 0};
        int[] h2 = {0, 1, 0, 1};
        assertEquals(0.0, VectorOps.dPrime(h1, h2, h1.length), EPS);
    }

    @Test
    void dPrimeEmpty() {
        assertEquals(0.0, VectorOps.dPrime(new int[]{}, new int[]{}, 0), EPS);
    }

    @Test
    void dPrimeIntermediate() {
        // 6 haplotypes: (0,0),(0,0),(0,1),(1,0),(1,1),(1,1)
        // pA = 3/6 = 0.5, pB = 3/6 = 0.5, pAB = 2/6 = 1/3
        // D = 1/3 - 0.25 = 1/12 > 0
        // Dmax = min(0.5*0.5, 0.5*0.5) = 0.25
        // D' = (1/12)/0.25 = 1/3
        int[] h1 = {0, 0, 0, 1, 1, 1};
        int[] h2 = {0, 0, 1, 0, 1, 1};
        assertEquals(1.0 / 3.0, VectorOps.dPrime(h1, h2, h1.length), 1e-10);
    }

    @Test
    void dPrimeNegativeD() {
        // Repulsion: alt at 1 tends to co-occur with ref at 2
        // (0,1),(0,1),(1,0),(1,0) → pA=0.5, pB=0.5, pAB=0 → D=-0.25
        // Dmax(D<0) = min(pA*pB, (1-pA)*(1-pB)) = min(0.25,0.25) = 0.25
        // |D'| = 0.25/0.25 = 1.0
        int[] h1 = {0, 0, 1, 1};
        int[] h2 = {1, 1, 0, 0};
        assertEquals(1.0, VectorOps.dPrime(h1, h2, h1.length), EPS);
    }
}
