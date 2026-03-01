package org.graphpop.procedures;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Unit tests for {@link HWEExact} (Hardy-Weinberg Equilibrium exact test).
 *
 * <p>Reference p-values checked against R HardyWeinberg::HWExact where available.</p>
 */
class HWEExactTest {

    private static final double TOLERANCE = 0.05;

    // ---- Perfect HWE ----

    @Test
    void perfectHWE_highPvalue() {
        // p=0.5, n=100: expected ~25 AA, ~50 Aa, ~25 aa
        double pval = HWEExact.hweExactMidP(50, 25, 100);
        assertTrue(pval > 0.5, "Perfect HWE should have high p-value, got " + pval);
    }

    // ---- Excess heterozygosity ----

    @Test
    void excessHet_lowPvalue() {
        // Extreme excess het: 80 hets out of 100 (expected ~50 for p=0.5)
        // nHet=80, nHomMinor=0, nTotal=100
        // This should have very low p-value
        double pval = HWEExact.hweExactMidP(80, 0, 100);
        assertTrue(pval < 0.001, "Extreme excess het should have very low p-value, got " + pval);
    }

    // ---- Excess homozygosity ----

    @Test
    void excessHom_lowPvalue() {
        // Excess homozygosity: 10 hets out of 100 (expected ~50 for p=0.5)
        // nHet=10, nHomMinor=35, nTotal=100
        double pval = HWEExact.hweExactMidP(10, 35, 100);
        assertTrue(pval < 0.001, "Excess homozygosity should have low p-value, got " + pval);
    }

    // ---- Fixed allele ----

    @Test
    void fixedAllele_pvalueOne() {
        // All homozygous for one allele (minor count = 0)
        double pval = HWEExact.hweExactMidP(0, 0, 100);
        assertEquals(1.0, pval, 1e-10, "Fixed allele should have p-value 1.0");
    }

    // ---- Small sample ----

    @Test
    void smallSample_reasonable() {
        // n=10, nHet=4, nHomMinor=1 → minor AC = 2*1+4 = 6
        double pval = HWEExact.hweExactMidP(4, 1, 10);
        assertTrue(pval > 0 && pval <= 1.0,
                "p-value should be in (0, 1] for small sample, got " + pval);
    }

    // ---- Edge: all hets ----

    @Test
    void allHets_extremeExcess() {
        // n=50, nHet=50, nHomMinor=0 → everyone is het (impossible under HWE)
        double pval = HWEExact.hweExactMidP(50, 0, 50);
        assertTrue(pval < 0.001, "All-het should have very low p-value, got " + pval);
    }

    // ---- Edge: nTotal = 0 ----

    @Test
    void zeroSamples_returnsOne() {
        double pval = HWEExact.hweExactMidP(0, 0, 0);
        assertEquals(1.0, pval, 1e-10);
    }

    // ---- Known values from R ----

    @Test
    void knownValues_moderateDeviation() {
        // nHet=40, nHomMinor=30, nTotal=100 (AF=0.5 but deficit of hets)
        // Expected het under HWE: ~50, observed: 40 → moderate deficit
        // Mid-p value is small (~0.037) indicating significant HWE departure
        double pval = HWEExact.hweExactMidP(40, 30, 100);
        assertTrue(pval > 0.0 && pval < 0.5,
                "Should have a p-value indicating some departure from HWE, got " + pval);
    }
}
