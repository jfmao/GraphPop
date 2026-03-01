package org.graphpop.procedures;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Unit tests for {@link FayWuH} (Fay and Wu's H statistic).
 */
class FayWuHTest {

    private static final double EPS = 1e-10;

    // ---- thetaHPerSite ----

    @Test
    void thetaHPerSite_knownValue() {
        // derivedCount=3, n=10: theta_H = 2 * 9 / (10 * 9) = 18/90 = 0.2
        assertEquals(0.2, FayWuH.thetaHPerSite(3, 10), EPS);
    }

    @Test
    void thetaHPerSite_singleDerived() {
        // derivedCount=1, n=10: theta_H = 2 * 1 / 90 = 0.02222...
        assertEquals(2.0 / 90.0, FayWuH.thetaHPerSite(1, 10), EPS);
    }

    @Test
    void thetaHPerSite_highDerived() {
        // derivedCount=9, n=10: theta_H = 2 * 81 / 90 = 1.8
        assertEquals(1.8, FayWuH.thetaHPerSite(9, 10), EPS);
    }

    @Test
    void thetaHPerSite_edgeCases() {
        // derivedCount=0 → 0
        assertEquals(0.0, FayWuH.thetaHPerSite(0, 10), EPS);
        // derivedCount=n → 0 (fixed, not segregating)
        assertEquals(0.0, FayWuH.thetaHPerSite(10, 10), EPS);
        // n < 2 → 0
        assertEquals(0.0, FayWuH.thetaHPerSite(1, 1), EPS);
    }

    // ---- compute (unnormalized H) ----

    @Test
    void compute_excessHighFreqDerived() {
        // If many high-frequency derived alleles → theta_H > pi → H < 0
        // Simulate: pi = 5.0, theta_H = 8.0 → H = 5 - 8 = -3
        assertEquals(-3.0, FayWuH.compute(5.0, 8.0), EPS);
    }

    @Test
    void compute_neutralExpectation() {
        // Under neutrality, pi ≈ theta_H → H ≈ 0
        assertEquals(0.0, FayWuH.compute(5.0, 5.0), EPS);
    }

    @Test
    void compute_deficitHighFreqDerived() {
        // If low-frequency derived alleles dominate → theta_H < pi → H > 0
        assertEquals(3.0, FayWuH.compute(8.0, 5.0), EPS);
    }

    // ---- normalizedH ----

    @Test
    void normalizedH_returnsZeroForNoData() {
        assertEquals(0.0, FayWuH.normalizedH(0, 0, 0, 10, 2.93), EPS);
        assertEquals(0.0, FayWuH.normalizedH(5, 8, 10, 3, 1.5), EPS); // n < 4
    }

    @Test
    void normalizedH_hasCorrectSign() {
        // For excess high-freq derived: H < 0, normalized should be < 0
        double piTotal = 20.0;
        double thetaHTotal = 40.0;
        long S = 100;
        int n = 100;
        double a_n = 0;
        for (int i = 1; i < n; i++) a_n += 1.0 / i;

        double hNorm = FayWuH.normalizedH(piTotal, thetaHTotal, S, n, a_n);
        assertTrue(hNorm < 0, "Normalized H should be negative for excess high-freq derived");
    }
}
