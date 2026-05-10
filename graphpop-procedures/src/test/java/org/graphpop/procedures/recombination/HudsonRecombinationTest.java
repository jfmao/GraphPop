package org.graphpop.procedures.recombination;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class HudsonRecombinationTest {

    @Test
    void expected_r2_at_zero_matches_10_over_22() {
        assertEquals(10.0 / 22.0, HudsonRecombination.expectedR2(0.0), 1e-12);
        assertEquals(HudsonRecombination.R2_AT_ZERO,
                HudsonRecombination.expectedR2(0.0), 1e-12);
    }

    @Test
    void expected_r2_monotonically_decreasing() {
        double prev = HudsonRecombination.expectedR2(0.0);
        for (double c : new double[]{0.1, 1.0, 5.0, 50.0, 500.0, 5_000.0}) {
            double curr = HudsonRecombination.expectedR2(c);
            assertTrue(curr < prev,
                    "expectedR2(" + c + ") = " + curr
                        + " not < previous " + prev);
            prev = curr;
        }
    }

    @Test
    void expected_r2_negative_C_throws() {
        assertThrows(IllegalArgumentException.class,
                () -> HudsonRecombination.expectedR2(-0.1));
    }

    @Test
    void bisection_recovers_known_rho() {
        // Pick ρ = 1e-4 per bp, distance = 1000 bp → C = 0.1.
        double rho = 1e-4;
        double d = 1000.0;
        double c = rho * d;
        double r2 = HudsonRecombination.expectedR2(c);
        double recovered = HudsonRecombination.solveRhoFromMeanR2(r2, d);
        assertEquals(rho, recovered, 1e-7,
                "should recover ρ=1e-4 (C=0.1) from its E[r²]");
    }

    @Test
    void bisection_handles_high_recombination() {
        // ρ = 1e-2 per bp, distance = 100 bp → C = 1.0.
        double rho = 1e-2;
        double d = 100.0;
        double r2 = HudsonRecombination.expectedR2(rho * d);
        double recovered = HudsonRecombination.solveRhoFromMeanR2(r2, d);
        assertEquals(rho, recovered, 1e-7);
    }

    @Test
    void bisection_clamps_at_zero_when_r2_above_asymptote() {
        // Observed r² above the C=0 ceiling → ρ_hat = 0 (no signal).
        assertEquals(0.0, HudsonRecombination.solveRhoFromMeanR2(0.5, 1000.0));
        assertEquals(0.0, HudsonRecombination.solveRhoFromMeanR2(0.99, 1000.0));
    }

    @Test
    void bisection_returns_NaN_for_zero_r2() {
        // Numerically expectedR2(C→∞)→0; bisection clamps at rhoMax.
        // For an exact zero we return NaN (no valid information).
        assertTrue(Double.isNaN(
                HudsonRecombination.solveRhoFromMeanR2(0.0, 1000.0)));
    }

    @Test
    void bisection_validates_distance() {
        assertThrows(IllegalArgumentException.class,
                () -> HudsonRecombination.solveRhoFromMeanR2(0.1, 0.0));
        assertThrows(IllegalArgumentException.class,
                () -> HudsonRecombination.solveRhoFromMeanR2(0.1, -10.0));
    }
}
