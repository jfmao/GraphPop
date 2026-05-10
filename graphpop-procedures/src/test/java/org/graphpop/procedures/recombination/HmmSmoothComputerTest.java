package org.graphpop.procedures.recombination;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class HmmSmoothComputerTest {

    @Test
    void hmm_smooth_constant_input_returns_close_to_input() {
        // 10 windows all observing log10(ρ) = -7.0; smoothed values
        // should also be near -7.0.
        double[] obs = new double[10];
        for (int i = 0; i < obs.length; i++) obs[i] = -7.0;
        HmmSmoothComputer.Result r = HmmSmoothComputer.smooth(
                obs, 20, -10.0, -4.0, 0.5, 0.1);
        for (double s : r.smoothedLog10Rho) {
            assertEquals(-7.0, s, 0.4,
                    "smoothed " + s + " too far from -7.0");
        }
    }

    @Test
    void hmm_smooth_one_outlier_pulled_toward_neighbours() {
        // 10 windows at -7.0, one at -4.0; smoothing should pull the
        // outlier toward -7.0.
        double[] obs = new double[10];
        for (int i = 0; i < obs.length; i++) obs[i] = -7.0;
        obs[5] = -4.0;
        HmmSmoothComputer.Result r = HmmSmoothComputer.smooth(
                obs, 30, -10.0, -3.0, 0.3, 0.05);
        // Smoothed outlier should be closer to -7.0 than -4.0 was.
        assertTrue(r.smoothedLog10Rho[5] < -5.0,
                "outlier smoothed to " + r.smoothedLog10Rho[5]
                    + "; expected < -5.0");
        assertTrue(r.smoothedLog10Rho[5] > -7.5,
                "outlier smoothed too aggressively: "
                    + r.smoothedLog10Rho[5]);
    }

    @Test
    void hmm_smooth_nan_observations_handled_as_uniform_prior() {
        double[] obs = new double[6];
        obs[0] = -7.0;
        obs[1] = Double.NaN;
        obs[2] = Double.NaN;
        obs[3] = Double.NaN;
        obs[4] = -7.0;
        obs[5] = -7.0;
        HmmSmoothComputer.Result r = HmmSmoothComputer.smooth(
                obs, 20, -10.0, -4.0, 0.5, 0.1);
        // NaN windows should interpolate toward -7.0 via transition prior.
        for (int i = 1; i <= 3; i++) {
            assertTrue(r.smoothedLog10Rho[i] < -5.5,
                    "NaN window " + i + " = " + r.smoothedLog10Rho[i]);
            assertTrue(r.smoothedLog10Rho[i] > -8.5);
        }
    }

    @Test
    void hmm_smooth_invalid_n_states_throws() {
        assertThrows(IllegalArgumentException.class,
                () -> HmmSmoothComputer.smooth(
                        new double[]{-7.0}, 1, -10.0, -4.0, 0.5, 0.1));
    }

    @Test
    void hmm_smooth_invalid_bounds_throws() {
        assertThrows(IllegalArgumentException.class,
                () -> HmmSmoothComputer.smooth(
                        new double[]{-7.0}, 5, 0.0, -2.0, 0.5, 0.1));
    }

    @Test
    void hmm_smooth_invalid_emission_sd_throws() {
        assertThrows(IllegalArgumentException.class,
                () -> HmmSmoothComputer.smooth(
                        new double[]{-7.0}, 5, -10.0, -4.0, 0.0, 0.1));
    }

    @Test
    void hmm_smooth_invalid_switch_rate_throws() {
        assertThrows(IllegalArgumentException.class,
                () -> HmmSmoothComputer.smooth(
                        new double[]{-7.0}, 5, -10.0, -4.0, 0.5, 1.5));
        assertThrows(IllegalArgumentException.class,
                () -> HmmSmoothComputer.smooth(
                        new double[]{-7.0}, 5, -10.0, -4.0, 0.5, -0.1));
    }
}
