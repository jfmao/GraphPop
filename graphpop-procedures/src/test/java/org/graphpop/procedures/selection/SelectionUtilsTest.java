package org.graphpop.procedures.selection;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class SelectionUtilsTest {

    @Test
    void welford_running_mean_and_sd() {
        SelectionUtils.Welford w = new SelectionUtils.Welford();
        for (double x : new double[]{1, 2, 3, 4, 5}) w.add(x);
        assertEquals(5L, w.n());
        assertEquals(3.0, w.mean(), 1e-12);
        // sample sd of {1..5} = sqrt(2.5)
        assertEquals(Math.sqrt(2.5), w.sd(), 1e-12);
    }

    @Test
    void welford_sd_is_NaN_below_two_samples() {
        SelectionUtils.Welford w = new SelectionUtils.Welford();
        assertTrue(Double.isNaN(w.sd()));
        w.add(42.0);
        assertTrue(Double.isNaN(w.sd()));
    }

    @Test
    void logit_freq_edges_endpoints_are_exact() {
        double[] edges = SelectionUtils.logitFreqEdges(5, 0.05, 0.95);
        assertEquals(0.05, edges[0], 1e-12);
        assertEquals(0.95, edges[edges.length - 1], 1e-12);
        assertEquals(6, edges.length);
        for (int i = 1; i < edges.length; i++) {
            assertTrue(edges[i] > edges[i - 1],
                    "edges must be strictly increasing");
        }
    }

    @Test
    void bin_index_partitions_correctly() {
        double[] edges = {0.05, 0.20, 0.50, 0.80, 0.95};
        // 4 bins: [0.05,0.20), [0.20,0.50), [0.50,0.80), [0.80,0.95)
        assertEquals(0, SelectionUtils.binIndex(edges, 0.10));
        assertEquals(1, SelectionUtils.binIndex(edges, 0.30));
        assertEquals(2, SelectionUtils.binIndex(edges, 0.60));
        assertEquals(3, SelectionUtils.binIndex(edges, 0.85));
        // Below first edge → bin 0
        assertEquals(0, SelectionUtils.binIndex(edges, 0.001));
        // At or above last edge → last bin
        assertEquals(3, SelectionUtils.binIndex(edges, 0.95));
        assertEquals(3, SelectionUtils.binIndex(edges, 0.99));
        // Edge values themselves go into the bin starting at that edge.
        assertEquals(0, SelectionUtils.binIndex(edges, 0.05));
        assertEquals(1, SelectionUtils.binIndex(edges, 0.20));
    }

    @Test
    void z_score_returns_zero_when_sd_is_undefined() {
        SelectionUtils.Welford w = new SelectionUtils.Welford();
        w.add(5.0);  // single sample → sd = NaN
        assertEquals(0.0, SelectionUtils.zScore(7.0, w));

        SelectionUtils.Welford zero = new SelectionUtils.Welford();
        zero.add(3.0);
        zero.add(3.0);  // sd = 0
        assertEquals(0.0, SelectionUtils.zScore(5.0, zero));
    }

    @Test
    void z_score_matches_classical_formula() {
        SelectionUtils.Welford w = new SelectionUtils.Welford();
        for (double x : new double[]{1, 2, 3, 4, 5}) w.add(x);
        // (5 - 3) / sqrt(2.5)
        assertEquals(2.0 / Math.sqrt(2.5),
                SelectionUtils.zScore(5.0, w), 1e-12);
    }

    @Test
    void sliding_windows_full_coverage() {
        long[][] wins = SelectionUtils.slidingWindows(100L, 30L, 30L);
        // [0,30), [30,60), [60,90), [90,100) -- last shorter, kept since
        // its length 10 >= window/2 = 15? No, 10 < 15, so dropped.
        assertEquals(3, wins.length);
        assertArrayEquals(new long[]{0, 30}, wins[0]);
        assertArrayEquals(new long[]{30, 60}, wins[1]);
        assertArrayEquals(new long[]{60, 90}, wins[2]);
    }

    @Test
    void sliding_windows_step_smaller_than_window() {
        long[][] wins = SelectionUtils.slidingWindows(100L, 50L, 25L);
        // Starts: 0, 25, 50, 75. Last: [75,100), length 25 >= 25. Kept.
        // Beyond start=100 stops.
        assertEquals(4, wins.length);
        assertArrayEquals(new long[]{0, 50}, wins[0]);
        assertArrayEquals(new long[]{75, 100}, wins[3]);
    }

    @Test
    void sliding_windows_zero_input_safe() {
        assertEquals(0, SelectionUtils.slidingWindows(0L, 10L, 10L).length);
        assertEquals(0, SelectionUtils.slidingWindows(100L, 0L, 10L).length);
        assertEquals(0, SelectionUtils.slidingWindows(100L, 10L, 0L).length);
    }
}
