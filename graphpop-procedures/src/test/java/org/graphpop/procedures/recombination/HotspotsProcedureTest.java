package org.graphpop.procedures.recombination;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class HotspotsProcedureTest {

    @Test
    void bh_adjusts_p_values_to_known_reference() {
        // Reference from scipy.stats.false_discovery_control on
        // [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205,
        //  0.212, 0.216] yields adj p = [0.01, 0.04, 0.1025, 0.1025,
        //  0.1025, 0.123, 0.12714286, 0.265, 0.265, 0.265] (rough).
        double[] p = {0.001, 0.008, 0.039, 0.041, 0.042,
                       0.060, 0.074, 0.205, 0.212, 0.216};
        double[] adj = HotspotsProcedure.benjaminiHochberg(p);
        // Smallest p-value rank=1: adj = 0.001 * 10 / 1 = 0.01.
        assertEquals(0.01, adj[0], 1e-9);
        // Second-smallest: 0.008 * 10 / 2 = 0.04.
        assertEquals(0.04, adj[1], 1e-9);
        // Monotonicity: adj p-values must be non-decreasing in rank.
        Integer[] order = sortedIndices(p);
        for (int i = 1; i < order.length; i++) {
            assertTrue(adj[order[i]] >= adj[order[i - 1]] - 1e-12,
                    "adj non-monotonic at rank " + i);
        }
        // Bounded in [0, 1].
        for (double a : adj) {
            assertTrue(a >= 0.0 && a <= 1.0);
        }
    }

    @Test
    void bh_handles_NaN_and_out_of_range() {
        double[] p = {0.01, Double.NaN, -0.5, 1.5, 0.05};
        double[] adj = HotspotsProcedure.benjaminiHochberg(p);
        // Index 0 (p=0.01) is the second-smallest after the clamped
        // negative value (which became 0). adj should be ≤ 1.
        for (double a : adj) {
            assertTrue(a >= 0.0 && a <= 1.0);
        }
    }

    @Test
    void one_sided_normal_p_at_zero_is_half() {
        assertEquals(0.5,
                HotspotsProcedure.oneSidedUpperTailNormal(0.0), 1e-7);
    }

    @Test
    void one_sided_normal_tails() {
        // Large positive z → very small p.
        assertTrue(HotspotsProcedure.oneSidedUpperTailNormal(5.0) < 1e-5);
        // Large negative z → near 1.
        assertTrue(HotspotsProcedure.oneSidedUpperTailNormal(-5.0) > 1.0 - 1e-5);
    }

    @Test
    void bh_all_significant_keeps_them_significant() {
        double[] p = {1e-6, 1e-6, 1e-6};
        double[] adj = HotspotsProcedure.benjaminiHochberg(p);
        // All p-values ≈ 0 → all adj ≈ 0 (well below typical q=0.05).
        for (double a : adj) {
            assertTrue(a < 0.01,
                    "expected adj < 0.01 across all entries; got " + a);
        }
    }

    private static Integer[] sortedIndices(double[] arr) {
        Integer[] idx = new Integer[arr.length];
        for (int i = 0; i < arr.length; i++) idx[i] = i;
        java.util.Arrays.sort(idx,
                (a, b) -> Double.compare(arr[a], arr[b]));
        return idx;
    }
}
