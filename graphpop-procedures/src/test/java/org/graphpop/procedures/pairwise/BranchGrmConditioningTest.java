package org.graphpop.procedures.pairwise;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Pure-Java tests for non-Cypher pieces of {@link BranchGrmConditioning}:
 * {@link BranchGrmConditioning#timeWindowWeight}, {@link BranchGrmConditioning#product},
 * and {@link BranchWeightFn#UNIT}. The Cypher-driven pathway and mutation
 * filters are exercised by {@code BranchGrmProcedureTest}.
 */
class BranchGrmConditioningTest {

    private static final double EPS = 1e-12;

    private static double w(BranchWeightFn fn,
                            double parentTime, double childTime) {
        return fn.weight(0, 1, parentTime, childTime, 0L, 100L);
    }

    @Test
    void timeWindow_windowContainsBranch_returnsOne() {
        BranchWeightFn fn = BranchGrmConditioning.timeWindowWeight(0.0, 10.0);
        // Branch [1, 5]; window [0, 10] strictly contains it.
        assertEquals(1.0, w(fn, 5.0, 1.0), EPS);
    }

    @Test
    void timeWindow_windowSplitsBranch_returnsFraction() {
        BranchWeightFn fn = BranchGrmConditioning.timeWindowWeight(2.0, 4.0);
        // Branch [1, 5]; window [2, 4] => overlap = 4-2 = 2; denom = 4 => 0.5.
        assertEquals(0.5, w(fn, 5.0, 1.0), EPS);
    }

    @Test
    void timeWindow_windowMissesBranch_returnsZero() {
        BranchWeightFn fn = BranchGrmConditioning.timeWindowWeight(10.0, 20.0);
        // Branch [1, 5] entirely below window.
        assertEquals(0.0, w(fn, 5.0, 1.0), EPS);
    }

    @Test
    void timeWindow_windowTouchesBranchOnLeft_returnsFraction() {
        BranchWeightFn fn = BranchGrmConditioning.timeWindowWeight(0.0, 3.0);
        // Branch [1, 5]; window [0, 3]; overlap = 3-1 = 2; denom = 4 => 0.5.
        assertEquals(0.5, w(fn, 5.0, 1.0), EPS);
    }

    @Test
    void timeWindow_zeroLengthBranch_returnsZero() {
        BranchWeightFn fn = BranchGrmConditioning.timeWindowWeight(0.0, 10.0);
        assertEquals(0.0, w(fn, 5.0, 5.0), EPS);
    }

    @Test
    void timeWindow_negativeRange_throws() {
        assertThrows(IllegalArgumentException.class,
            () -> BranchGrmConditioning.timeWindowWeight(5.0, 3.0));
    }

    @Test
    void product_isMultiplicative() {
        BranchWeightFn half = (a, b, c, d, e, f) -> 0.5;
        BranchWeightFn third = (a, b, c, d, e, f) -> 1.0 / 3.0;
        BranchWeightFn composed = BranchGrmConditioning.product(half, third);
        assertEquals(1.0 / 6.0, w(composed, 1.0, 0.0), EPS);
    }

    @Test
    void product_emptyReturnsUnit() {
        BranchWeightFn fn = BranchGrmConditioning.product();
        assertEquals(1.0, w(fn, 1.0, 0.0), EPS);
    }

    @Test
    void product_shortCircuitsOnZero() {
        // The second function would return NaN, but the first returns 0
        // and product short-circuits.
        BranchWeightFn zero = (a, b, c, d, e, f) -> 0.0;
        BranchWeightFn nan = (a, b, c, d, e, f) -> Double.NaN;
        BranchWeightFn composed = BranchGrmConditioning.product(zero, nan);
        assertEquals(0.0, w(composed, 1.0, 0.0), EPS);
    }

    @Test
    void unit_alwaysOne() {
        assertEquals(1.0, w(BranchWeightFn.UNIT, 1.0, 0.0), EPS);
        assertEquals(1.0, w(BranchWeightFn.UNIT, 0.0, 0.0), EPS);
    }
}
