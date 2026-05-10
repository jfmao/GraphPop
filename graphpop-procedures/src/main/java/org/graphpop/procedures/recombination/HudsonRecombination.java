package org.graphpop.procedures.recombination;

/**
 * Hudson 1985 closed-form expected r² between two segregating
 * sites as a function of the population-scaled recombination
 * distance C = ρ · d:
 *
 * <pre>
 * E[r² | C] = (10 + C) / ((2 + C) · (11 + 13·C + C²))
 * </pre>
 *
 * <p>and a bisection routine that solves for ρ given an observed
 * mean r² and a physical distance d. The expectation is
 * monotonically decreasing in C with limits</p>
 *
 * <ul>
 *   <li>{@code C → 0}     → 10/22 ≈ 0.4545</li>
 *   <li>{@code C → ∞}     → 0</li>
 * </ul>
 *
 * <p>so {@link #solveRhoFromMeanR2} returns 0 when the observed
 * mean r² lies above the {@code C = 0} bound (small-sample
 * sampling noise dominates).</p>
 */
public final class HudsonRecombination {

    private HudsonRecombination() {}

    /** Asymptotic value at {@code C → 0}: {@code 10/22}. */
    public static final double R2_AT_ZERO = 10.0 / 22.0;

    /**
     * Hudson 1985 expected r² for a population-scaled recombination
     * distance {@code C = ρ · d}.
     */
    public static double expectedR2(double c) {
        if (c < 0.0) {
            throw new IllegalArgumentException("C must be ≥ 0; got " + c);
        }
        double num = 10.0 + c;
        double denom = (2.0 + c) * (11.0 + 13.0 * c + c * c);
        return num / denom;
    }

    /**
     * Solve for ρ (per-bp recombination probability) given an observed
     * mean r² and physical distance d. Returns 0 when the observed
     * r² exceeds the {@code C = 0} bound.
     *
     * @param meanR2  observed mean r² (must be in (0, 1])
     * @param distance physical distance d in bp (must be > 0)
     * @param rhoMin  lower bracket for ρ-per-bp (default 1e-12)
     * @param rhoMax  upper bracket for ρ-per-bp (default 1.0)
     * @param tol     relative tolerance for bisection
     * @param maxIter cap on bisection iterations
     */
    public static double solveRhoFromMeanR2(double meanR2, double distance,
                                              double rhoMin, double rhoMax,
                                              double tol, int maxIter) {
        if (!(meanR2 > 0.0)) {
            return Double.NaN;
        }
        if (!(distance > 0.0)) {
            throw new IllegalArgumentException("distance must be > 0");
        }
        if (meanR2 >= R2_AT_ZERO) {
            // Observed r² above the asymptotic bound — no signal.
            return 0.0;
        }
        // Convert bracket to C-space and bisect there.
        double lo = rhoMin * distance;
        double hi = rhoMax * distance;
        double fLo = expectedR2(lo) - meanR2;
        double fHi = expectedR2(hi) - meanR2;
        if (fLo < 0.0) {
            // Even C=lo gives r² below observation → ρ even smaller.
            return rhoMin;
        }
        if (fHi > 0.0) {
            // Even C=hi gives r² above observation → ρ larger; return ceiling.
            return rhoMax;
        }
        for (int i = 0; i < maxIter; i++) {
            double mid = 0.5 * (lo + hi);
            double fMid = expectedR2(mid) - meanR2;
            if (fMid > 0.0) lo = mid;
            else hi = mid;
            if (hi - lo < tol * Math.max(hi, 1e-12)) break;
        }
        return 0.5 * (lo + hi) / distance;
    }

    /** Convenience overload with sensible defaults. */
    public static double solveRhoFromMeanR2(double meanR2, double distance) {
        return solveRhoFromMeanR2(meanR2, distance,
                1e-12, 1.0, 1e-9, 200);
    }
}
