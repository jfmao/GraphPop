package org.graphpop.procedures;

/**
 * Fay and Wu's H statistic for detecting positive selection.
 *
 * <p>H = π − θ_H, where θ_H is based on the frequency of derived alleles.
 * A negative H indicates an excess of high-frequency derived alleles,
 * consistent with a selective sweep carrying a derived allele to high frequency.</p>
 *
 * <p>Requires ancestral/derived allele polarization. Variants without
 * ancestral allele annotation should be excluded from the computation.</p>
 *
 * <p>Reference: Fay JC, Wu CI. Hitchhiking under positive Darwinian selection.
 * Genetics. 2000;155(3):1405-13.</p>
 *
 * <p>Follows the same static utility class pattern as {@link TajimaD}.</p>
 */
final class FayWuH {

    private FayWuH() {}

    /**
     * Compute the per-site contribution to θ_H for a single variant.
     *
     * <p>θ_H per site = 2 * i² / (n * (n - 1)), where i is the derived allele count
     * and n is the number of sequences.</p>
     *
     * @param derivedCount number of derived alleles at this site
     * @param n            total number of sequences (= 2 * n_diploid)
     * @return per-site θ_H contribution
     */
    static double thetaHPerSite(int derivedCount, int n) {
        if (n < 2 || derivedCount <= 0 || derivedCount >= n) return 0.0;
        double i = derivedCount;
        return 2.0 * i * i / (n * (double) (n - 1));
    }

    /**
     * Compute Fay and Wu's H (unnormalized).
     *
     * <p>H = π − θ_H (both summed over sites, NOT divided by L).</p>
     *
     * @param piTotal     sum of per-site pi values over polarized sites
     * @param thetaHTotal sum of per-site θ_H values over polarized sites
     * @return H statistic
     */
    static double compute(double piTotal, double thetaHTotal) {
        return piTotal - thetaHTotal;
    }

    /**
     * Compute normalized Fay and Wu's H.
     *
     * <p>H_norm = H / sqrt(Var(H)), where the variance approximation follows
     * Zeng et al. (2006).</p>
     *
     * <p>Reference: Zeng K, Fu YX, Shi S, Wu CI. Statistical tests for detecting
     * positive selection by utilizing high-frequency variants.
     * Genetics. 2006;174(3):1431-9.</p>
     *
     * @param piTotal     sum of per-site pi values over polarized sites
     * @param thetaHTotal sum of per-site θ_H values over polarized sites
     * @param S           number of polarized segregating sites
     * @param n           number of sequences (= 2 * n_diploid)
     * @param a_n         harmonic number sum(1/i, i=1..n-1)
     * @return normalized H, or 0 if not computable
     */
    static double normalizedH(double piTotal, double thetaHTotal,
                              long S, int n, double a_n) {
        if (S == 0 || n < 4 || a_n == 0.0) return 0.0;

        double H = compute(piTotal, thetaHTotal);

        // Variance approximation from Zeng et al. (2006)
        // Var(H) ≈ theta * f(n) where theta = S / a_n
        double theta = S / a_n;

        double nD = n;
        // b_n = sum(1/i^2, i=1..n-1)
        double b_n = 0.0;
        for (int i = 1; i < n; i++) b_n += 1.0 / ((double) i * i);

        // Variance terms from Zeng et al. (2006) Eq. 9
        // thetaH coefficient:
        double nMinus1 = nD - 1.0;
        double term1 = theta * (nD - 2.0) / (6.0 * nMinus1);
        double term2 = theta * theta * (
                (18.0 * nD * nD * (3.0 * nD + 2.0) * b_n
                 - (88.0 * nD * nD * nD + 9.0 * nD * nD - 13.0 * nD + 6.0))
                / (9.0 * nD * nMinus1 * nMinus1)
        );

        double variance = term1 + term2;
        if (variance <= 0) return 0.0;

        return H / Math.sqrt(variance);
    }
}
