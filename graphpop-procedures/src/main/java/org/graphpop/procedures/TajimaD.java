package org.graphpop.procedures;

/**
 * Tajima's D computation following Tajima (1989).
 *
 * <p>Expects the number of sequences (n = 2 * n_diploid_samples) and uses
 * the pre-computed harmonic numbers a_n = H(n-1) and a_n2 = H2(n-1).</p>
 */
final class TajimaD {

    private TajimaD() {}

    /**
     * Compute Tajima's D statistic.
     *
     * @param piTotal sum of per-site pi values (NOT divided by L)
     * @param S       number of segregating sites
     * @param n       number of sequences (= 2 * n_diploid for diploids)
     * @param a1      harmonic number sum(1/i, i=1..n-1)
     * @param a2      sum(1/i^2, i=1..n-1)
     * @return Tajima's D, or 0 if not computable
     */
    static double compute(double piTotal, long S, int n, double a1, double a2) {
        if (S == 0 || a1 == 0 || n < 4) return 0.0;

        double thetaW = S / a1;
        double d = piTotal - thetaW;

        double b1 = (n + 1.0) / (3.0 * (n - 1.0));
        double b2 = 2.0 * (n * (double) n + n + 3.0) / (9.0 * n * (n - 1.0));

        double c1 = b1 - 1.0 / a1;
        double c2 = b2 - (n + 2.0) / (a1 * n) + a2 / (a1 * a1);

        double e1 = c1 / a1;
        double e2 = c2 / (a1 * a1 + a2);

        double variance = e1 * S + e2 * S * (S - 1);
        if (variance <= 0) return 0.0;

        return d / Math.sqrt(variance);
    }

    /**
     * Compute Tajima's D using a PopulationContext.
     * Derives n from the stored n_samples (n = 2 * n_samples).
     */
    static double compute(double piTotal, long S, PopulationContext ctx) {
        int n = 2 * ctx.nSamples; // number of sequences
        return compute(piTotal, S, n, ctx.a_n, ctx.a_n2);
    }
}
