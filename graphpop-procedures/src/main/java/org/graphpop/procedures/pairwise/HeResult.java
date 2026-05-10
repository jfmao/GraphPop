package org.graphpop.procedures.pairwise;

/**
 * Result row for {@code graphpop.kinship.branch_grm_he}: a single
 * heritability estimate for the supplied phenotype.
 */
public class HeResult {

    public double h2;
    public double se;
    public double num;          // y' · G · y
    public double tr_g_sq;      // tr(G²) (Hutchinson estimate)
    public long n_samples;
    public long n_hutchinson;
    public String method;       // "he_grm"

    public HeResult() {}

    public HeResult(double h2, double se, double num, double trGsq,
                    long nSamples, long nHutchinson, String method) {
        this.h2 = h2;
        this.se = se;
        this.num = num;
        this.tr_g_sq = trGsq;
        this.n_samples = nSamples;
        this.n_hutchinson = nHutchinson;
        this.method = method;
    }
}
