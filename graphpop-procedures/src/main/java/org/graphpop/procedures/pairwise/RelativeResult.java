package org.graphpop.procedures.pairwise;

/**
 * Result row for {@code graphpop.relate.classify}: per-pair
 * categorical relationship.
 */
public class RelativeResult {

    public String sample_a;
    public String sample_b;
    public String relationship;
    public long degree;
    public double phi;
    public double ibs0_frac;
    public String source;

    public RelativeResult() {}

    public RelativeResult(String a, String b, String relationship,
                          long degree, double phi, double ibs0Frac,
                          String source) {
        this.sample_a = a;
        this.sample_b = b;
        this.relationship = relationship;
        this.degree = degree;
        this.phi = phi;
        this.ibs0_frac = ibs0Frac;
        this.source = source;
    }
}
