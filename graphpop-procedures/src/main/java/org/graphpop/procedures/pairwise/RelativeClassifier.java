package org.graphpop.procedures.pairwise;

/**
 * Pairwise relationship classifier from KING-robust kinship + IBS0
 * fraction (Manichaikul et al. 2010 cutoffs).
 *
 * <pre>
 *   identical (MZ twin / dup) : phi >= 0.354
 *   parent_child              : 0.177 <= phi < 0.354 AND ibs0_frac < 0.0050
 *   full_sibling              : 0.177 <= phi < 0.354 AND ibs0_frac >= 0.0050
 *   second_degree             : 0.0884 <= phi < 0.177
 *   third_degree              : 0.0442 <= phi < 0.0884
 *   unrelated                 : phi < 0.0442
 * </pre>
 *
 * <p>Cutoffs may be overridden via {@link Cutoffs}. Pure function;
 * no Neo4j dependency.</p>
 */
public final class RelativeClassifier {

    private RelativeClassifier() {}

    /** Threshold knobs; defaults match Manichaikul 2010. */
    public static final class Cutoffs {
        public final double identical;       // phi >= identical -> "identical"
        public final double firstDegree;     // phi >= firstDegree (and < identical) -> first-degree
        public final double secondDegree;    // phi >= secondDegree (and < firstDegree) -> 2nd
        public final double thirdDegree;     // phi >= thirdDegree (and < secondDegree) -> 3rd
        public final double ibs0Threshold;   // ibs0_frac < this -> parent-child (else full-sib)

        public Cutoffs(double identical, double firstDegree,
                       double secondDegree, double thirdDegree,
                       double ibs0Threshold) {
            this.identical = identical;
            this.firstDegree = firstDegree;
            this.secondDegree = secondDegree;
            this.thirdDegree = thirdDegree;
            this.ibs0Threshold = ibs0Threshold;
        }

        public static Cutoffs defaults() {
            return new Cutoffs(0.354, 0.177, 0.0884, 0.0442, 0.0050);
        }
    }

    public static final class Verdict {
        public final String relationship;
        public final int degree;       // 0=identical, 1=PC/FS, 2, 3, 99=unrelated

        Verdict(String relationship, int degree) {
            this.relationship = relationship;
            this.degree = degree;
        }
    }

    /**
     * Classify one pair from {@code phi} + {@code ibs0_frac}. Pass
     * {@link Double#NaN} for {@code ibs0_frac} when no IBS0 signal is
     * available (e.g. IBD source); the classifier returns
     * {@code "first_degree"} as the umbrella label rather than
     * disambiguating parent-child vs full-sib.
     */
    public static Verdict classify(double phi, double ibs0Frac, Cutoffs c) {
        if (phi >= c.identical) return new Verdict("identical", 0);
        if (phi >= c.firstDegree) {
            if (Double.isNaN(ibs0Frac)) {
                return new Verdict("first_degree", 1);
            }
            if (ibs0Frac < c.ibs0Threshold) {
                return new Verdict("parent_child", 1);
            }
            return new Verdict("full_sibling", 1);
        }
        if (phi >= c.secondDegree) return new Verdict("second_degree", 2);
        if (phi >= c.thirdDegree) return new Verdict("third_degree", 3);
        return new Verdict("unrelated", 99);
    }

    public static Verdict classify(double phi, double ibs0Frac) {
        return classify(phi, ibs0Frac, Cutoffs.defaults());
    }
}
