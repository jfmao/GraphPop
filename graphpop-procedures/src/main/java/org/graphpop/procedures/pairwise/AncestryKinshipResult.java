package org.graphpop.procedures.pairwise;

/**
 * Per-(sample, sample, ancestry) row from
 * {@code graphpop.kinship.branch_grm_by_ancestry}. The
 * {@code b_ij_component} is the share of the unconditional
 * {@code B_ij} attributable to branches painted as this ancestry; a
 * sum across ancestries reproduces the unconditional matrix exactly
 * (modulo unpainted nodes).
 */
public class AncestryKinshipResult {

    public String sample_a;
    public String sample_b;
    public String ancestry;
    public double b_ij_component;
    public long n_branches;
    public String method;

    public AncestryKinshipResult() {}

    public AncestryKinshipResult(String a, String b, String ancestry,
                                  double bij, long nBranches, String method) {
        this.sample_a = a;
        this.sample_b = b;
        this.ancestry = ancestry;
        this.b_ij_component = bij;
        this.n_branches = nBranches;
        this.method = method;
    }
}
