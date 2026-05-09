package org.graphpop.procedures.pairwise;

/**
 * Per-branch weight multiplier in {@code [0, 1]} used to apply a
 * conditional predicate to the branch GRM accumulator.
 *
 * <p>{@link BranchGrmComputer} multiplies a branch's contribution to
 * {@code mu} (and therefore to both the matrix update and
 * {@code total_mu}) by this weight. {@code 0} skips the branch
 * entirely; {@code 1} is the unconditional case
 * (see {@link #UNIT}).</p>
 *
 * <p>Implementations may inspect the branch's tskit ids, times, and
 * the active marginal-tree interval. They must be pure / referentially
 * transparent — {@link BranchGrmComputer#compute} may call them in any
 * order and any number of times.</p>
 */
@FunctionalInterface
public interface BranchWeightFn {

    /**
     * @param parentTskitId tskit local id of the edge's parent node
     * @param childTskitId  tskit local id of the edge's child node
     * @param parentTime    parent time (generations)
     * @param childTime     child time (generations)
     * @param intervalStart genomic interval start (bp, inclusive)
     * @param intervalEnd   genomic interval end (bp, exclusive)
     * @return weight in {@code [0, 1]}
     */
    double weight(int parentTskitId, int childTskitId,
                  double parentTime, double childTime,
                  long intervalStart, long intervalEnd);

    /** Unit weight: every branch counts at full strength. */
    BranchWeightFn UNIT = (p, c, pt, ct, s, e) -> 1.0;
}
