package org.graphpop.procedures.pairwise;

import org.neo4j.graphdb.Result;
import org.neo4j.graphdb.Transaction;

import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Factory of {@link BranchWeightFn}s for the conditional predicates of
 * {@code graphpop.kinship.branch_grm} (M4.1 step 4).
 *
 * <p>Predicates:</p>
 * <ul>
 *   <li>{@link #pathwayWeight pathwayWeight} — branches whose mutations
 *       lie in a {@code :Pathway}. Closes G2.</li>
 *   <li>{@link #mutationFilterWeight mutationFilterWeight} — branches
 *       whose mutations have a given consequence type. Closes G2.</li>
 *   <li>{@link #timeWindowWeight timeWindowWeight} — geometric branch
 *       truncation to {@code [tLo, tHi]}.</li>
 *   <li>{@link #product product} — multiplicative composition.</li>
 * </ul>
 *
 * <p>{@link #fromOptions} reads the predicate keys from the procedure's
 * options Map and returns a single composed {@link BranchWeightFn}
 * (which is {@link BranchWeightFn#UNIT} when no predicate is set).</p>
 */
public final class BranchGrmConditioning {

    private BranchGrmConditioning() {}

    // ------------------------------------------------------------------
    // Time window
    // ------------------------------------------------------------------

    /**
     * Geometric branch-time truncation. For a branch of duration
     * {@code [childTime, parentTime]} and window {@code [tLo, tHi]},
     * returns the overlap fraction.
     *
     * @throws IllegalArgumentException if {@code tHi <= tLo}
     */
    public static BranchWeightFn timeWindowWeight(double tLo, double tHi) {
        if (tHi <= tLo) {
            throw new IllegalArgumentException(
                "time_window: tHi (" + tHi + ") must be > tLo (" + tLo + ")");
        }
        return (parentTskit, childTskit, parentTime, childTime, s, e) -> {
            double denom = parentTime - childTime;
            if (denom <= 0) return 0.0;
            double overlap = Math.min(parentTime, tHi) - Math.max(childTime, tLo);
            if (overlap <= 0) return 0.0;
            return overlap / denom;
        };
    }

    // ------------------------------------------------------------------
    // Pathway / mutation-class predicates ("lit-child" mask)
    // ------------------------------------------------------------------

    /**
     * Pathway-restricted weight: branch counts iff its child node carries
     * at least one {@code :MUTATED_ON} edge from a {@code :Variant} that
     * reaches the named pathway via
     * {@code :HAS_CONSEQUENCE -> :Gene -> :IN_PATHWAY -> :Pathway}.
     */
    public static BranchWeightFn pathwayWeight(Transaction tx, String runId,
                                               String pathwayId) {
        BitSet litChildren = loadLitChildren(tx, runId,
                "MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(:Gene)-[:IN_PATHWAY]->"
              + "(:Pathway {pathwayId: $pathway}) "
              + "MATCH (v)-[:MUTATED_ON {runId: $runId}]->(c:TreeNode) "
              + "RETURN DISTINCT c.nodeId AS nodeId",
                Map.of("runId", runId, "pathway", pathwayId));
        return litChildBranchWeight(litChildren);
    }

    /**
     * Consequence-restricted weight: branch counts iff its child node
     * carries at least one mutation whose
     * {@code :HAS_CONSEQUENCE.consequence} matches the supplied type.
     */
    public static BranchWeightFn mutationFilterWeight(Transaction tx,
                                                      String runId,
                                                      String consequence) {
        BitSet litChildren = loadLitChildren(tx, runId,
                "MATCH (v:Variant)-[c:HAS_CONSEQUENCE]->(:Gene) "
              + "WHERE c.consequence = $consequence "
              + "MATCH (v)-[:MUTATED_ON {runId: $runId}]->(tn:TreeNode) "
              + "RETURN DISTINCT tn.nodeId AS nodeId",
                Map.of("runId", runId, "consequence", consequence));
        return litChildBranchWeight(litChildren);
    }

    private static BitSet loadLitChildren(Transaction tx, String runId,
                                          String cypher,
                                          Map<String, Object> params) {
        BitSet bits = new BitSet();
        Result r = tx.execute(cypher, params);
        try {
            while (r.hasNext()) {
                Map<String, Object> row = r.next();
                int nid = ((Number) row.get("nodeId")).intValue();
                if (nid >= 0) bits.set(nid);
            }
        } finally {
            r.close();
        }
        return bits;
    }

    private static BranchWeightFn litChildBranchWeight(BitSet litChildren) {
        return (parentTskit, childTskit, pt, ct, s, e) ->
                litChildren.get(childTskit) ? 1.0 : 0.0;
    }

    // ------------------------------------------------------------------
    // Composition
    // ------------------------------------------------------------------

    /** Multiplicative composition of weight functions. */
    public static BranchWeightFn product(BranchWeightFn... fns) {
        if (fns.length == 0) return BranchWeightFn.UNIT;
        if (fns.length == 1) return fns[0];
        return (parentTskit, childTskit, pt, ct, s, e) -> {
            double w = 1.0;
            for (BranchWeightFn fn : fns) {
                w *= fn.weight(parentTskit, childTskit, pt, ct, s, e);
                if (w <= 0.0) return 0.0;  // short-circuit
            }
            return w;
        };
    }

    // ------------------------------------------------------------------
    // Options parsing
    // ------------------------------------------------------------------

    /**
     * Build the composed weight function from the procedure's options
     * Map. Recognised keys:
     *
     * <ul>
     *   <li>{@code restrict_to_pathway} (String) — pathwayId</li>
     *   <li>{@code mutation_filter} (String) — consequence type</li>
     *   <li>{@code time_window} (List of two Numbers, generations)</li>
     * </ul>
     */
    @SuppressWarnings("unchecked")
    public static BranchWeightFn fromOptions(Transaction tx, String runId,
                                             Map<String, Object> options) {
        if (options == null || options.isEmpty()) return BranchWeightFn.UNIT;

        java.util.List<BranchWeightFn> fns = new java.util.ArrayList<>();
        Object pathway = options.get("restrict_to_pathway");
        if (pathway instanceof String s && !s.isEmpty()) {
            fns.add(pathwayWeight(tx, runId, s));
        }
        Object consequence = options.get("mutation_filter");
        if (consequence instanceof String s && !s.isEmpty()) {
            fns.add(mutationFilterWeight(tx, runId, s));
        }
        Object timeWindow = options.get("time_window");
        if (timeWindow instanceof List<?> tw && tw.size() == 2) {
            double tLo = ((Number) tw.get(0)).doubleValue();
            double tHi = ((Number) tw.get(1)).doubleValue();
            fns.add(timeWindowWeight(tLo, tHi));
        }

        if (fns.isEmpty()) return BranchWeightFn.UNIT;
        return product(fns.toArray(new BranchWeightFn[0]));
    }
}
