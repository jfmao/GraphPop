/**
 * Recombination-map inference (Phase 5, layer 5, M11).
 *
 * <p>Two complementary procedures: an ARG-derived per-window
 * breakpoint-density estimator
 * ({@link org.graphpop.procedures.recombination.ArgBreakpointsProcedure})
 * and an LD-based Hudson-Kaplan moment estimator
 * ({@link org.graphpop.procedures.recombination.LdDecayProcedure}).
 * The two are cross-validating: agreement signals a robust map,
 * disagreement signals model violations.</p>
 */
package org.graphpop.procedures.recombination;
