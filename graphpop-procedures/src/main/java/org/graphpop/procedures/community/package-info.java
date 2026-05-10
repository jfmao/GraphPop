/**
 * Community-detection procedures (Phase 4 layer 4 / deferred from M5).
 *
 * <p>Pure-Java Louvain modularity optimisation
 * ({@link org.graphpop.procedures.community.LouvainProcedure},
 * Blondel et al. 2008). No external GDS dependency — keeps the
 * deployment self-contained and the test harness uniform with the
 * rest of {@code graphpop-procedures}.</p>
 */
package org.graphpop.procedures.community;
