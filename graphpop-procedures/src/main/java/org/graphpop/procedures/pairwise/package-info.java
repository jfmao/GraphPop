/**
 * Pairwise sample-level statistics (Phase 4).
 *
 * <p>Procedures in this package compute relationships between pairs of
 * samples by traversing the sparse {@code CARRIES} graph (full path).
 * Targets:</p>
 * <ul>
 *   <li>{@code graphpop.kinship} — KING-robust and GRM kinship coefficients.</li>
 *   <li>{@code graphpop.ibs} — identity-by-state matrices over windows.</li>
 *   <li>{@code graphpop.ibd} — IBD segment detection (hap-IBD-style).</li>
 * </ul>
 *
 * <p>Validation baselines: PLINK2 {@code --make-king}, hap-ibd, GERMLINE.</p>
 */
package org.graphpop.procedures.pairwise;
