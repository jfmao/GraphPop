/**
 * Selection-inference procedures (Phase 5, layer 5).
 *
 * <p>ARG-aware selection scans: per-variant allele-age z-scores
 * conditional on derived-allele frequency
 * ({@link org.graphpop.procedures.selection.AlleleAgeScanProcedure}),
 * and per-window branch-length outlier scans
 * ({@link org.graphpop.procedures.selection.BranchOutlierScanProcedure}).
 * Both build on the M6 ARG-statistics primitives.</p>
 */
package org.graphpop.procedures.selection;
