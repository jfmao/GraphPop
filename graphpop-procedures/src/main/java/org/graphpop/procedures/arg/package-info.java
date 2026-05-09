/**
 * Ancestral recombination graph (ARG) procedures (Phase 4).
 *
 * <p>ARGs are stored natively in Neo4j as {@code (:TreeNode)} and
 * {@code (:TreeNode)-[:PARENT_OF {start, end}]->(:TreeNode)} edges, with
 * {@code (:TreeNode)-[:REPRESENTS]->(:Sample)} for sample leaves and
 * {@code (:Variant)-[:LIES_ON]->(:ARGEdge)} for mutation placement.
 * See "Phase 4 — ARG layer" in
 * {@code docs/GraphPop_Compiled_Holistic_Design.md}.</p>
 *
 * <p>Targets:</p>
 * <ul>
 *   <li>{@code graphpop.arg.tmrca} — pairwise time-to-most-recent-common-ancestor.</li>
 *   <li>{@code graphpop.arg.coalescence_rate} — coalescence rate over time bins.</li>
 *   <li>{@code graphpop.arg.branch_diversity} — branch-length-based diversity.</li>
 * </ul>
 *
 * <p>Validation baselines: tskit branch-length statistics on
 * msprime-simulated ground truth; tsinfer / Relate / SINGER on empirical data.</p>
 */
package org.graphpop.procedures.arg;
