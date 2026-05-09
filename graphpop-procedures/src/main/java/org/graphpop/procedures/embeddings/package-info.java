/**
 * Graph embedding procedures (Phase 4).
 *
 * <p>Thin wrappers over Neo4j Graph Data Science (GDS) algorithms applied to
 * Sample–Variant and Sample–Sample subgraphs: Leiden / Louvain community
 * detection, node2vec, FastRP, GraphSAGE.</p>
 *
 * <p>No SIMD here — heavy lifting belongs to GDS.</p>
 */
package org.graphpop.procedures.embeddings;
