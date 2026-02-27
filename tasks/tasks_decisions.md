# GraphPop — Architectural Decision Log

## Format
Each decision: date, context, options considered, decision, rationale.

---

## Decisions

### ADR-001: Graph Schema v3 — Stratified Aggregation Model
**Date:** Pre-project (design phase)
**Context:** Storing genotype data for 1M+ samples
**Options:** (A) Dosage vectors on Variant nodes, (B) Sparse CARRIES relationships
**Decision:** B — Sparse CARRIES + allele count arrays on Variant nodes
**Rationale:** Dosage vectors → 1 PB at All of Us scale. Sparse CARRIES → ~135 TB. Fast-path stats use allele count arrays (O(V×K), independent of N).

### ADR-002: Neo4j as sole storage engine
**Date:** Pre-project (design phase)
**Context:** Whether to use external columnar store (Parquet, HDF5) alongside Neo4j
**Decision:** All computation graph-native. No external store.
**Rationale:** Core thesis of the project. If we need external stores, the thesis fails.

### ADR-003: Java Vector API for SIMD
**Date:** Pre-project (design phase)
**Context:** Implementing VectorOps for numerical computation
**Decision:** Java 21+ Vector API (incubator module)
**Rationale:** Native to Neo4j plugin ecosystem. Avoids JNI overhead of C++. Benchmarkable. Fallback to C++ only if Java performance is insufficient.
