# GraphPop — Phase 1: Data Foundation + Summary Statistics

## Status: In Progress
## Target: Months 1–8

---

## Milestone 1.1 — Import Pipeline (Month 1–2)

### Setup
- [x] Initialize project structure (pyproject.toml, pom.xml, packages)
- [~] Install and configure Neo4j 5.x in WSL2 (neo4j 2026.01.4 installed; run `sudo bash scripts/configure-neo4j.sh` to apply config)
- [x] Verify Java 21+ with Vector API support
- [x] Download benchmark datasets (1000 Genomes NYGC 30x chr22 — 496 MB; Rice 3K deferred, SNP-Seek offline)

### VCF Parser (graphpop-import)
- [x] Implement VCF reader using cyvcf2
- [x] Compute per-population allele counts (ac, an, af, het_count, hom_alt_count)
- [x] Handle missing genotypes correctly
- [x] Emit Variant node CSV (with allele count arrays)
- [x] Emit CARRIES relationship CSV (sparse: only non-ref genotypes)
- [x] Emit Sample node CSV from VCF header + metadata file
- [x] Emit Population node CSV
- [x] Emit Chromosome node CSV

### Annotation Lift
- [x] Parse VEP CSQ field → Gene nodes + HAS_CONSEQUENCE edges
- [x] Pathway/GOTerm import from external files

### NEXT Edge Builder
- [x] Sort variants by (chr, pos) → emit NEXT relationship CSV
- [~] Include distance_cM from genetic map if available (deferred — requires genetic map file)

### Bulk Import
- [x] Write neo4j-admin import script with all CSVs
- [x] Chunking: per-chromosome, 100K variant blocks
- [x] Test: Import 1000G chr22 into Neo4j and validate
- [x] Post-import validation: orphaned nodes, NEXT chain completeness, AF consistency

### Milestone 1.1 — Complete ✓

---

## Milestone 1.2 — VectorOps + Core Statistics (Month 2–4)

### Neo4j Stored Procedures (graphpop-procedures)
- [ ] Set up Maven build with Neo4j 2026.x procedure API
- [ ] `graphpop.vectorAF(variantId)` — return allele-frequency vector across populations
- [ ] `graphpop.vectorHet(variantId)` — return heterozygosity vector
- [ ] `graphpop.cosineSimilarity(v1, v2)` — cosine similarity between two AF vectors
- [ ] `graphpop.euclideanDistance(v1, v2)` — Euclidean distance between AF vectors
- [ ] Use Java Vector API (jdk.incubator.vector) for SIMD-accelerated array ops
- [ ] Unit tests for all procedures (JUnit 5 + Neo4j harness)

### Summary Statistics (per-variant)
- [ ] Expected heterozygosity: He = 1 − Σ(pi²) per population
- [ ] Observed heterozygosity: Ho = het_count / n per population
- [ ] Nucleotide diversity (π) via NEXT-edge window traversal
- [ ] Store computed stats as Variant node properties or separate :Statistic nodes

### Population Differentiation
- [ ] Fst (Weir & Cockerham 1984) — pairwise between all population pairs
- [ ] Store Fst as DIFFERENTIATES relationship or Population node property
- [ ] Fst confidence intervals via jackknife over NEXT-linked variant blocks

### Windowed Statistics
- [ ] Sliding-window π using NEXT chain traversal (Cypher procedure)
- [ ] Sliding-window Fst across NEXT chain
- [ ] Window size configurable (default: 50 kb, step: 10 kb)

---

## Milestone 1.3 — LD + Haplotype Statistics (Month 4–6)
- [ ] ... (to be expanded)

## Milestone 1.4 — Genome Scan + Benchmark (Month 6–7)
- [ ] ... (to be expanded)

## Milestone 1.5 — Scaling Validation (Month 7–8)
- [ ] ... (to be expanded)

---

## Review
(To be filled after milestone completion)
