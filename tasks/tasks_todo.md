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

### VectorOps SIMD Core (graphpop-procedures)
- [ ] Implement `VectorOps.dotProduct` — SIMD dot product via jdk.incubator.vector
- [ ] `VectorOps.cosineSimilarity` — cosine similarity between AF vectors
- [ ] `VectorOps.euclideanDistance` — Euclidean distance between AF vectors
- [ ] `VectorOps.alleleCounts` — sum allele counts across vector lanes
- [ ] Unit tests for all VectorOps methods (JUnit 5)
- [ ] Maven build verified: `mvn package` produces graphpop-procedures.jar

### Stored Procedure: `graphpop.diversity(chr, start, end, pop, [options])`
- [ ] Register as @Procedure, query Variant nodes in [start, end] range by chr
- [ ] Compute π (nucleotide diversity): Σ 2·p·(1-p) / L
- [ ] Compute θ_W (Watterson's theta): S / a_n
- [ ] Compute Tajima's D: (π - θ_W) / sqrt(Var)
- [ ] Compute H_e (expected heterozygosity): 2·p·(1-p) per pop
- [ ] Compute H_o (observed heterozygosity): het_count / (an/2) per pop
- [ ] Compute F_IS (inbreeding coefficient): 1 - H_o/H_e
- [ ] Support conditioned queries: filter by consequence, pathway, annotation
- [ ] Return results as Stream<Record> with all statistics

### Stored Procedure: `graphpop.sfs(chr, start, end, pop, [folded])`
- [ ] Compute site frequency spectrum: histogram of derived allele counts
- [ ] Support folded SFS (fold symmetrically)

### Stored Procedure: `graphpop.divergence(chr, start, end, pop1, pop2)`
- [ ] Compute Hudson's F_ST from allele count arrays
- [ ] Compute D_xy (net divergence): Σ (p1·(1-p2) + p2·(1-p1)) / L
- [ ] Compute D_a (absolute divergence)
- [ ] Return fst_hudson, dxy, da as Map

### Stored Procedure: `graphpop.joint_sfs(chr, start, end, pop1, pop2, [folded])`
- [ ] Compute 2D joint site frequency spectrum
- [ ] Dimensions: ac[pop1] × ac[pop2]

### Stored Procedure: `graphpop.genome_scan(chr, pop, window, step, [options])`
- [ ] Sliding window along NEXT chain (default: 100 kb window, 50 kb step)
- [ ] Materialize GenomicWindow nodes with computed stats (π, θ_W, D, H_e, H_o, F_IS)
- [ ] Link GenomicWindow → Chromosome via ON_CHROMOSOME
- [ ] Analysis versioning: each run tagged with run_id
- [ ] Comparative mode: include F_ST, D_xy if second population specified

### Index Strategy
- [ ] Composite index on (Variant.chr, Variant.pos) for range queries
- [ ] Full-text index on GOTerm.name for enrichment filtering

### Validation
- [ ] Reproduce vcftools/scikit-allel π for 1000G chr22 to <0.01% error
- [ ] Reproduce Tajima's D to <0.01% error
- [ ] Reproduce Hudson's F_ST to <0.01% error
- [ ] Benchmark: conditioned Tajima's D (missense only in pathway) vs classical pipeline

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
