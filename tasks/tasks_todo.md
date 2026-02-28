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
- [x] Implement `VectorOps.dotProduct` — SIMD dot product via jdk.incubator.vector
- [x] `VectorOps.cosineSimilarity` — cosine similarity between AF vectors
- [x] `VectorOps.euclideanDistance` — Euclidean distance between AF vectors
- [x] `VectorOps.sum`, `sumOfSquares`, `subtract`, `expectedHeterozygosity`
- [x] `VectorOps.piPerSite`, `hudsonFstComponents`, `dxyPerSite`
- [x] Unit tests for all VectorOps methods (31 tests, JUnit 5)
- [x] Maven build verified: `mvnw package` produces graphpop-procedures.jar (22 KB)

### Stored Procedure: `graphpop.diversity(chr, start, end, pop, [options])`
- [x] Register as @Procedure, query Variant nodes in [start, end] range by chr
- [x] Compute π, θ_W, Tajima's D, H_e, H_o, F_IS
- [x] Support conditioned queries: filter by consequence, pathway
- [x] Return results as Stream<DiversityResult>

### Stored Procedure: `graphpop.sfs(chr, start, end, pop, [folded])`
- [x] Compute site frequency spectrum: histogram of derived allele counts
- [x] Support folded SFS (fold symmetrically)

### Stored Procedure: `graphpop.divergence(chr, start, end, pop1, pop2)`
- [x] Compute Hudson's F_ST, D_xy, D_a from allele count arrays

### Stored Procedure: `graphpop.joint_sfs(chr, start, end, pop1, pop2, [folded])`
- [x] Compute 2D joint site frequency spectrum

### Stored Procedure: `graphpop.genome_scan(chr, pop, window, step, [options])`
- [x] Sliding window (default: 100 kb window, 50 kb step)
- [x] Materialize GenomicWindow nodes with computed stats (π, θ_W, D, H_e, H_o, F_IS)
- [x] Link GenomicWindow → Chromosome via ON_CHROMOSOME
- [x] Analysis versioning: each run tagged with run_id
- [x] Comparative mode: include F_ST, D_xy, D_a if pop2 provided

### Deployment
- [x] Deploy script: `sudo bash scripts/deploy-procedures.sh`
- [x] Verify procedures registered in Neo4j (all 5 confirmed)
- [x] Create composite index on (Variant.chr, Variant.pos) for range queries
- [x] Full-text index on GOTerm.name for enrichment filtering
- [x] GenomicWindow constraint + indexes

### Validation
- [x] Reproduce scikit-allel π for 1000G chr22 to <0.01% error (0.000000% achieved)
- [x] Reproduce Tajima's D to <0.01% error (0.000001% achieved)
- [x] Reproduce Hudson's F_ST to <0.01% error (0.000000% achieved)
- [x] Benchmark: conditioned Tajima's D (missense in pathway) vs classical — 0.000000% error, 9x warm-cache speedup

### Milestone 1.2 — Complete ✓

---

## Milestone 1.3 — LD + Haplotype Statistics (Month 4–6)

### VectorOps LD Primitives
- [x] `VectorOps.pearsonR2` — SIMD-accelerated r² from dosage vectors
- [x] `VectorOps.dPrime` — D' from haplotype vectors (2×2 table)
- [x] Unit tests: perfect/inverse/zero correlation, monomorphic, known values, large arrays (14 new tests)

### GenotypeLoader Helper
- [x] `buildSampleIndex` — population-filtered sample→index mapping
- [x] `loadDosage` — single variant CARRIES→dosage array
- [x] `loadDosageBatch` — batch multi-variant dosage load (single Cypher query)
- [x] `loadHaplotypes` — phased haplotype loading (hap0/hap1 per sample)

### Stored Procedure: `graphpop.ld(chr, start, end, pop, max_dist, r2_threshold)`
- [x] Pairwise r² and D' computation within genomic window
- [x] Sub-window processing (200 variants) to bound memory
- [x] Batch dosage loading for efficiency
- [x] Allele frequency pre-filter (skip monomorphic sites)
- [x] Write LD edges: `(:Variant)-[:LD {r2, dprime, population}]->(:Variant)`

### EHHComputer Shared Helper
- [x] Bidirectional NEXT chain walk with haplotype group splitting
- [x] EHH computation: nCr(k_i,2) / nCr(n,2) per group
- [x] Trapezoidal integration (iHH) over physical distance
- [x] Configurable EHH cutoff (default 0.05)

### Stored Procedure: `graphpop.ihs(chr, pop, [options])`
- [x] Ancestral/derived haplotype partitioning from phased CARRIES data
- [x] iHH_ancestral and iHH_derived via EHHComputer
- [x] Unstandardized iHS = ln(iHH_A / iHH_D)
- [x] Frequency-bin standardization (20 bins, z-score)
- [x] Write ihs_{pop} and ihs_unstd_{pop} properties on Variant nodes
- [x] Configurable min_af filter (default 0.05)

### Stored Procedure: `graphpop.xpehh(chr, pop1, pop2, [options])`
- [x] Cross-population EHH comparison
- [x] Whole-population iHH computation per population
- [x] XP-EHH = ln(iHH_pop1 / iHH_pop2), genome-wide standardization
- [x] Write xpehh_{pop1}_{pop2} properties on Variant nodes

### Deployment
- [x] LD relationship index: `CREATE INDEX ld_pop FOR ()-[r:LD]-() ON (r.population)`
- [x] Smoke test for graphpop.ld in post-deploy-setup.sh
- [x] Build: `mvnw package` — 45 tests pass, JAR produced

### Validation
- [x] validate-ld.py: compare r² against scikit-allel rogers_huff_r
- [ ] Deploy and run smoke tests on live Neo4j instance
- [ ] Run validate-ld.py against live database

## Milestone 1.4 — Genome Scan + Benchmark (Month 6–7)
- [ ] ... (to be expanded)

## Milestone 1.5 — Scaling Validation (Month 7–8)
- [ ] ... (to be expanded)

---

## Review
(To be filled after milestone completion)
