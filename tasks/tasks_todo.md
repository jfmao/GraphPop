# GraphPop — Phase 1: Data Foundation + Summary Statistics

## Status: Not Started
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
- [ ] Parse VEP CSQ field → Gene nodes + HAS_CONSEQUENCE edges
- [ ] Pathway/GOTerm import from external files

### NEXT Edge Builder
- [x] Sort variants by (chr, pos) → emit NEXT relationship CSV
- [ ] Include distance_cM from genetic map if available

### Bulk Import
- [x] Write neo4j-admin import script with all CSVs
- [x] Chunking: per-chromosome, 100K variant blocks
- [ ] Test: Import 1000G chr22 into Neo4j and validate
- [ ] Post-import validation: orphaned nodes, NEXT chain completeness, AF consistency

---

## Milestone 1.2 — VectorOps + Core Statistics (Month 2–4)
- [ ] ... (to be expanded when 1.1 is complete)

## Milestone 1.3 — LD + Haplotype Statistics (Month 4–6)
- [ ] ... (to be expanded)

## Milestone 1.4 — Genome Scan + Benchmark (Month 6–7)
- [ ] ... (to be expanded)

## Milestone 1.5 — Scaling Validation (Month 7–8)
- [ ] ... (to be expanded)

---

## Review
(To be filled after milestone completion)
