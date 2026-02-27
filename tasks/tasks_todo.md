# GraphPop — Phase 1: Data Foundation + Summary Statistics

## Status: Not Started
## Target: Months 1–8

---

## Milestone 1.1 — Import Pipeline (Month 1–2)

### Setup
- [x] Initialize project structure (pyproject.toml, pom.xml, packages)
- [ ] Install and configure Neo4j 5.x in WSL2
- [ ] Verify Java 21+ with Vector API support
- [ ] Download benchmark datasets (1000 Genomes chr22, Rice 3K chr1 subset)

### VCF Parser (graphpop-import)
- [ ] Implement VCF reader using cyvcf2
- [ ] Compute per-population allele counts (ac, an, af, het_count, hom_alt_count)
- [ ] Handle missing genotypes correctly
- [ ] Emit Variant node CSV (with allele count arrays)
- [ ] Emit CARRIES relationship CSV (sparse: only non-ref genotypes)
- [ ] Emit Sample node CSV from VCF header + metadata file
- [ ] Emit Population node CSV
- [ ] Emit Chromosome node CSV

### Annotation Lift
- [ ] Parse VEP CSQ field → Gene nodes + HAS_CONSEQUENCE edges
- [ ] Pathway/GOTerm import from external files

### NEXT Edge Builder
- [ ] Sort variants by (chr, pos) → emit NEXT relationship CSV
- [ ] Include distance_bp; distance_cM from genetic map if available

### Bulk Import
- [ ] Write neo4j-admin import script with all CSVs
- [ ] Chunking: per-chromosome, 100K variant blocks
- [ ] Test: Import Rice 3K chr1 subset (100K variants) in <10 min
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
