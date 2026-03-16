# GraphPop — Task Tracker

## Status: Memory Management Complete, Ready for Re-run
## Updated: 2026-03-16

---

## Phase 1 — Data Foundation + Summary Statistics: COMPLETE

### Milestone 1.1 — Import Pipeline: COMPLETE
- [x] VCF parser (cyvcf2), per-population allele counts
- [x] CSV emitter (Variant, Sample, Population, Chromosome, NEXT)
- [x] VEP annotation lift (Gene + HAS_CONSEQUENCE)
- [x] Pathway/GOTerm import (Reactome + BioMart)
- [x] neo4j-admin bulk import script + validation

### Milestone 1.2 — VectorOps + Core Statistics: COMPLETE
- [x] VectorOps SIMD core (11 methods, 31 tests)
- [x] 5 stored procedures: diversity, divergence, sfs, joint_sfs, genome_scan
- [x] Validated against scikit-allel (<0.000001% error)
- [x] Conditioned query benchmark: 9x warm-cache speedup

### Milestone 1.3 — LD + Haplotype Statistics: COMPLETE
- [x] VectorOps: pearsonR2, dPrime (14 new tests)
- [x] 3 new procedures: graphpop.ld, graphpop.ihs, graphpop.xpehh (8 total)
- [x] HaplotypeMatrix: dense bit-packed matrix from packed arrays
- [x] EHHComputer: position-sorted walks, trapezoidal iHH
- [x] LD validated: r=0.996 vs scikit-allel

### Milestone 1.4 — Performance Optimization & Benchmark: COMPLETE
- [x] EHHComputer: flat array tracking + incremental sumPairs
- [x] VectorOps: dPrime(byte[]) branchless overload + dPrimePacked (bitCount-based)
- [x] GenomeScanProcedure: pre-computed per-variant stat arrays
- [x] Benchmark suite: scripts/benchmark-all.py, benchmark-compare.py

---

## Phase 2 — Extended Statistics Toolkit: COMPLETE

### Milestone 2.1 — Filtering + Ancestral + FayWuH + nSL + ROH: COMPLETE
- [x] VariantFilter: shared per-variant QC (AF, call rate, HWE, variant_type)
- [x] VariantQuery: shared Cypher builder (consequence/pathway joins)
- [x] HWEExact: Wigginton 2005 mid-p exact test
- [x] SampleSubsetComputer: on-the-fly stats from packed arrays
- [x] FayWuH: Fay & Wu's H + normalized H (Zeng 2006)
- [x] NSLProcedure: nSL (Ferrer-Admetlla 2014), r=1.000 vs scikit-allel
- [x] ROHProcedure: HMM (Narasimhan 2016) + sliding-window, r=0.96 vs bcftools
- [x] PopSummaryProcedure: whole-chromosome summary → Population nodes
- [x] Ancestral allele annotation (EPO FASTA)
- [x] All 8 existing procedures updated with filter + samples support

### Milestone 2.2 — Complete Statistics Toolkit: COMPLETE
- [x] VectorOps: wcFstComponents (Weir & Cockerham 1984)
- [x] DivergenceProcedure: +fst_wc, +pbs (3-population branch statistic)
- [x] GenomeScanProcedure: +fst_wc, +pbs per window
- [x] GarudH + GarudHProcedure: H1/H12/H2_H1/hap_diversity (sliding-window)
- [x] EHHComputer: instance-configurable (min_ehh, max_gap, gap_scale)
- [x] IHS/NSL: configurable n_af_bins
- [x] Import pipeline: species-agnostic chr lengths from VCF ##contig headers
- [x] 12 procedures, 120 tests

---

## Phase 3 — Storage Optimization + Ploidy + Memory Management: COMPLETE

### Milestone 3.1 — Packed Genotypes (Strategy 4): COMPLETE
- [x] Replaced ~6.3B CARRIES edges with packed byte arrays on Variant nodes
- [x] gt_packed (2 bits/sample), phase_packed (1 bit/sample)
- [x] PackedGenotypeReader: bit-extraction + encoding constants (16 tests)
- [x] All 12 procedures updated for packed arrays
- [x] Python: vcf_parser/csv_emitter updated, CarriesRecord removed
- [x] DB size: ~200 GB → ~46 GB

### Milestone 3.2 — Ploidy-Aware Import: COMPLETE
- [x] VCFParser: ploidy auto-detection (diploid/haploid/mixed)
- [x] ploidy_packed byte[] on Variant nodes
- [x] Sample.sex from panel PED
- [x] Ploidy-aware allele counting (AN/AC)

### Milestone 3.3 — Memory Management Improvements: COMPLETE
- [x] IHS chunking: 5 Mb core + 2 Mb EHH margin (same pattern as XP-EHH)
- [x] IHS unstd_only mode: write only ihs_unstd_{pop}, return SUMMARY row
- [x] HaplotypeMatrix bit-packing: 87% memory reduction (1 bit/haplotype vs 1 byte)
- [x] VectorOps.dPrimePacked: D' from bit-packed arrays using Integer.bitCount
- [x] Python auto-retry: exponential backoff (5s/10s/20s) for transient errors
- [x] Python iHS per-region: 50 Mb regions + AF-bin standardization in Python
- [x] 147 tests passing

---

## Full-Genome Analysis — 1000 Genomes (22 autosomes, 26 populations)

### Phase 1 (per-pop stats): COMPLETE
- [x] 26 populations × 22 chromosomes: diversity, SFS, genome scan, iHS, nSL, ROH

### Phase 2 (pairwise stats): COMPLETE
- [x] 18 key population pairs × 22 chromosomes = 396 XP-EHH tasks (all OK)
- [x] 325 population pairs × 22 chromosomes = 7150 divergence tasks (all OK)

### Phase 3 (ancestral allele analyses): COMPLETE
- [x] Fay & Wu's H, polarized SFS, genome scan

### Phase 4 (population-specific analyses): COMPLETE
- [x] PBS genome scans, Garud's H scans

### Results
- `human_full_results.json` — accumulated results across all phases

---

## Immediate Next Steps

- All 4 analysis phases complete. New JAR (M3.3 memory improvements) deployed 2026-03-16.
- Next: commit all uncommitted changes, then proceed to Phase 4 (kinship/IBS/IBD) or paper writing.

---

## Next Phases (from design doc)

| Phase | Topic | Status |
|-------|-------|--------|
| 4 | Pairwise sample statistics (kinship, IBS/IBD) | Not started |
| 5 | Structure & dimensionality reduction (PCA, community detection) | Not started |
| 6 | Selection inference (composite likelihood, ML classifiers) | Not started |
| 7 | Demographic inference (dadi/moments, ABC) | Not started |
| 8 | Simulation integration (msprime/SLiM) | Not started |
| 9 | AI-augmented inference (GraphRAG, agentic loops) | Not started |

---

## Summary Stats

| Metric | Count |
|--------|-------|
| Stored procedures | 12 |
| Java source files | 31 |
| Java test files | 11 |
| Tests passing | 147 |
| VectorOps methods | 12+ |
| Benchmark datasets | 2 (human 1000G, rice 3K) |
