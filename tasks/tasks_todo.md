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

## Deep Functional Investigations (Novel Graph-Native) — COMPLETE 2026-03-21

### Investigations 1–10: COMPLETE
- [x] Inv.01: pathway_fst — cardiac repolarization top pathway
- [x] Inv.02: convergent_sweeps — 9 genes with multi-population sweeps
- [x] Inv.03: roh_sweep_correlation — ρ=-0.75 π vs sweeps
- [x] Inv.04: kcne1_deepdive — ancient pre-OoA sweep (h12>0.9, FST ratio=0.296)
- [x] Inv.05: go_enrichment — per-population GO enrichment
- [x] Inv.06: synthesis — multi-evidence synthesis (HIGH/MED/LOW tiers)
- [x] Inv.07: ooa_gradient — PCA gradient (PC1=56%, PC2=30%)
- [x] Inv.08: gene_fst — 21,194 genes; SLC24A5 rank 10, max_FST=0.843
- [x] Inv.09: pbs_scan — 6 pops; YRI=HLA, GIH=SLC24A5
- [x] Inv.10: temporal_selection — FROH, H12, iHS temporal stratification

### Investigations 11–15: COMPLETE
- [x] Inv.11: consequence_bias — HIGH FST=0.124 < LOW FST=0.143 (purifying selection on HIGH-impact); 891k records, 214s
- [x] Inv.12: pathway_coselection — 2232 pathways, 1.18M co-selection edges (|ρ|>0.7), 3 communities; 41s
- [x] Inv.13: cross_species — module-level convergence (development_growth top); rice genes not in Neo4j, redesigned to keyword-matching; 1s
- [x] Inv.14: rare_burden — 284,945 private-rare events, 27,899 genes; TTN top; FROH correlation ρ=0.11 (ns); 69s
- [x] Inv.15: sweep_extent — 9 sweep genes, KCNE1 largest 50kb; FST decay from variant_cache; 3.4s
- [x] Figures fig14–fig18 added to generate_supp_figures.py

### Key findings (Inv.11–15)
- HAS_CONSEQUENCE schema: only impact/consequence/feature_type/feature (NO CADD/SIFT props)
- Min AC in 1000G coding variants = 2 (singletons excluded); private rare threshold must be AF<2%
- Rice genes (LOC_Os*) are NOT in Neo4j; cross-species requires module-level functional analysis
- ihs_{POP} sparse on Variant nodes for 19 non-African pops only

---

## Investigation 16b — Individual-Level Evolutionary Trajectory: COMPLETE 2026-03-21

### Inv.16b: individual_trajectory — COMPLETE
- [x] Genome-wide FROH per sample from roh_hmm checkpoint (3094/3202 samples with ROH)
- [x] Chr22 individual features: het_rate, rare_burden, rare_missense (chunked hap_cache, 88s)
- [x] PCA (PC1=57%, PC2=31%): PC1=diversity axis, PC2=inbreeding axis
- [x] UMAP: 3202 individuals in clear continental group clusters
- [x] Velocity components: drift_inbreeding, functional_burden, population_diversity
- [x] Panel F: individual PC1 vs population PC1 correlation ρ=0.74 (validates approach)
- [x] Outputs: individual_trajectory.tsv/.json, fig20_individual_trajectory.png

### Key findings (Inv.16b)
- African: highest het_rate (0.0453), lowest FROH (0.0047), highest rare_missense (34.5/indiv)
- South Asian: highest FROH (0.0136), moderate het_rate (0.0354)
- FROH vs het_rate: ρ=−0.394, p=1.2e−29 (genome drift ↔ local chr22 diversity)
- Note: chr22 has no H12 GenomicWindow nodes → sweep_burden replaced by rare_missense_chr22
- Note: hap_cache only has chr22; chr22-based features are representative proxies

---

## Investigation 17 — Novel Graph-Native Individual Analyses: COMPLETE 2026-03-21

### Inv.17: novel_ind_analyses (Analysis 4 + Analysis 6) — COMPLETE
- [x] Analysis 4: Outlier individuals — z-distance from population centroid in PC space (3202 samples)
  - Top outlier: HG01816 (CDX, z=9.50, FROH=0.125 = 18× CDX mean — likely recent consanguinity)
  - 4 UNKNOWN-pop samples in top 10; z-distance weakly predicts FROH (ρ=-0.087)
- [x] Analysis 6: Purifying selection efficiency gradient — pathway-constrained vs divergent rare missense
  - 55,446 Neo4j variant-pathway rows → 5,833 chr22 positions classified
  - Constrained:Divergent ratio: African=0.187, East Asian=0.086 (2.2× difference)
  - ρ(FROH, C:D ratio)=-0.065 p=2e-4 — bottleneck predicts preferential constrained pathway depletion
- [x] Figure fig21: Purifying selection gradient (3 panels)
- [x] Figure fig22: Outlier individual analysis (2 panels)
- [x] Source data files in paper/source_data/
- [x] paper/novel_graph_analyses.md updated with confirmed findings

### Paper preparation document: novel_graph_analyses.md
- [x] All 6 analyses documented with research question, method, Cypher queries, findings, novelty
- [x] "What's Genuinely Novel" table (main text candidate)
- [x] Analysis 4 findings updated with actual confirmed results
- [x] Analysis 6 findings updated with actual confirmed results
- [ ] Analyses 1, 2, 3, 5 remain planned/feasible (need full hap_cache or further development)

### Figures status
- [x] fig19_evo_trajectory.png — population-level (Inv.16)
- [x] fig20_individual_trajectory.png — individual UMAP (Inv.16b)
- [x] fig21_purifying_selection_gradient.png — main text candidate (Inv.17/Analysis 6)
- [x] fig22_outlier_analysis.png — supplementary (Inv.17/Analysis 4)
- [ ] fig23 — pathway trajectory embedding (Analysis 5, planned)

### Paper section placement
- Main text (application.tex, human 1000G section): fig21 + "What's Novel" table
- Supplementary: fig20, fig22, individual_trajectory source data

---

## Immediate Next Steps

- All 4 analysis phases complete. New JAR (M3.3 memory improvements) deployed 2026-03-16.
- Next: commit all uncommitted changes, then proceed to Phase 4 (kinship/IBS/IBD) or paper writing.

### [x] Redo benchmark (scripts/benchmark-vs-tools.py) — DONE 2026-03-20
Results: benchmarks/results/benchmark_v4.jsonl (all 4 regions, exit 0)
Added Group F: Joint nSL+iHS+XP-EHH (sel_stats_numpy_joint) — 3.6× speedup on full chr22.
Key full-chr22 speedups vs scikit-allel:
  π/θ_W/D: 145×  |  Fst/Dxy: 245×  |  SFS: 327×  |  iHS: 179×  |  XP-EHH: 136×
GraphPop-numpy vs scikit-allel:
  iHS: 90×  |  XP-EHH: 53×  |  Joint parallel: 3.6× vs independent numpy
Note: pairwise LD gap (PLINK2-pgen 0.5–1.1s vs GraphPop 1.3–8.8s) remains; deferred
until numpy LD cache is implemented. selscan skipped (binary not returning results).
TODO: update paper/manuscript/benchmarks.tex + regenerate paper/scripts/fig2_benchmarks.py.

---

## Deferred Work

### GenomicWindow write-back (deferred to Phase 4)
The numpy gscan bypass writes results to the JSON checkpoint only (not to Neo4j).
Writing ~400k `GenomicWindow` nodes back to the graph is deferred until Phase 4+,
when Cypher-based traversals over window results are actually needed (e.g. GNN
embeddings, graph-based selection signal queries). At that point, implement a
standalone `import_gscan_windows.py` that reads from `human_interpretation_results.json`
and bulk-imports via batched `UNWIND ... MERGE`. Do NOT add this to `human_interpret.py`.

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
