# Deep Integrative Analysis Plan

## Overview

GraphPop's core value proposition is that all statistics live as persistent properties
on the same graph nodes, enabling multi-statistic queries impossible with flat-file tools.
The current rice interpretation uses ~30% of this potential. This plan describes the
deeper integrative analyses, first for rice (data exists), then replicated for human.

---

## Part A: Rice — Analyses from Existing Results (no Neo4j needed)

These use `results/rice/rice_full_results.json` (91 MB) and `rice_interpretation_report.md`.

### A1. Evolutionary Fingerprint Clustering

**Goal**: Cluster 12 rice populations by multi-statistic evolutionary profiles, not
by genetic distance (Fst). Reveals populations under similar evolutionary regimes
even if genetically distant.

**Statistics per population** (genome-wide means, all already computed):
- π (nucleotide diversity)
- θ_W (Watterson's estimator)
- Tajima's D
- Fay & Wu's H
- piN/piS (from interpretation results)
- F_IS (inbreeding coefficient)
- FROH (ROH fraction)
- n_sweeps (from Garud's H)
- hard_sweep_fraction (H2/H1 < threshold)
- mean |XP-EHH| variance (from existing XP-EHH results)

**Method**: Extract vectors → standardize → PCA + hierarchical clustering → compare
to Fst-based UPGMA tree.

**Prediction**: GJ-tmp and cA-Aus may cluster (both specialist ecotypes with strong
bottlenecks) despite high Fst. Admixed populations should be outliers.

**Script**: `scripts/rice_deep_integration.py integrate_fingerprints`

### A2. Hard vs Soft Sweep Classification at Known Genes

**Goal**: For each known domestication gene (GW5, Hd1, PROG1, Wx, OsC1, DRO1),
determine sweep type from existing Garud's H results.

**Method**: Extract H1, H12, H2/H1 from windows overlapping each gene.
- H2/H1 < 0.05 → hard sweep (single haplotype)
- H2/H1 > 0.10 → soft sweep (multiple haplotypes)

**Biological relevance**: Wx (waxy gene) is known to have been selected independently
in japonica and indica → should show soft sweep. GW5 may show hard sweep.

**Script**: `scripts/rice_deep_integration.py sweep_classification`

### A3. ROH Length Distribution by Population

**Goal**: Decompose FROH into short (<1 Mb), medium (1-5 Mb), and long (>5 Mb) ROH
segments per population. Different lengths indicate different demographic events.

**Method**: Parse ROH segments from existing results, bin by length class.
- Long ROH → recent inbreeding/self-pollination
- Medium ROH → recent bottleneck (tens of generations)
- Short ROH → ancient bottleneck/founder effect

**Prediction**: GJ-tmp should have predominantly medium ROH (temperate adaptation
bottleneck ~4-5 kya), not long ROH (would indicate recent selfing).

**Script**: `scripts/rice_deep_integration.py roh_distribution`

### A4. Cross-Statistic Correlation Matrix

**Goal**: Compute correlations between all pairs of population-level statistics.
Reveals which evolutionary forces co-vary.

**Expected patterns**:
- π and piN/piS: negative (larger Ne → more efficient purifying selection → lower piN/piS)
- FROH and Tajima's D: both reflect bottlenecks, should correlate
- piN/piS and FROH: positive (bottleneck + drift → more load)
- n_sweeps and |Tajima's D|: positive (more sweeps → more negative D)

**Script**: `scripts/rice_deep_integration.py stat_correlations`

---

## Part B: Rice — New Neo4j Queries Required

These require the rice database loaded in Neo4j.

### B1. Annotation-Conditioned Tajima's D and Fay & Wu's H

**Goal**: Is purifying selection detectable BEYOND the piN/piS ratio? If Tajima's D
is more negative for missense than synonymous variants, it directly shows purifying
selection acting on protein-coding variants independent of demographic effects.

**Cypher** (per population, per chromosome):
```cypher
// Missense diversity (includes Tajima's D, Fay & Wu's H)
CALL graphpop.diversity('Chr1', 1, 43270923, 'GJ-tmp',
    {consequence: 'missense_variant'})

// Synonymous diversity
CALL graphpop.diversity('Chr1', 1, 43270923, 'GJ-tmp',
    {consequence: 'synonymous_variant'})
```

**Analysis**: For each population, compare:
- Tajima's D (missense) vs Tajima's D (synonymous) → Δ_D
- Fay & Wu's H (missense) vs Fay & Wu's H (synonymous) → Δ_H

If Δ_D < 0 and Δ_H < 0: purifying selection preferentially removes missense variants.
GJ-tmp should show the LARGEST Δ (strongest selection pressure).

**Scale**: 12 pops × 12 chr × 2 consequences = 288 calls (same as piN/piS, but we
already have the synonymous/missense results — we just need to extract Tajima's D
and Fay & Wu's H from them, not just π).

**IMPORTANT**: The piN/piS analysis already called `graphpop.diversity` with
consequence filters. Those results include Tajima's D and Fay & Wu's H. We should
CHECK if these fields are already in the existing rice interpretation results.
If so, no new Neo4j queries needed — just deeper extraction from existing JSON.

**Script**: `scripts/rice_deep_integration.py annotation_conditioned_neutrality`

### B2. Multi-Statistic Evidence Table at Selection Loci

**Goal**: For the top 20 selection candidate loci, query ALL statistics from the graph
simultaneously. This is THE showcase of GraphPop's persistent analytical record.

**Cypher** (for each candidate locus, e.g., Hd1 at Chr6:9.0 Mb):
```cypher
// Get all variant-level statistics at a locus
MATCH (v:Variant)
WHERE v.chr = 'Chr6' AND v.pos >= 8900000 AND v.pos <= 9600000
OPTIONAL MATCH (v)-[c:HAS_CONSEQUENCE]->(g:Gene)
RETURN v.variantId AS variant,
       v.pos AS pos,
       v.ihs_GJ_tmp AS ihs,
       v.xpehh_GJ_trp_XI_1A AS xpehh,
       v.nsl_GJ_tmp AS nsl,
       g.symbol AS gene,
       c.consequence AS consequence,
       c.impact AS impact
ORDER BY v.pos
```

```cypher
// Get window-level statistics at the same locus
MATCH (w:GenomicWindow)
WHERE w.chr = 'Chr6' AND w.start >= 8900000 AND w.end <= 9600000
RETURN w.start AS start, w.end AS end,
       w.pi_GJ_tmp AS pi,
       w.fst_wc_GJ_tmp_XI_1A AS fst,
       w.tajima_d_GJ_tmp AS tajima_d,
       w.fay_wu_h_GJ_tmp AS fay_wu_h,
       w.pbs_GJ_trp AS pbs,
       w.h12_GJ_tmp AS h12,
       w.h2_h1_GJ_tmp AS h2_h1
ORDER BY w.start
```

**Output**: A table per locus with columns for every statistic, gene annotations,
and sweep classification. This table goes directly into Supplementary Materials
and demonstrates GraphPop's unique capability.

**Candidate loci** (from existing results):
1. GW5 (Chr5:5.26 Mb) — grain width
2. Hd1 (Chr6:8.95 Mb) — heading date
3. PROG1 (Chr7:3.18 Mb) — prostrate growth
4. Wx (Chr6:1.45 Mb) — waxy/amylose
5. OsC1 (Chr6:4.95 Mb) — hull color
6. DRO1 (Chr9:22.04 Mb) — deep rooting
7. Top missense-Fst outlier (Chr1:7.45 Mb)
8-20. Top XP-EHH / PBS peaks not overlapping known genes

**Script**: `scripts/rice_deep_integration.py convergence_table`

### B3. Annotation-Conditioned Selection Scans (iHS, XP-EHH, nSL)

**Goal**: Do missense variants show stronger selection signals than synonymous?
This would directly demonstrate protein-level selection.

**Cypher**:
```cypher
// iHS restricted to missense variants
CALL graphpop.ihs('Chr1', 'GJ-tmp', {consequence: 'missense_variant'})

// iHS restricted to synonymous variants
CALL graphpop.ihs('Chr1', 'GJ-tmp', {consequence: 'synonymous_variant'})
```

**Analysis**: Compare distributions of |iHS|, |XP-EHH|, |nSL| between missense
and synonymous variants. If mean |iHS| for missense > synonymous, selection acts
preferentially on protein function.

**Caveat**: Annotation-conditioned EHH statistics may have reduced power (fewer
variants per haplotype). Need to verify sufficient variant density.

**Scale**: Start with 2-3 populations × 2-3 chromosomes as pilot.

**Script**: `scripts/rice_deep_integration.py annotation_conditioned_selection`

### B4. LD Decay by Population

**Goal**: Characterize LD landscape per population. Faster LD decay = larger
historical Ne = more efficient selection.

**Cypher**:
```cypher
// LD for a region (writes LD edges, returns r² and distance)
CALL graphpop.ld('Chr1', 1, 1000000, 'GJ-tmp',
    {max_distance: 500000, min_r2: 0.01})
```

**Analysis**: Bin r² by physical distance → fit LD decay curve → estimate half-decay
distance per population. Compare: indica (fast decay) vs japonica (slow decay).

**Scale**: Sample 3-5 regions per chromosome × 12 populations. LD is slow — limit
to representative regions.

**Script**: `scripts/rice_deep_integration.py ld_decay`

### B5. Joint SFS: Demographic History × Annotation

**Goal**: Compare the joint SFS for missense vs synonymous variants between key
population pairs. Different shapes reveal whether selection has distorted the
demographic signal.

**Cypher**:
```cypher
// Joint SFS for all variants
CALL graphpop.joint_sfs('Chr1', 1, 43270923, 'GJ-tmp', {pop2: 'XI-1A'})

// Joint SFS for missense only
CALL graphpop.joint_sfs('Chr1', 1, 43270923, 'GJ-tmp',
    {pop2: 'XI-1A', consequence: 'missense_variant'})

// Joint SFS for synonymous only
CALL graphpop.joint_sfs('Chr1', 1, 43270923, 'GJ-tmp',
    {pop2: 'XI-1A', consequence: 'synonymous_variant'})
```

**Analysis**: Compare shapes. Synonymous jSFS reflects demography. Missense jSFS
shows distortion from selection. The DIFFERENCE is the selection signal, free of
demographic confounding.

**Scale**: 3 key pairs × 12 chr × 3 consequence types = 108 calls.

**Script**: `scripts/rice_deep_integration.py joint_sfs_annotation`

---

## Part C: Human — Primary Annotation-Conditioned Analyses

These are the BLOCKING analyses needed before the paper. Uses existing
`scripts/human_interpret.py` subcommands.

### C1. piN/piS for All 26 Populations (human_interpret.py pinsps)

**Goal**: First systematic genome-wide quantification of purifying selection
efficacy across all 26 human subpopulations.

**Method**: Already implemented in human_interpret.py cmd_pinsps().
288 → 1,144 calls (26 pops × 22 chr × 2 consequences).

**Predictions**:
- African populations (YRI, LWK, GWD): lowest piN/piS (largest Ne)
- Out-of-Africa populations (CEU, CHB): higher piN/piS (bottleneck)
- Admixed (MXL, PUR, CLM, PEL): intermediate or elevated

**Command**: `python scripts/human_interpret.py pinsps`

### C2. Consequence-Stratified Genome Scans (human_interpret.py gscan)

**Goal**: Missense vs synonymous Fst at continental boundaries.

**Already implemented** for 3 key pairs (YRI-CEU, CEU-CHB, YRI-CHB).

**Command**: `python scripts/human_interpret.py gscan`

### C3. Gene Annotation at Selection Peaks (human_interpret.py annotate)

**Goal**: Cross-reference XP-EHH and Garud's H peaks with Gene nodes.

**Command**: `python scripts/human_interpret.py annotate`

### C4. PBS Genome Scans (human_interpret.py pbs)

**Goal**: Population branch-specific selection for 6 focal populations.

**Command**: `python scripts/human_interpret.py pbs`

### C5. Additional Subcommands

```bash
python scripts/human_interpret.py fay_wu      # Fay & Wu's H
python scripts/human_interpret.py usfs         # Unfolded SFS
python scripts/human_interpret.py hwscan       # Garud's H genome scan
python scripts/human_interpret.py roh_hmm      # ROH HMM
python scripts/human_interpret.py daf_enrichment  # DAF at peaks
python scripts/human_interpret.py report       # Generate Markdown report
```

**Execution order**: pinsps → gscan → annotate → pbs → fay_wu → usfs → hwscan →
roh_hmm → daf_enrichment → report

**Estimated time**: 8-16 hours total (most time in pinsps and gscan).

---

## Part D: Human — Deeper Integrative Analyses

Same framework as rice Part A/B, replicated for human.

### D1. Out-of-Africa Load Gradient

**Goal**: Plot piN/piS vs distance-from-Africa for all 26 populations. Classic
population genetics prediction: serial bottlenecks during range expansion
should create a gradient of increasing genetic load.

**Data needed**: piN/piS (from C1) + geographic coordinates (from 1000G metadata).

**Prediction**: Strong positive correlation between piN/piS and geographic distance
from East Africa.

### D2. Evolutionary Fingerprints (26 Human Populations)

**Goal**: Same as rice A1 but for 26 populations. Cluster by evolutionary regime.

**Predictions**:
- AMR populations cluster separately (admixed, unique demographic history)
- AFR populations form tight cluster (large Ne, similar regimes)
- EUR/EAS may overlap (similar bottleneck history despite different geography)

### D3. Multi-Statistic Evidence at Known Selection Loci

**Target loci** (from existing KNOWN_GENES in human_interpret.py):
1. LCT (Chr2:136.6 Mb) — lactase persistence (EUR)
2. SLC24A5 (Chr15:48.4 Mb) — skin pigmentation (EUR)
3. EDAR (Chr2:109.5 Mb) — hair/teeth (EAS)
4. DARC (Chr1:159.2 Mb) — malaria resistance (AFR)
5. HBB (Chr11:5.2 Mb) — sickle cell (AFR)
6. FADS1/FADS2 (Chr11:61.6 Mb) — fatty acid metabolism
7. HLA region (Chr6:29-33 Mb) — immune function (balancing selection)
8. ADH1B (Chr4:99.3 Mb) — alcohol metabolism (EAS)

**For each**: full evidence table (iHS, XP-EHH, nSL, PBS, Fay & Wu's H,
Garud's H, local π, Tajima's D, gene annotation, consequence type).

### D4. Annotation-Conditioned Neutrality Tests

Same as rice B1: compare Tajima's D and Fay & Wu's H for missense vs synonymous
across all 26 populations.

### D5. Immune Pathway Analysis

**Goal**: Are immune-related genes under stronger selection?

**Cypher**:
```cypher
// Diversity restricted to immune pathway genes
CALL graphpop.diversity('chr6', 28000000, 34000000, 'CEU',
    {pathway: 'Immune System'})
```

Compare piN/piS for immune pathway vs genome background.

---

## Implementation: Single Script Architecture

All analyses implemented in one script with subcommands:

```
scripts/rice_deep_integration.py
    integrate_fingerprints   — Part A1
    sweep_classification     — Part A2
    roh_distribution         — Part A3
    stat_correlations        — Part A4
    annotation_conditioned_neutrality  — Part B1
    convergence_table        — Part B2
    annotation_conditioned_selection   — Part B3
    ld_decay                 — Part B4
    joint_sfs_annotation     — Part B5
    report                   — Generate Markdown report

scripts/human_deep_integration.py
    (same subcommands, adapted for human populations)
```

---

## Execution Order

### Phase 1: Rice from existing data (no Neo4j, hours)
1. A1: Evolutionary fingerprints
2. A2: Sweep classification at known genes
3. A3: ROH length distribution
4. A4: Cross-statistic correlations

### Phase 2: Human primary analyses (Neo4j required, 8-16 hours)
5. C1: piN/piS (1,144 calls)
6. C2: Genome scans (annotation-conditioned)
7. C3-C5: Annotation, PBS, Fay & Wu's H, etc.

### Phase 3: Rice new Neo4j queries (switch to rice DB, 4-8 hours)
8. B1: Annotation-conditioned Tajima's D / Fay & Wu's H
9. B2: Multi-statistic convergence table
10. B3: Annotation-conditioned iHS/XP-EHH (pilot)
11. B4: LD decay (sampled regions)
12. B5: Joint SFS annotation-conditioned

### Phase 4: Human deeper integration (switch to human DB, 4-8 hours)
13. D1-D5: Same framework as rice

### Phase 5: Interpretation reports
14. Rice deep integration report
15. Human deep integration report
16. Update manuscript with new findings

---

## Key Design Principle

Every analysis in this plan exploits GraphPop's unique capability:
**querying multiple statistics and annotations on the same graph nodes**.
No analysis here is achievable with a single existing tool.
The interpretation should always highlight WHY graph-native computation
enabled the discovery, not just WHAT was discovered.
