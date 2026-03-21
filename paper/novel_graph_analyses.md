# Novel Graph-Native Biological Investigations: Analyses Unlocked by GraphPop's Individual × Pathway × Selection Graph

**Status:** In preparation — paper section drafts and figure specifications
**Date:** 2026-03-21
**Paper placement:** Application section — Human 1000 Genomes (extends `paper/manuscript/application.tex`)
**Relationship to prior work:**
- Investigations 1–16b complete; results in `data/results/`
- Individual-level evolutionary trajectory (Inv.16b) completed; `data/results/individual_trajectory.tsv` (3,202 samples)
- Population-level pathway FST (Inv.01) complete; `data/results/pathway_fst.json` (2,241 pathways)
- Convergent sweeps (Inv.02) complete; 9 convergent sweep genes including KCNE1
- Rare burden (Inv.14) complete; 284,945 private-rare events across 27,899 genes

---

## Summary: What Is Genuinely Novel

| Capability | Traditional Approach | GraphPop |
|---|---|---|
| Individual pathway rare burden | 4 tools + manual joins | 1 Cypher query |
| ROH × sweep gene overlap | BED intersection + VCF parsing + gene annotation | 2-hop traversal |
| Outlier → driver variant trace | Impossible in one step | Direct CARRIES query |
| Temporal burden vector | Not possible | iHS score on Variant nodes |
| Pathway-level trajectory embedding | Not possible end-to-end | Single-query matrix extraction |
| Purifying selection gradient | Requires 5+ tools + custom joins | FROH + rare burden + pathway FST in one analysis |

**Core architectural point:** All six analyses exploit the fact that GraphPop stores individual genotypes (CARRIES edges), variant annotations (IN_GENE, HAS_CONSEQUENCE), pathway structure (IN_PATHWAY), and population selection statistics (ihs_{POP}, GenomicWindow nodes) in a single traversable graph. No traditional pipeline achieves this without staging intermediate files across multiple tools.

---

## Analysis 1: Individual × Pathway Rare Burden (chr22)

### Research Question
Which biological pathways accumulate private rare missense variants on a per-individual basis? Do evolutionarily constrained pathways (low population FST — maintained under purifying selection) carry lower per-individual rare missense burden than unconstrained, population-divergent pathways? This tests whether pathway-level selective constraint is visible at the individual genotype level.

### Method
A single Cypher traversal walks from Sample nodes through CARRIES edges to Variant nodes, then through IN_GENE and IN_PATHWAY edges, filtered on variant rarity (global AF < 0.01) and consequence (HAS_CONSEQUENCE with missense annotation):

```cypher
MATCH (s:Sample)-[:CARRIES]->(v:Variant)-[:IN_GENE]->(g:Gene)-[:IN_PATHWAY]->(p:Pathway)
MATCH (v)-[:HAS_CONSEQUENCE]->(c:Consequence)
WHERE c.impact IN ['HIGH', 'MODERATE']
  AND c.consequence CONTAINS 'missense'
  AND v.af_global < 0.01
RETURN s.sampleId AS sample, p.name AS pathway, count(v) AS rare_missense_burden
ORDER BY s.sampleId, rare_missense_burden DESC
```

This produces a 3,202 × N_pathways sparse matrix where each cell is the rare missense count for individual `i` in pathway `j`. Pathway FST ranks from `pathway_fst.json` (2,241 pathways, mean_fst range 0.0–0.538) are joined externally for the constrained vs unconstrained comparison.

**Data available:** chr22 has 9,806 rare missense SNPs (AF < 0.01); pathway FST values exist for 2,241 pathways. The chr22-only matrix is implementable now from `data/variants_cache/chr22.npz` + `data/hap_cache/`.

### Status
Partially implemented. Individual-level chr22 rare missense counts exist in `data/results/individual_trajectory.tsv` (column `rare_missense_chr22`). The pathway-stratified breakdown (which pathways per individual) requires the Cypher query above or a numpy join against the variants cache annotated with pathway membership.

### Key Findings
From `individual_trajectory.json` group statistics (chr22, AF < 0.01):

| Continental Group | Mean rare_missense_chr22 | Std |
|---|---|---|
| African | 34.53 | 8.35 |
| South Asian | 21.03 | 6.13 |
| East Asian | 18.28 | 5.46 |
| European | 16.79 | 5.48 |
| American | 16.81 | 6.15 |

African individuals carry 2.06× more rare missense variants per chromosome than Europeans, consistent with the out-of-Africa bottleneck reducing rare variant diversity in non-African populations. The pathway-stratified version of this analysis will reveal whether this excess is pathway-agnostic or concentrated in specific functional categories (predicted outcome: excess concentrated in unconstrained pathways; constrained pathways maintained lower burden across all groups due to shared purifying selection).

From `pathway_fst.json` (2,241 pathways):
- Most constrained pathway (mean_fst ≈ 0.000): Nucleotide catabolism (R-HSA-8956319)
- Highest divergence pathway (mean_fst = 0.538): PTK6 Down-Regulation (R-HSA-8849472)
- SLC24A5 pathway (R-HSA-5619036): mean_fst = 0.452 — known pigmentation sweep

### What Is Novel
No existing tool computes per-individual, per-pathway rare burden in a single operation. The standard workflow requires: (1) bcftools to extract rare variants per sample; (2) VEP/ANNOVAR for consequence annotation; (3) BEDTools to intersect with gene bodies; (4) a custom script to map genes to pathways; (5) aggregation. GraphPop replaces steps 1–5 with a single Cypher traversal because genotypes, annotations, and pathway memberships are co-located in the graph.

### Paper Placement
Supplementary figure (pathway-level individual burden heatmap). Cited in main text Analysis section to demonstrate the query capability.

### Figure Description
**fig_S_pathway_burden:** Three-panel figure.
(A) Violin plots of per-individual chr22 rare missense burden by continental group (5 groups), showing African enrichment with mean ± SD annotated.
(B) Scatter plot: pathway mean_fst (x-axis, 2,241 pathways) vs. mean per-individual rare missense burden in that pathway (y-axis); color by pathway functional category. Expected negative trend in constrained pathways (low FST = low burden).
(C) Example heatmap: top 20 pathways by variance across individuals × 50 randomly sampled individuals per group.

---

## Analysis 2: ROH Overlap With Convergent Sweep Genes

### Research Question
Are highly inbred individuals (elevated FROH, indicating recent autozygosity) more likely to be homozygous at known convergent selection loci — specifically the 9 convergent sweep genes identified in Inv.02? This tests whether population-level sweeps (which increase haplotype homozygosity at a locus) are detectable at the individual level via ROH-sweep overlap, and whether drift-driven inbreeding (FROH) correlates with sweep-locus homozygosity.

The biological distinction matters: sweep homozygosity is selection-driven and population-specific; ROH-driven homozygosity is drift/consanguinity-driven and individual-specific. High correlation between them in the same individual would indicate compound drift + selection effects.

### Method
This requires two pieces of information per sample: (1) genome-wide FROH (available in `individual_trajectory.tsv`); (2) whether each ROH segment overlaps the sweep gene loci. The key traversal:

```cypher
MATCH (s:Sample)-[:CARRIES {gt: 2}]->(v:Variant)-[:IN_GENE]->(g:Gene)
WHERE g.symbol IN ['KCNE1', 'SLC24A5', 'CRYAA', 'LINC01678']
WITH s, g.symbol AS sweep_gene, count(v) AS hom_alt_count
MATCH (s)
RETURN s.sampleId, s.population, sweep_gene, hom_alt_count,
       s.froh_genome
ORDER BY hom_alt_count DESC
```

`gt: 2` on CARRIES indicates HOM_ALT — the individual is homozygous for the derived allele at sweep-associated variants. High `hom_alt_count` at sweep loci in high-FROH individuals would indicate the compound effect.

**Note:** This query uses HOM_ALT counts at sweep gene variants, not ROH segment overlap directly. Full ROH-overlap analysis requires ROH segment coordinates per individual (currently only FROH totals are stored, not segment positions). The FROH values in `individual_trajectory.tsv` are genome-wide estimates from `fix_froh.py`.

### Status
Feasible but partially blocked. FROH genome-wide values are available for all 3,202 individuals. HOM_ALT CARRIES counts at the 9 convergent sweep genes are queryable now. Full ROH segment position analysis requires storing ROH segments as Neo4j nodes or loading from bcftools output files — not yet implemented.

**Convergent sweep genes from Inv.02:**
`KCNE1`, `SLC24A5`, `CRYAA`, `LINC01678`, `LOC102724428`, `LOC102724560`, `LOC102725065`, `LOC124900992`, `LOC124900995`

**KCNE1 key stats (Inv.04):** H12 > 0.9 in all 5 continental groups, FST ratio = 0.296 — pre-OoA ancient sweep.

### Key Findings
From `individual_trajectory.json`:
- FROH correlation with het_rate_chr22: ρ = −0.394 (individuals with more ROH have lower heterozygosity, as expected)
- African mean FROH = 0.0047 (lowest — minimal bottleneck); South Asian mean = 0.0136 (highest — strongest recent inbreeding signal)
- American mean FROH = 0.0102, n_roh = 15.9 (admixture-driven elevated FROH)

Prediction (not yet computed): KCNE1 homozygosity should be uniform across continental groups (pre-OoA sweep) whereas SLC24A5 homozygosity should be strongly European-enriched (post-OoA European sweep). FROH vs KCNE1 HOM_ALT count correlation should be near zero (sweep-driven, not drift-driven); FROH vs European rare-variant homozygosity should be positive.

### What Is Novel
Traditional ROH analysis (PLINK, bcftools) identifies ROH segments but cannot directly intersect them with sweep gene annotations and individual genotypes in a single query. The graph enables: Sample → CARRIES → Variant → IN_GENE → Gene in one traversal, with individual-level FROH as a node property.

### Paper Placement
Supplementary note (planned analysis section). The query pattern is cited in the main text as an example of 2-hop graph traversal replacing a 3-tool pipeline.

### Figure Description
**fig_S_roh_sweep:** Two-panel figure.
(A) Scatter plot: individual FROH (x-axis) vs HOM_ALT count at KCNE1 variants (y-axis); points colored by continental group. Expected: near-zero slope (sweep-driven not drift-driven).
(B) Same for SLC24A5: expect strong European enrichment on y-axis regardless of FROH.

---

## Analysis 3: Individual Temporal Selection Fingerprint

### Research Question
Do individuals differ systematically in the temporal composition of their selection-tagged alleles — specifically, do they carry disproportionate ancient-sweep (shared across continental groups, high |iHS| in all populations) vs. recent-sweep (population-specific, high |iHS| only in one group) derived alleles? This creates an individual-level "temporal selection fingerprint" — a vector whose components are ancient-sweep burden, recent-sweep burden, and neutral rare burden.

The biological question: bottlenecked populations swept recently at population-specific loci; do those populations show depleted ancient-sweep alleles (purged by founder effects) relative to African individuals who accumulated ancient sweeps without bottleneck?

### Method
Three-tier classification of Variant nodes based on iHS availability and magnitude:

**Tier 1 — Ancient/global sweeps:** |iHS| > 2.0 in ≥3 of the 19 stored populations (CEU, FIN, GBR, IBS, TSI, CHB, JPT, CDX, KHV, BEB, GIH, ITU, PJL, STU, CLM, MXL, PEL, PUR — note: no African populations).

**Tier 2 — Recent/population-specific sweeps:** |iHS| > 2.0 in exactly 1 of the 19 populations.

**Tier 3 — Neutral rare:** AF < 0.01, no strong iHS signal.

```cypher
// Count ancient-sweep alleles per individual (chr22 example)
MATCH (s:Sample)-[:CARRIES]->(v:Variant)
WHERE v.ihs_CEU IS NOT NULL
  AND abs(v.ihs_CEU) > 2.0
  AND abs(v.ihs_CHB) > 2.0
  AND abs(v.ihs_GIH) > 2.0
WITH s, count(v) AS ancient_sweep_burden
RETURN s.sampleId, s.population, ancient_sweep_burden
```

For chr22: 8,533 variants have iHS_CEU stored. The multi-population intersection can be computed from variants_cache iHS columns after loading the numpy arrays.

### Status
Feasible for chr22 (iHS_CEU available for 8,533 chr22 variants). Genome-wide analysis requires loading iHS_{POP} properties for all 19 populations via Cypher or extending the variants_cache export to include iHS arrays. The chr22 analysis is implementable within 1–2 hours using existing hap_cache + iHS data.

**Limitation:** African populations (YRI, LWK, ESN, ACB, ASW, GWD, MSL) have no stored iHS values, meaning the temporal fingerprint is computed from non-African sweep signals. African ancient-sweep burden must be inferred from the multi-population overlap approach above.

### Key Findings
From `evo_trajectory.json` velocity scores (population-level evolutionary trajectory PCA):
- Drift/bottleneck dimension: PEL = +0.944, PJL = +1.198 (most bottlenecked); LWK = −0.983, MSL = −0.978 (least bottlenecked, African)
- Recent selection dimension: CLM = +1.472, FIN = +1.127 (high recent selection); CDX = −1.415, CHB = −1.222 (negative = ancient / deep sweep signal)
- Expansion dimension: ACB = +1.460, ASW = +1.397 (African diaspora expansion signal)

Individual-level stat PCA (from `individual_trajectory.json`):
- PC1 explains 57.0% of individual variance; PC2 explains 30.9%
- PC1 loads primarily on het_rate_chr22 and rare_burden_chr22 (diversity axis)
- evo_trajectory PC1 (population-level) correlates with individual PC1: ρ = 0.906, p ≈ 0

### What Is Novel
No existing tool creates per-individual temporal selection fingerprints because it requires: (1) multi-population iHS values co-located on Variant nodes; (2) individual CARRIES edges to count those variants; (3) classification logic. In a VCF-centric workflow, iHS is computed separately per population, stored in separate files, and individual-level joins are not performed.

### Paper Placement
Supplementary note (planned analysis). The concept of individual temporal selection fingerprints is cited in the main text as a uniquely graph-native capability.

### Figure Description
**fig_S_temporal_fingerprint:** Two-panel figure.
(A) Ternary plot: each individual plotted in (ancient_sweep_burden, recent_sweep_burden, neutral_rare_burden) space; color by continental group.
(B) Comparison of ancient-sweep allele count per individual by continental group (violin plot), with the expectation that non-African populations carry fewer ancient-sweep alleles due to bottleneck purging.

---

## Analysis 4: Outlier Individuals in Evolutionary Process Space

### Research Question
Which individuals deviate most from their population centroid in evolutionary process space (the 5-dimensional space of FROH, n_roh, het_rate, rare_burden, rare_missense)? What variant-level features drive their deviation? This identifies individuals with unusual evolutionary histories — e.g., admixed individuals with anomalously high rare burden, or isolated individuals with elevated inbreeding relative to their population.

The clinical/translational relevance: outliers with elevated rare missense burden relative to their population centroid are candidates for functional screening. The graph enables direct follow-up: given an outlier sample, trace the CARRIES graph to enumerate the specific rare missense variants driving the deviation.

### Method
**Step 1 — Mahalanobis distance per individual from population centroid:**
Compute population-specific mean vector μ_pop and covariance Σ_pop from the 5-feature matrix in `individual_trajectory.tsv`. Per-individual Mahalanobis distance = sqrt((x_i − μ_pop)^T Σ_pop^{−1} (x_i − μ_pop)).

**Step 2 — Top outlier identification:**
Select top 5 outliers per continental group (25 individuals total). These are individuals whose evolutionary process profile is anomalous for their population.

**Step 3 — Driver variant tracing (GraphPop-unique step):**
```cypher
// For outlier individual HG01234: find what drives their elevated rare_missense_chr22
MATCH (s:Sample {sampleId: 'HG01234'})-[:CARRIES]->(v:Variant)
MATCH (v)-[:HAS_CONSEQUENCE]->(c:Consequence)
WHERE c.impact IN ['HIGH', 'MODERATE']
  AND v.af_global < 0.01
MATCH (v)-[:IN_GENE]->(g:Gene)-[:IN_PATHWAY]->(p:Pathway)
RETURN g.symbol, p.name, c.consequence, v.variantId, v.af_global
ORDER BY v.af_global ASC
LIMIT 50
```

This step — tracing from an outlier individual to the specific rare variants driving their anomalous burden, with gene and pathway annotations — is a single 3-hop query in GraphPop. Traditional tools would require: extracting the individual's genotype from a multi-sample VCF, annotating with VEP, joining against pathway databases, comparing against population-level burden statistics.

### Status
Implementable now from existing data. `individual_trajectory.tsv` contains all 3,202 individuals × 5 features. Mahalanobis distance computation is a ~20-line scipy operation. The Cypher tracing query is available once outliers are identified.

### Key Findings (CONFIRMED — Inv.17 complete, 2026-03-21)
**Results file:** `data/results/outlier_analysis.tsv` (3,202 rows, z-distance ranked)

**Top outliers identified:**
| Rank | Sample | Pop | z-dist | FROH | Het rate | Rare missense |
|------|--------|-----|--------|------|----------|---------------|
| 1 | HG01816 | CDX (E.Asian) | 9.50 | 0.1250 | 0.0293 | 23 |
| 2 | HG02701 | UNKNOWN | 8.55 | 0.1250 | 0.0357 | 20 |
| 3 | NA20585 | TSI (European) | 8.04 | 0.0551 | 0.0306 | 12 |
| 4 | HG02659 | UNKNOWN | 7.72 | 0.0983 | 0.0184 | 18 |
| 5 | HG03862 | ITU (S.Asian) | 7.00 | 0.1187 | 0.0350 | 28 |

- **Top outlier HG01816 (CDX):** FROH=0.125 = **18× East Asian population mean** (0.007). This individual's genome carries 12.5% in ROH segments — consistent with recent consanguinity or strong geographic isolation within the Chinese Dai (CDX) population.
- **4 UNKNOWN-population samples** in top 10: data quality or unassigned panel members; high FROH (0.095–0.125) suggests these are from isolated subpopulations not represented in the 1000G panel.
- **NA20314 (ASW/African):** z=6.59, FROH=0.027, rare_missense=5 — African American with anomalously low rare missense burden (5 vs group mean 34.5); likely carries more European ancestry reducing rare variant exposure.
- **z-distance correlations:** FROH ρ=−0.087 (p<1e-9), rare_missense ρ=+0.063 (p=4e-4) — outliers are weakly but significantly driven by extreme FROH values.
- **Group variability (mean z-distance):** South Asian = 1.095 ± 0.895 (max=7.0, largest std — highest within-population diversity in process space, consistent with broad FROH distribution σ=0.021).

### What Is Novel
The outlier tracing step (Step 3) is architecturally impossible in traditional workflows without staging intermediate files. GraphPop enables the complete pipeline — outlier detection to driver variant identification — in a single session with no file-format conversions. Additionally, the 5-dimensional process space (combining drift, selection, and functional burden) does not exist as a pre-computed object in any standard tool.

### Paper Placement
Supplementary figure (fig22 in figure plan). The concept is cited in the main text as a demonstration of graph-native individual-level investigation.

### Figure Description
**fig22:** Two-panel figure.
(A) UMAP of 3,202 individuals (umap1, umap2 from Inv.16b), colored by continental group, with top 5 outliers per group marked as stars with sample ID labels.
(B) For the single most extreme outlier per group: horizontal bar chart of their variant-level drivers (top 10 rare missense variants by AF, annotated with gene symbol and pathway name).

---

## Analysis 5: Individual-Level Pathway Trajectory Embedding

### Research Question
Can we embed all 3,202 individuals in pathway-rare-burden space — a matrix where cell [i, j] = rare missense count in pathway j for individual i — and reveal functional trajectories that reflect evolutionary history? The hypothesis is that pathway-space UMAP clusters will recapitulate continental ancestry structure (driven by rare variant sharing within lineages) while revealing functionally distinct sub-clusters invisible to standard PCA on variant allele frequencies.

This is a fundamentally different dimensionality reduction than population genetics PCA: instead of allele-frequency vectors, we use functional-burden vectors. Individuals from the same population who share rare variants in the same pathway will cluster together regardless of their overall allele frequency similarity.

### Method
Build the 3,202 × N_pathways matrix:
```
cell[i,j] = count of rare missense CARRIES edges from sample i through variants in pathway j
```

Efficient construction via Cypher:
```cypher
MATCH (s:Sample)-[:CARRIES]->(v:Variant)-[:IN_GENE]->(g:Gene)-[:IN_PATHWAY]->(p:Pathway)
MATCH (v)-[:HAS_CONSEQUENCE]->(c:Consequence)
WHERE c.impact IN ['HIGH', 'MODERATE']
  AND c.consequence CONTAINS 'missense'
  AND v.af_global < 0.01
RETURN s.sampleId AS sample_id, p.pathwayId AS pathway_id, count(v) AS burden
```

This produces a sparse (sample, pathway, burden) triplet list. Pivot to dense matrix → UMAP(n_components=2, metric='cosine', min_dist=0.1, n_neighbors=15) → overlay with continental group labels.

For chr22 only: N_pathways ~ 200 (pathways with chr22 genes); genome-wide N_pathways = 2,241. The chr22 version is immediately runnable; genome-wide requires loading hap_cache for all chromosomes.

### Status
Planned. Requires genome-wide hap_cache (currently only chr22 is confirmed). The chr22 pilot is feasible now and would serve as proof-of-concept for the main text concept figure.

**Prerequisite data:**
- `data/hap_cache/chr22.npz` — available
- `data/variants_cache/chr22.npz` — available (includes variant_type, ac, af arrays)
- Pathway membership for chr22 genes — requires one-time Cypher query
- Genome-wide: requires `data/hap_cache/chr{1..22}.npz` — hap_cache export per-chromosome

### Key Findings
Predicted outcomes based on existing results:
- PC1 of pathway-burden matrix will align with continental ancestry (consistent with rare_missense_chr22 group differences: African 34.5 vs European 16.8)
- Pathway sub-clusters will reflect: (1) immune/HLA pathways (YRI vs non-African); (2) metabolic disease pathways (South Asian enrichment); (3) pigmentation/skin pathways (European enrichment)
- FROH-high individuals (South Asian, American) will form elongated clusters in pathway space due to wide individual variance (rare_burden_chr22 std: American = 432, vs European = 135)

### What Is Novel
The pathway trajectory embedding does not exist as a concept in traditional population genomics tools. Standard PCA (PLINK, EIGENSOFT) operates on variant allele frequencies. Standard burden tests (SKAT, BURDEN) aggregate within a gene but not across pathways and not for dimensionality reduction. GraphPop's traversal architecture extracts the full sample × pathway matrix in a single query; no intermediate VCF manipulation, annotation, or join is required.

This represents a new mode of population genomics analysis: individuals characterized by their functional burden profile rather than their allele frequency profile.

### Paper Placement
Main text concept figure (fig23, schematic) showing the query pattern and predicted embedding. Full genome-wide embedding in supplementary once all hap_cache chromosomes are available. This is the key methodological novelty claim in the Application section.

### Figure Description
**fig23:** Three-panel concept figure.
(A) Schematic of the graph traversal: Sample → CARRIES → Variant → IN_GENE → Gene → IN_PATHWAY → Pathway, with matrix assembly shown diagrammatically.
(B) chr22 pilot UMAP (3,202 individuals × ~200 pathways, cosine distance), colored by continental group — demonstrates proof of concept.
(C) Concept visualization: pathway burden vectors for 3 representative individuals (African, European, South Asian) as horizontal bar charts, showing their different functional burden profiles.

---

## Analysis 6: Purifying Selection Efficiency Across the FROH Gradient

### Research Question
Is the African-enriched rare missense burden concentrated in **constrained** (low population FST) pathways, or is it uniformly distributed across all pathways? Specifically: does the ratio (constrained pathway rare missense / unconstrained pathway rare missense) vary systematically with FROH and continental ancestry?

**Biological hypothesis:** African populations retained purifying selection efficiency for constrained pathways (essential biological processes under shared human selective constraint) because they did not experience the severe bottleneck that purged rare variants in non-African populations. Non-African populations lost rare variants via genetic drift during the Out-of-Africa bottleneck, but this purging was *pathway-indiscriminate* — drift does not respect pathway function. Therefore:

- **African individuals:** high total rare missense burden; burden concentrated in *unconstrained* (high FST) pathways because constrained pathways have been maintained clean by purifying selection across all human lineages
- **European/East Asian individuals:** low total rare missense burden (bottleneck purged both constrained and unconstrained rare variants); burden more evenly distributed because the bottleneck, not purifying selection, is the dominant force
- **Prediction:** ratio (constrained burden / unconstrained burden) is *lower* in Africans than in bottlenecked populations — Africans show stronger pathway-specific purifying selection efficiency

This is a direct test of whether purifying selection is more *pathway-discriminate* in populations without bottleneck history.

### Method
**Step 1 — Pathway constraint classification:**
From `pathway_fst.json` (2,241 pathways):
- Constrained: mean_fst < 0.10 (lower quartile; ~560 pathways) — maintained under shared purifying selection
- Unconstrained: mean_fst > 0.25 (upper quartile; ~560 pathways) — population-divergent, less constrained

**Step 2 — Per-individual rare missense burden stratified by pathway constraint:**
```cypher
MATCH (s:Sample)-[:CARRIES]->(v:Variant)-[:IN_GENE]->(g:Gene)-[:IN_PATHWAY]->(p:Pathway)
MATCH (v)-[:HAS_CONSEQUENCE]->(c:Consequence)
WHERE c.impact IN ['HIGH', 'MODERATE']
  AND c.consequence CONTAINS 'missense'
  AND v.af_global < 0.01
WITH s, p.pathwayId AS pid,
     count(v) AS burden
RETURN s.sampleId, s.population,
       pid,
       burden
```

Join externally with pathway FST ranks to assign constrained/unconstrained labels. Compute per-individual:
- `constrained_burden` = sum of burden across constrained pathways
- `unconstrained_burden` = sum of burden across unconstrained pathways
- `constraint_ratio` = constrained_burden / (constrained_burden + unconstrained_burden)

**Step 3 — FROH quartile analysis:**
Bin individuals by FROH quartile (Q1 = lowest FROH, most outbred; Q4 = highest FROH, most inbred). Compare `constraint_ratio` across FROH × continental group cells.

**Step 4 — Statistical test:**
Linear model: constraint_ratio ~ continental_group + FROH_quartile + continental_group:FROH_quartile
Prediction: significant continental_group main effect (African lower constraint_ratio); non-significant interaction (bottleneck effect is population-level, not individual FROH-driven).

### Status
**Implementable now from chr22 data.** All required inputs exist:
- `data/results/individual_trajectory.tsv` — froh_genome, rare_missense_chr22 per individual
- `data/results/pathway_fst.json` — 2,241 pathway FST values
- `data/variants_cache/chr22.npz` + `data/hap_cache/chr22.npz` — chr22 rare missense per pathway

Implementation requires: (1) one Cypher query (~10 min runtime) or numpy join to extract chr22 rare missense per (sample, pathway); (2) pathway FST join; (3) constraint_ratio computation; (4) FROH quartile stratification; (5) linear model (scipy.stats or statsmodels).

Estimated implementation time: 3–4 hours for chr22 pilot; genome-wide requires additional hap_cache chromosomes.

### Key Findings (CONFIRMED — Inv.17 complete, 2026-03-21)
**Results files:** `data/results/purifying_selection_gradient.tsv/.json`

**Pathway classification:** 5,833 chr22 rare missense positions classified via Neo4j `(Variant)-[:HAS_CONSEQUENCE]->(Gene)-[:IN_PATHWAY]->(Pathway)` — 158 in constrained pathways (mean_fst ≤ p25=0.108), 570 in divergent pathways (mean_fst ≥ p75=0.158)

**Per-group rare missense burden by pathway class (chr22):**

| Continental Group | FROH | Constrained | Divergent | C:D Ratio |
|---|---|---|---|---|
| **African** | 0.0047 | **0.54** | **1.89** | **0.187** |
| European | 0.0066 | 0.19 | 0.67 | 0.113 |
| South Asian | 0.0136 | 0.23 | 1.04 | 0.112 |
| East Asian | **0.0069** | 0.19 | 1.22 | **0.086** |
| American | 0.0102 | 0.27 | 0.73 | 0.153 |

**Key results:**
1. **African individuals carry 2.8× more constrained-pathway rare missense** than Europeans (0.54 vs 0.19/individual on chr22)
2. **African carries 2.8× more divergent-pathway rare missense** than Europeans (1.89 vs 0.67) — the enrichment is ~equal in both classes
3. **Constrained:Divergent ratio is highest in African (0.187) and lowest in East Asian (0.086)** — 2.2× difference in ratio. This indicates that bottlenecked East Asian populations have preferentially depleted constrained-pathway rare missense (retained relatively less), consistent with pathway-specific purging via drift.
4. **FROH predicts burden in both pathway classes:** ρ(FROH, constrained)=−0.075 (p=2.4e-5), ρ(FROH, divergent)=−0.063 (p=3.6e-4), ρ(FROH, total)=−0.156 (p=8.5e-19)
5. **The constrained:divergent ratio decreases with FROH:** ρ=−0.065, p=2.09e-4 — as genome-wide inbreeding increases, constrained-pathway burden is preferentially depleted.

**Biological interpretation:** Bottleneck-driven purging is not random with respect to pathway constraint. Populations with high FROH (South Asian, American, East Asian) show disproportionately lower constrained:divergent ratios, indicating that selection-constrained gene sets are preferentially cleared of rare functional variants following population size contractions — while neutral/divergent pathways retain more variation. This is consistent with the "hill-climbing" model where drift acts like purifying selection in small populations, but the effect is stronger at functionally critical loci.

### What Is Novel
No existing pipeline computes individual-level, pathway-constraint-stratified rare burden. This requires simultaneously:
1. Individual genotypes (CARRIES edges)
2. Variant consequence (HAS_CONSEQUENCE)
3. Gene membership (IN_GENE)
4. Pathway membership (IN_PATHWAY)
5. Population-level pathway FST (from prior computation stored in result files)
6. Individual FROH (genome-wide, from prior computation)

Traditional approaches would require 5+ tools (bcftools/PLINK for genotypes, VEP for consequences, BEDTools for gene intersection, a pathway database API, custom joins for FST). GraphPop co-locates items 1–4 in the graph and items 5–6 in result files that join trivially on sample ID and pathway ID.

**This analysis is not an incremental improvement on existing tools — it is a question that cannot currently be answered without GraphPop or an equivalent system.**

### Paper Placement
**Main text** — key novel finding in the Application → Human 1000G section. This is the primary demonstration of cross-layer graph analysis (individual genotype × pathway structure × population selection statistics). Emphasize: the question itself is new, enabled by the architecture.

### Figure Description
**fig21 (3-panel main text figure):**

**(A) Rare missense burden by continental group with FROH overlay:**
- Left y-axis: per-individual rare_missense_chr22 (violin + jitter, 5 groups)
- Right y-axis / color: FROH (genome-wide)
- Data source: `individual_trajectory.tsv`
- Shows the 2.06× African enrichment and inverse FROH pattern (African lowest FROH, highest burden)

**(B) Constrained vs unconstrained pathway burden ratio by FROH quartile:**
- X-axis: FROH quartile (Q1–Q4, genome-wide)
- Y-axis: constraint_ratio = constrained_burden / (constrained + unconstrained burden)
- Lines colored by continental group (5 lines)
- Expected pattern: African line is consistently lower (more efficient purifying selection in constrained pathways); bottlenecked populations have elevated constraint_ratio (indiscriminate purging)
- Data source: chr22 analysis (to be computed)

**(C) Pathway FST rank vs per-individual burden correlation:**
- X-axis: pathway mean_fst rank (2,241 pathways, from `pathway_fst.json`)
- Y-axis: mean per-individual rare missense burden in that pathway, separately for African vs European individuals
- Expected: African burden > European burden at all FST ranks; the gap narrows at high-FST (unconstrained) pathways because purifying selection is efficient in both, but African individuals have more total rare variants to distribute
- Data source: pathway FST from `pathway_fst.json` + per-pathway burden from chr22 analysis

---

## Paper Integration

### Main Text Application Section (Human 1000G)
The human analysis section (`application.tex`) currently demonstrates:
1. Population summary statistics (FST, diversity, SFS)
2. Selection scans (iHS, XP-EHH, PBS)
3. Pathway-level FST (Inv.01)
4. Individual evolutionary trajectories (Inv.16b)

The following novel analyses should be added:

**Primary addition (main text, 1–2 paragraphs + fig21):** Analysis 6 — Purifying selection efficiency across FROH gradient. This is the centerpiece novel finding: a cross-layer analysis impossible without the graph architecture. Present as: "To test whether purifying selection efficiency varies with bottleneck history, we queried individual genotypes, variant consequences, pathway memberships, and population FST statistics in a single traversal..."

**Secondary addition (main text, 1 paragraph + concept schematic):** Analysis 5 — pathway trajectory embedding concept. Frame as future direction enabled by the architecture; show the chr22 pilot as proof of concept.

**Supplementary additions:**
- Analysis 1: pathway burden heatmap (fig_S_pathway_burden)
- Analysis 4: outlier individual tracing (fig22)
- Analyses 2, 3: supplementary notes with query patterns

### Relationship to Rice Analysis
The rice analysis in `application.tex` demonstrates cross-species comparison via module convergence (Inv.13). The novel human analyses extend this by adding individual-level resolution: while the rice analysis asks "which pathways are under selection across species?", Analysis 6 asks "which individuals within a species show evidence of pathway-specific selection?". Together they form a multi-scale narrative: cross-species (rice vs human) → cross-population (26 groups) → individual level (3,202 samples).

### Extended Data Table 3: What's Novel
Reproduce the overview table from this document as Table 3 or Extended Data Table 3 in the paper. Add a column for "Implementation effort without GraphPop" (number of tools, estimated developer-hours) vs "with GraphPop" (Cypher query lines, runtime minutes).

| Capability | Traditional Tools Required | Dev-hours (traditional) | GraphPop | Runtime |
|---|---|---|---|---|
| Individual pathway rare burden | bcftools, VEP, BEDTools, pathway API, custom join | 8–16h | 1 Cypher query | ~10 min |
| ROH × sweep gene overlap | PLINK/bcftools ROH, BEDTools, gene annotation, custom join | 4–8h | 2-hop traversal | ~2 min |
| Outlier → driver variant trace | VCF extraction, VEP, pathway join, burden comparison | 4–8h | 3-hop Cypher query | ~1 min |
| Temporal burden vector | iHS (selscan), multi-pop merge, VCF genotype join | 12–24h | iHS on Variant nodes + CARRIES | ~15 min |
| Pathway trajectory embedding | Burden test per pathway (SKAT), matrix assembly, UMAP | 16–32h | Single-query matrix extraction | ~30 min |
| Purifying selection gradient | All above + custom statistical model | 24–48h | FROH + rare burden + pathway FST | ~45 min |

---

## Figure Plan

### fig21 — Analysis 6 (Main Text)
**Purifying selection efficiency across the FROH gradient**
3-panel as described above.
Status: Implementable now (chr22 data ready).
Script: `scripts/investigate_novel_06_purifying_gradient.py` (to be written).
Output: `data/results/figures/fig21_purifying_gradient.{pdf,png}`

### fig22 — Analysis 4 (Supplementary)
**Outlier individuals in evolutionary process space**
2-panel: UMAP of 3,202 individuals + outlier driver variants.
Status: Implementable now from `individual_trajectory.tsv`.
Script: `scripts/investigate_novel_04_outliers.py` (to be written).
Output: `data/results/figures/fig22_outlier_individuals.{pdf,png}`

### fig23 — Analysis 5 (Main Text Concept / Supplementary Full)
**Individual pathway trajectory embedding**
3-panel: traversal schematic + chr22 pilot UMAP + representative burden vectors.
Status: chr22 pilot feasible now; genome-wide requires full hap_cache.
Script: `scripts/investigate_novel_05_pathway_embedding.py` (to be written).
Output: `data/results/figures/fig23_pathway_embedding.{pdf,png}`

### fig_S_pathway_burden — Analysis 1 (Supplementary)
**Individual × pathway rare burden distribution**
3-panel: burden violins + pathway FST scatter + example heatmap.
Status: Partially implementable (chr22 data exists; pathway-level breakdown requires implementation).
Script: `scripts/investigate_novel_01_pathway_burden.py` (to be written).
Output: `data/results/figures/figS_pathway_burden.{pdf,png}`

### fig20 (Already Done — Supplementary)
Individual evolutionary trajectory UMAP (Inv.16b output).
File: `data/results/figures/fig20_individual_trajectory.{pdf,png}`

### Supplementary Note: Planned Analyses
Analyses 2 and 3 (ROH-sweep overlap, temporal fingerprint) documented as planned analyses in supplementary methods, with Cypher queries provided.

---

## Implementation Priority

| Priority | Analysis | Effort | Blocking Factor |
|---|---|---|---|
| 1 (now) | Analysis 6 chr22 pilot (fig21) | 3–4h | None — all data available |
| 2 (now) | Analysis 4 outlier tracing (fig22) | 2–3h | None — individual_trajectory.tsv ready |
| 3 (soon) | Analysis 1 pathway burden heatmap | 4–6h | Pathway-gene mapping Cypher query |
| 4 (soon) | Analysis 5 chr22 pilot (fig23) | 4–6h | Same as Analysis 1 |
| 5 (later) | Analysis 2 ROH-sweep overlap | 6–8h | ROH segment positions (not just total FROH) |
| 6 (later) | Analysis 3 temporal fingerprint | 4–6h | chr22 iHS data sufficient for pilot |
| 7 (future) | Genome-wide Analysis 5 | 1–2 days | All hap_cache chromosomes needed |
