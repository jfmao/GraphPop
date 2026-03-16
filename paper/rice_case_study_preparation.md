# Rice 3K Case Study — Preparation for GraphPop Platform Paper

> **Organizing principle**: Every result is presented through the lens of *what GraphPop uniquely makes possible*. Standard population genetics results (pi, Fst, tree) establish correctness by matching published values; the novel analyses that only a graph-native platform can deliver are the core contribution. Biology validates the tool; the tool's unique capabilities justify the platform.

---

## Part A: Main Text Section — "Application to the Rice 3K Genome Dataset"

### A.0 Dataset and Graph Construction

- **Source**: 3,000 Rice Genomes Project (Wang et al. 2018), ~29.6M SNPs, 3,024 accessions, 12 subpopulations, 12 chromosomes
- **Graph**: Variant nodes carry packed genotype arrays (2 bits/sample), allele-count arrays per population, and ancestral allele annotations. Gene nodes are connected via HAS_CONSEQUENCE edges (VEP), Pathway/GOTerm nodes via further edges. Total database: ~46 GB.
- **Import**: Single pipeline run (VCF → CSV → neo4j-admin import), species-agnostic. All downstream analyses are pure Cypher procedure calls — no file export, no format conversion, no external tool invocation.

### A.1 Correctness Validation: Standard Statistics Match Published Results

*Purpose: establish that GraphPop produces correct results before presenting novel capabilities.*

GraphPop's 12 stored procedures were applied to all 12 rice subpopulations across all 12 chromosomes. Standard summary statistics reproduce the established biology of Asian rice:

- **Population structure**: UPGMA tree from pairwise W&C Fst recovers the indica–japonica split and all 9 subpopulations of Wang et al. 2018, with Fst=0.71 between GJ-tmp and XI-1A (the most divergent pair) and Fst=0.06 between XI-adm and XI-3 (the most closely related).
- **Diversity gradient**: Indica (XI) 2–3x more diverse than japonica (GJ); GJ-tmp is the least diverse (pi=0.019) with strongest bottleneck signatures (Tajima's D = -2.41, F_IS = 0.72, FROH = 0.080) — matching the known temperate japonica demographic history.
- **Selection signals**: XP-EHH and Garud's H recover known domestication genes (GW5, Hd1, PROG1, Wx, OsC1, DRO1) at expected genomic positions.
- **Fay & Wu's H**: Indica populations show more negative H than japonica (XI-1A: H = -0.081 vs GJ-tmp: H = -0.044), indicating stronger/more frequent sweeps — consistent with indica's larger effective population where selection is more efficient.

These results confirm that GraphPop matches established tools (validated to <0.000001% error against scikit-allel on 1000 Genomes chr22) while operating entirely within the graph database. However, the real power of the platform lies in what follows.

### A.2 GraphPop-Unique Capability 1: Annotation-Conditioned Population Statistics

**What is new**: GraphPop computes any population statistic restricted to variants matching functional annotation criteria — consequence type, gene, pathway — in a single procedure call. This exploits the graph topology: `Variant→HAS_CONSEQUENCE→Gene` edges are traversed at query time to filter the variant set before computing statistics. No classical tool offers this; the equivalent requires a multi-step pipeline: (1) extract variant coordinates, (2) intersect with annotation BED via bedtools, (3) filter VCF, (4) run statistics tool, (5) repeat for each condition.

**Demonstration: piN/piS across all populations reveals the cost of domestication**

By calling `graphpop.diversity(chr, start, end, pop, {consequence: 'missense_variant'})` and the same with `'synonymous_variant'`, we computed piN/piS ratios for all 12 populations × 12 chromosomes (288 procedure calls). Result:

| Population | piN/piS | Population | piN/piS |
|-----------|---------|-----------|---------|
| admix | 1.146 | cA-Aus | 1.076 |
| GJ-adm | 1.145 | XI-2 | 1.061 |
| GJ-trp | 1.107 | XI-3 | 1.058 |
| XI-adm | 1.079 | XI-1B | 1.054 |
| cB-Bas | 1.078 | XI-1A | 1.053 |
| GJ-sbtrp | 1.077 | GJ-tmp | 1.018 |

All populations show piN/piS > 1.0 — genome-wide relaxation of purifying selection, the hallmark "cost of domestication" (Lu et al. 2006). Admixed populations show the highest ratios (accumulating slightly deleterious variants from multiple bottlenecked lineages), while GJ-tmp (temperate japonica) is paradoxically the most constrained despite having the strongest bottleneck, possibly reflecting strong artificial selection during temperate adaptation.

**Why this matters for GraphPop**: In classical pipelines, computing piN/piS for 12 populations × 12 chromosomes would require 288 separate VCF-filtering + statistics runs. In GraphPop, each is a one-line procedure call with an `{consequence: '...'}` option. The annotation filter is not a preprocessing step — it is an integral part of the query, resolved by graph traversal in real time.

**Demonstration: Annotation-conditioned sliding-window Fst detects adaptive protein evolution**

By running `graphpop.genome_scan()` with `{consequence: 'missense_variant'}` and `{consequence: 'synonymous_variant'}` separately, we produced sliding-window Fst profiles restricted to each variant class. No existing tool can do this.

| Population pair | Fst (missense) | Fst (synonymous) | Ratio | Interpretation |
|----------------|---------------|-----------------|-------|---------------|
| GJ-tmp vs XI-1A | 0.604 | 0.579 | 1.044 | Adaptive protein divergence |
| GJ-tmp vs GJ-trp | 0.252 | 0.228 | 1.106 | Adaptive protein divergence |
| XI-1A vs cA-Aus | 0.366 | 0.356 | 1.028 | Adaptive protein divergence |

All three pairs show missense Fst exceeding synonymous Fst (ratio > 1.0), indicating that protein-coding changes are more differentiated than neutral sites — a signature of adaptive evolution at the protein level. The strongest signal (1.106) occurs between GJ-tmp and GJ-trp, the temperate–tropical japonica split that occurred ~4.3–5.5 kya (Alam et al. 2025), suggesting rapid adaptive protein divergence during ecological specialization.

The top missense-Fst outlier windows pinpoint specific adaptive divergence hotspots (e.g., Chr1:7.45 Mb, Fst=0.966 for XI-1A vs cA-Aus) — candidate loci for adaptive protein evolution that would be invisible to tools lacking annotation conditioning.

### A.3 GraphPop-Unique Capability 2: Seamless Cross-Layer Queries (Genotype → Statistic → Gene → Pathway)

**What is new**: In a graph database, the path from a genomic variant to its functional annotation is a direct edge traversal, not a coordinate join. After computing any selection statistic (XP-EHH, PBS, iHS, Garud's H), GraphPop can immediately query the genes at peak regions via the existing `Variant→HAS_CONSEQUENCE→Gene` edges. This eliminates the "last mile" problem in classical pipelines: the gap between identifying a significant window and knowing which genes and pathways are affected.

**Demonstration: Known domestication genes recovered at selection peaks**

| Gene | Function | Statistic | Value | Comparison | Recovery method |
|------|----------|----------|-------|------------|----------------|
| GW5 | Grain width | XP-EHH | 4.91 | XI-1A vs cA-Aus | Graph traversal |
| Hd1 | Heading date | XP-EHH | -4.15 | GJ-trp vs XI-1A | Graph traversal |
| PROG1 | Prostrate growth | XP-EHH | 4.53 | GJ-tmp vs GJ-trp | Graph traversal |
| Hd1 | Heading date | XP-EHH | 4.28 | GJ-tmp vs GJ-trp | Graph traversal |
| GW5 | Grain width | XP-EHH | -3.56 | GJ-trp vs XI-1A | Graph traversal |
| OsC1 | Hull color | Garud's H12 | 0.250 | GJ-tmp (hard sweep) | Graph traversal |
| Wx | Amylose (waxy) | Garud's H12 | 0.138 | GJ-adm (soft sweep) | Graph traversal |
| DRO1 | Deep rooting | Garud's H12 | 0.139 | GJ-tmp (soft sweep) | Graph traversal |

These are well-characterized rice domestication loci (Huang et al. 2012; Alam et al. 2025). Their recovery is expected and validates the procedures. The key point for GraphPop is *how* they were found: not by exporting coordinates and running bedtools intersect, but by a graph query that follows edges from significant variants to their annotated genes. This is zero additional steps for the user.

**Demonstration: Multi-statistic convergence at a single locus**

Several loci show convergent signals across multiple statistics — detectable because all statistics reside in the same graph:

- **Chr6:9.0–9.6 Mb** (near Hd1): XP-EHH peak in GJ-trp vs XI-1A (-4.15), PBS peak in GJ-trp (1.10), Fay & Wu's H window in GJ-trp (-0.492) — all accessible by querying the same Variant nodes with different procedures, no file merging required.
- **Chr5:5.3 Mb** (near GW5): XP-EHH peaks in two independent comparisons (XI-1A vs cA-Aus: 4.91; GJ-trp vs XI-1A: -3.56) — cross-comparison consistency verified by querying the same graph.

**Why this matters**: Classical pipelines require merging outputs from separate tools (scikit-allel for EHH, PLINK for ROH, custom scripts for Garud's H, bedtools for annotation) with coordinate-based joins that are error-prone and labor-intensive. In GraphPop, all statistics are properties on the same nodes and edges, and cross-referencing is a Cypher query.

### A.4 GraphPop-Unique Capability 3: Integrated Multi-Statistic Characterization from a Single Database

**What is new**: GraphPop houses genotypes, annotations, and computed statistics in one graph. A complete population-genetic characterization — diversity, divergence, SFS, selection scans (EHH, Garud's H, Fay & Wu's H, PBS), ROH — is executed as a series of procedure calls with no intermediate files, no format conversions, and no data export. Results are written back to graph nodes (GenomicWindow, Population, Variant), immediately available for downstream Cypher queries.

**Rice 3K characterization — scope**:

| Analysis | Procedure | Scale | Key finding |
|----------|----------|-------|-------------|
| Diversity (pi, theta_W, D, F_IS) | graphpop.diversity | 12 pops × 12 chr | GJ-tmp least diverse (pi=0.019) |
| **piN/piS** (annotation-conditioned) | graphpop.diversity | 12 pops × 12 chr × 2 | All pops piN/piS > 1.0 (**unique**) |
| Divergence (Fst, Dxy, PBS) | graphpop.divergence | 66 pairs + 6 PBS trios | GJ-tmp vs XI-1A: Fst=0.71 |
| **Missense vs synonymous Fst** | graphpop.genome_scan | 3 pairs × 12 chr × 2 | Ratio 1.03–1.11 (**unique**) |
| Unfolded SFS | graphpop.sfs | 12 pops | Polarized via ancestral alleles |
| Genome scan (windows) | graphpop.genome_scan | ~7,468 windows/pop | Sliding-window Fst, pi, D, H |
| XP-EHH | graphpop.xpehh | 6 pairs | GW5, Hd1, PROG1 recovered |
| Garud's H (sweeps) | graphpop.garud_h | 4 pops | 408 sweeps in GJ-tmp |
| nSL | graphpop.nsl | Selected pops | AF-bin standardized |
| ROH (HMM) | graphpop.roh | 12 pops (~3,024 samples) | GJ-tmp FROH=0.080 |
| Fay & Wu's H | graphpop.diversity | 12 pops × 12 chr | XI pops: stronger sweeps |
| Pop summary | graphpop.pop_summary | 12 pops | Written to Population nodes |

Rows marked **unique** represent analyses that no existing tool provides as a single operation.

The entire analysis — from import to final results — required no tool outside GraphPop and the Neo4j database. No VCF re-extraction, no BED intersection, no scikit-allel, no PLINK, no bedtools. This is the operational advantage of a graph-native platform.

### A.5 Biological Synthesis (Brief — for Platform Paper Context)

The rice 3K results, made possible by GraphPop's integrated analysis capabilities, reveal:

1. **The cost of domestication is universal but unequal**: piN/piS > 1.0 in all 12 subpopulations, with admixed groups accumulating the most deleterious load — a finding that required annotation-conditioned diversity (GraphPop-unique).

2. **Adaptive protein evolution accompanies ecological divergence**: Missense Fst exceeds synonymous Fst across all major population splits, with the strongest signal (1.106) at the temperate–tropical japonica boundary — detected via annotation-conditioned genome scans (GraphPop-unique).

3. **Temperate japonica bears the strongest bottleneck signature**: Lowest diversity (pi=0.019), highest inbreeding (FROH=0.080), most selective sweeps (408 Garud's H windows), and strongly negative Tajima's D (-2.41) — a coherent picture assembled from multiple procedures querying the same graph.

4. **Convergent selection signals validate across statistics**: Domestication genes (GW5, Hd1, PROG1, Wx, OsC1) appear as peaks in independent analyses (XP-EHH, Garud's H, PBS, Fay & Wu's H), all cross-referenced to gene annotations via graph traversal — demonstrating the power of having all data layers in one queryable structure.

---

## Part B: Supplementary Material — Full Rice 3K Analysis Report

> In each section below, we first state the GraphPop capability demonstrated, then present the full results. Sections marked with **(GraphPop-unique)** describe analyses that cannot be replicated with any single existing tool.

### Supplementary Note S1: Dataset, Graph Construction, and Import Pipeline

**Data source**: The 3K Rice Genomes Project (Wang et al. 2018) SNP dataset, comprising 3,024 accessions of Asian cultivated rice (*Oryza sativa* L.) with ~29.6 million SNPs called against the Nipponbare IRGSP 1.0 reference genome.

**Population classification**: 12 subpopulations following Wang et al. 2018 ADMIXTURE at K=9:

| Population | Full name | Geographic origin |
|-----------|-----------|-------------------|
| XI-1A | Xian/Indica 1A | East Asia |
| XI-1B | Xian/Indica 1B | Modern varieties, diverse |
| XI-2 | Xian/Indica 2 | South Asia |
| XI-3 | Xian/Indica 3 | Southeast Asia |
| XI-adm | Xian/Indica admixed | Mixed |
| GJ-tmp | Geng/Japonica temperate | East Asia (temperate) |
| GJ-trp | Geng/Japonica tropical | Southeast Asia |
| GJ-sbtrp | Geng/Japonica subtropical | Subtropical East/SE Asia |
| GJ-adm | Geng/Japonica admixed | Mixed |
| cA-Aus | circum-Aus | South Asia (Bangladesh, India) |
| cB-Bas | circum-Basmati | South Asia (aromatic varieties) |
| admix | Fully admixed | Mixed (between major groups) |

**Graph schema** (demonstrating GraphPop's multi-layer data model):

```
(:Variant {chr, pos, ref, alt, pop_ids[], ac[], an[], af[],
           gt_packed, phase_packed, ancestral_allele, is_polarized})
    -[:NEXT {distance}]-> (:Variant)     // genomic order
    -[:HAS_CONSEQUENCE {impact, biotype}]-> (:Gene {name, chr, start, end})
        -[:IN_PATHWAY]-> (:Pathway {name, source})
        -[:HAS_GO_TERM]-> (:GOTerm {go_id, name, namespace})
(:Sample {name, population, sex, packed_index})
    -[:IN_POPULATION]-> (:Population {name})
(:Chromosome {name, length})
(:GenomicWindow {chr, start, end, pi, theta_w, tajima_d, fst_wc, pbs, fay_wu_h, ...})
    -[:ON_CHROMOSOME]-> (:Chromosome)
```

This schema is the foundation of GraphPop's unique capabilities: annotation-conditioned statistics traverse `HAS_CONSEQUENCE` edges; cross-statistic gene lookups follow the same edges in reverse; all statistics reside on the same nodes, enabling integrated queries.

**Storage**: Packed genotype arrays (2 bits/sample for genotypes, 1 bit/sample for phase) on Variant nodes. Database size: ~46 GB (reduced from ~200 GB under the previous CARRIES-edge model).

### Supplementary Table S1: Per-Population Diversity Statistics

*GraphPop capability demonstrated: `graphpop.diversity` — standard summary statistics computed from packed genotype arrays on Variant nodes.*

| Rank | Population | Mean pi | Mean theta_W | Mean Tajima's D | Mean F_IS |
|------|-----------|---------|-------------|----------------|----------|
| 1 | admix | 0.07770 | 0.11278 | -1.117 | 0.616 |
| 2 | XI-adm | 0.06052 | 0.11325 | -1.493 | 0.632 |
| 3 | XI-2 | 0.05707 | 0.09338 | -1.528 | 0.632 |
| 4 | cA-Aus | 0.05488 | 0.12711 | -0.939 | 0.648 |
| 5 | XI-3 | 0.05418 | 0.08262 | -1.092 | 0.615 |
| 6 | XI-1B | 0.05064 | 0.09863 | -1.306 | 0.598 |
| 7 | XI-1A | 0.04691 | 0.08591 | -0.859 | 0.541 |
| 8 | cB-Bas | 0.04489 | 0.07592 | -2.043 | 0.663 |
| 9 | GJ-trp | 0.03239 | 0.05454 | -1.312 | 0.700 |
| 10 | GJ-adm | 0.03234 | 0.06210 | -0.988 | 0.706 |
| 11 | GJ-sbtrp | 0.02579 | 0.04686 | -0.890 | 0.673 |
| 12 | GJ-tmp | 0.01891 | 0.05244 | -2.406 | 0.720 |

Indica (XI) is 2–3x more diverse than japonica (GJ), matching Wang et al. 2018. All populations show negative Tajima's D and high F_IS, consistent with the selfing mating system of cultivated rice.

### Supplementary Table S2: piN/piS Ratios (GraphPop-unique)

*GraphPop capability demonstrated: annotation-conditioned diversity — `graphpop.diversity` with `{consequence: 'missense_variant'}` and `{consequence: 'synonymous_variant'}` options. This analysis requires 288 procedure calls (12 pops × 12 chr × 2 consequence types) but no VCF filtering, no external annotation tools, and no intermediate files.*

| Rank | Population | Mean piN/piS | Per-chromosome range |
|------|-----------|-------------|---------------------|
| 1 | admix | 1.146 | 0.952–1.302 |
| 2 | GJ-adm | 1.145 | 0.933–1.302 |
| 3 | GJ-trp | 1.107 | 0.933–1.302 |
| 4 | XI-adm | 1.079 | 0.952–1.302 |
| 5 | cB-Bas | 1.078 | 0.933–1.302 |
| 6 | GJ-sbtrp | 1.077 | 0.933–1.302 |
| 7 | cA-Aus | 1.076 | 0.952–1.302 |
| 8 | XI-2 | 1.061 | 0.952–1.302 |
| 9 | XI-3 | 1.058 | 0.952–1.302 |
| 10 | XI-1B | 1.054 | 0.952–1.302 |
| 11 | XI-1A | 1.053 | 0.952–1.302 |
| 12 | GJ-tmp | 1.018 | 0.933–1.105 |

All populations show piN/piS > 1.0 — genome-wide relaxation of purifying selection, the "cost of domestication." GJ-tmp, despite the strongest bottleneck, has the lowest piN/piS (1.018), possibly reflecting strong artificial selection for agronomic traits during temperate adaptation that partially counteracts genetic drift.

**Classical pipeline equivalent**: To compute piN/piS for one population on one chromosome, one would need to: (1) extract missense variants from VCF using SnpSift or bcftools with a GFF annotation, (2) compute pi on the filtered set with scikit-allel or vcftools, (3) repeat for synonymous variants, (4) compute ratio. For 288 combinations, this becomes a substantial scripting effort. GraphPop reduces each to a single procedure call with one changed parameter.

### Supplementary Table S3: Pairwise Population Divergence (W&C Fst)

*GraphPop capability demonstrated: `graphpop.divergence` — pairwise Fst for all 66 population pairs, computed from allele-count arrays on Variant nodes.*

**Top 10 most divergent pairs:**

| Pair | Mean Fst |
|------|---------|
| GJ-tmp vs XI-1A | 0.711 |
| GJ-tmp vs XI-1B | 0.694 |
| GJ-tmp vs cA-Aus | 0.678 |
| GJ-tmp vs XI-2 | 0.663 |
| XI-3 vs GJ-tmp | 0.652 |
| XI-1A vs GJ-sbtrp | 0.636 |
| GJ-trp vs XI-1A | 0.622 |
| XI-1B vs GJ-sbtrp | 0.615 |
| GJ-trp vs XI-1B | 0.603 |
| XI-1A vs GJ-adm | 0.596 |

**Top 5 most closely related pairs:**

| Pair | Mean Fst |
|------|---------|
| XI-adm vs XI-3 | 0.061 |
| XI-adm vs XI-1B | 0.062 |
| XI-adm vs XI-2 | 0.067 |
| XI-adm vs XI-1A | 0.096 |
| XI-3 vs XI-2 | 0.125 |

### Supplementary Table S4: Population Tree (UPGMA on Fst)

*GraphPop capability demonstrated: complete Fst matrix from `graphpop.divergence`, used to construct population phylogeny.*

Newick: `(((GJ-adm:0.0626,GJ-tmp:0.0626):0.0821,(GJ-sbtrp:0.1286,GJ-trp:0.1286):0.0161):0.1282,((XI-1A:0.0815,(XI-1B:0.0617,(XI-2:0.0481,(XI-3:0.0303,XI-adm:0.0303):0.0178):0.0136):0.0198):0.1004,(cA-Aus:0.1793,(admix:0.1047,cB-Bas:0.1047):0.0746):0.0026):0.0911);`

Two major clades: japonica (GJ) and indica (XI) + aus/basmati/admix, matching the established taxonomy.

### Supplementary Table S5: Fay & Wu's H

*GraphPop capability demonstrated: `graphpop.diversity` with Fay & Wu's H computed from ancestral-allele-polarized SFS, requiring `ancestral_allele` properties stored on Variant nodes during import.*

| Rank | Population | Mean H | Mean H_norm |
|------|-----------|--------|------------|
| 1 | XI-1A | -0.081 | -1.139 |
| 2 | XI-adm | -0.078 | -1.506 |
| 3 | XI-1B | -0.075 | -1.083 |
| 4 | XI-3 | -0.074 | -1.407 |
| 5 | XI-2 | -0.068 | -1.099 |
| 6 | cA-Aus | -0.061 | -0.810 |
| 7 | GJ-trp | -0.048 | -0.953 |
| 8 | GJ-tmp | -0.044 | -0.868 |
| 9 | cB-Bas | -0.038 | -0.416 |
| 10 | GJ-adm | -0.028 | -0.457 |
| 11 | GJ-sbtrp | -0.028 | -0.430 |
| 12 | admix | -0.002 | -0.044 |

Indica populations show stronger sweep signals (more negative H) than japonica, indicating more frequent recent sweeps in the larger, more diverse indica gene pool.

### Supplementary Table S6: Garud's H Selective Sweep Summary

*GraphPop capability demonstrated: `graphpop.garud_h` — sliding-window haplotype homozygosity, distinguishing hard from soft sweeps, with immediate gene annotation via graph traversal.*

| Population | Sweeps (H12>0.10) | Hard | Soft |
|-----------|-------------------|------|------|
| GJ-tmp | 408 | 135 | 273 |
| GJ-adm | 96 | 19 | 77 |
| GJ-sbtrp | 13 | 5 | 8 |
| GJ-trp | 10 | 3 | 7 |
| XI-1B | 9 | 6 | 3 |
| XI-1A | 7 | 5 | 2 |
| admix | 7 | 3 | 4 |
| cB-Bas | 7 | 2 | 5 |
| XI-adm | 6 | 4 | 2 |
| XI-3 | 6 | 4 | 2 |
| XI-2 | 6 | 4 | 2 |
| cA-Aus | 5 | 3 | 2 |

GJ-tmp has 40x more sweep windows than typical indica populations, with 2:1 soft-to-hard sweep ratio. Known domestication genes at sweep peaks (identified via graph traversal):
- **OsC1** (hull color): GJ-tmp Chr6:4.95–5.10 Mb, H12=0.250 (hard sweep)
- **Wx** (amylose/waxy): GJ-adm Chr6:1.45–1.60 Mb, H12=0.138 (soft sweep)
- **DRO1** (deep rooting): GJ-tmp Chr9:22.04–22.14 Mb, H12=0.139 (soft sweep)

### Supplementary Table S7: Runs of Homozygosity (HMM Method)

*GraphPop capability demonstrated: `graphpop.roh` — per-sample HMM-based ROH detection (Narasimhan et al. 2016), parallelized across all samples within the database.*

| Rank | Population | Total ROH (Mb) | Mean FROH |
|------|-----------|---------------|-----------|
| 1 | GJ-tmp | 3,371.7 | 0.080 |
| 2 | GJ-adm | 81.3 | 0.046 |
| 3 | GJ-trp | 78.2 | 0.036 |
| 4 | admix | 5.3 | 0.037 |
| 5 | GJ-sbtrp | 22.0 | 0.035 |
| 6 | cB-Bas | 2.1 | 0.029 |
| 7 | XI-1A | 1.0 | 0.028 |
| 8 | cA-Aus | 1.0 | 0.024 |
| 9–12 | XI-adm, XI-3, XI-2, XI-1B | 0.0 | 0.000 |

GJ-tmp: 153/287 samples show ROH on Chr1, 175/287 on Chr2.

### Supplementary Table S8: XP-EHH Selection Peaks (6 population pairs)

*GraphPop capability demonstrated: `graphpop.xpehh` — cross-population EHH, with results written as properties on Variant nodes. Gene annotations at peaks retrieved via graph traversal (no bedtools intersection).*

| Pair | Total hits | Peaks | Max |XP-EHH| | Known genes at peaks |
|------|-----------|-------|----------------|---------------------|
| cA-Aus vs cB-Bas | 29,661 | 398 | 5.33 | — |
| XI-1A vs cA-Aus | 99,624 | 621 | 7.95 | GW5 (Chr5:5.27 Mb) |
| GJ-tmp vs cA-Aus | 17,994 | 248 | 6.63 | — |
| GJ-tmp vs XI-1A | 15,242 | 258 | 5.48 | — |
| GJ-trp vs XI-1A | 18,477 | 348 | 5.48 | Hd1 (Chr6:8.95 Mb), GW5 |
| GJ-tmp vs GJ-trp | 42,018 | 339 | 5.51 | PROG1 (Chr7:3.18 Mb), Hd1 |

### Supplementary Table S9: Population Branch Statistic (PBS)

*GraphPop capability demonstrated: `graphpop.divergence` with `{pop3: outgroup}` — three-population PBS (Yi et al. 2010) computed from W&C Fst.*

| Focal | vs | Outgroup | Mean PBS | Max PBS | Top locus |
|-------|-----|---------|---------|---------|-----------|
| GJ-tmp | GJ-trp | XI-1A | 0.305 | 2.760 | Chr5:28.2 Mb |
| cA-Aus | XI-1A | GJ-tmp | 0.213 | 2.755 | Chr1:8.3 Mb |
| XI-1A | XI-1B | GJ-tmp | 0.122 | 1.730 | Chr4:24.2 Mb |
| cB-Bas | cA-Aus | GJ-tmp | 0.098 | 1.590 | Chr6:6.2 Mb |
| XI-2 | XI-3 | GJ-tmp | 0.073 | 0.923 | Chr7:1.9 Mb |
| GJ-trp | GJ-tmp | XI-1A | 0.051 | 1.369 | Chr7:24.5 Mb |

### Supplementary Table S10: Annotation-Conditioned Divergence (GraphPop-unique)

*GraphPop capability demonstrated: `graphpop.genome_scan` with `{consequence: 'missense_variant'}` — sliding-window Fst restricted to a functional annotation class. **This analysis has no equivalent in any existing population genetics tool.***

| Pair | Fst (missense) | Fst (synonymous) | Ratio | Interpretation |
|------|---------------|-----------------|-------|---------------|
| GJ-tmp vs XI-1A | 0.604 | 0.579 | 1.044 | Adaptive protein divergence |
| GJ-tmp vs GJ-trp | 0.252 | 0.228 | 1.106 | Adaptive protein divergence |
| XI-1A vs cA-Aus | 0.366 | 0.356 | 1.028 | Adaptive protein divergence |

**Top missense-Fst outlier windows** (candidate adaptive protein evolution hotspots):

| Pair | Chr | Position | Fst | #Variants |
|------|-----|----------|-----|-----------|
| XI-1A vs cA-Aus | Chr1 | 7,453,527–7,553,526 | 0.966 | 449 |
| XI-1A vs cA-Aus | Chr1 | 7,403,527–7,503,526 | 0.959 | 608 |
| GJ-tmp vs XI-1A | Chr5 | 25,004,007–25,104,006 | 0.958 | 356 |
| GJ-tmp vs XI-1A | Chr4 | 35,452,105–35,552,104 | 0.955 | 57 |
| GJ-tmp vs XI-1A | Chr2 | 12,003,895–12,103,894 | 0.954 | 770 |

These loci harbor extreme missense differentiation — candidate regions for adaptive protein evolution during the indica–japonica and indica–aus divergence. They are uniquely discoverable because GraphPop can restrict genome scans to a specific functional annotation class.

### Supplementary Table S11: Unfolded Site Frequency Spectrum

*GraphPop capability demonstrated: `graphpop.sfs` with ancestral allele polarization — unfolded SFS from `ancestral_allele` properties on Variant nodes.*

| Population | Total variants | Polarized | Singletons | High-freq derived (>50%) |
|-----------|---------------|-----------|------------|------------------------|
| XI-adm | 29,629,320 | 16,650,808 | 11,987,722 | 1,880,240 |
| GJ-tmp | 29,629,003 | 16,650,727 | 20,986,070 | 889,429 |
| XI-1A | 29,566,291 | 16,639,041 | 20,363,376 | 2,030,250 |
| cA-Aus | 29,505,350 | 16,620,990 | 18,698,580 | 1,939,737 |
| GJ-trp | 29,625,440 | 16,649,776 | 19,952,358 | 1,183,032 |

### Supplementary Table S12: Derived Allele Frequency Enrichment at Selection Peaks

*GraphPop capability demonstrated: integrated cross-referencing of selection signals with ancestral allele annotations. Because both selection statistics and ancestral allele data reside on the same Variant nodes, DAF enrichment analysis requires no external tool or file join.*

**PBS peak regions:**

| Focal population | Peak DAF | Background DAF | Enrichment |
|-----------------|----------|---------------|------------|
| GJ-trp | 0.136 | 0.070 | 1.93x |
| cB-Bas | 0.157 | 0.082 | 1.93x |
| cA-Aus | 0.140 | 0.101 | 1.38x |
| XI-1A | 0.139 | 0.101 | 1.37x |
| GJ-tmp | 0.062 | 0.059 | 1.06x |

**XP-EHH peak regions:**

| Pair | Focal pop | Peak DAF | Background DAF | Enrichment |
|------|----------|----------|---------------|------------|
| GJ-trp vs XI-1A | GJ-trp | 0.120 | 0.070 | 1.70x |
| GJ-tmp vs cA-Aus | GJ-tmp | 0.097 | 0.059 | 1.64x |
| GJ-tmp vs XI-1A | GJ-tmp | 0.094 | 0.059 | 1.61x |
| GJ-tmp vs GJ-trp | GJ-tmp | 0.091 | 0.059 | 1.55x |
| cA-Aus vs cB-Bas | cA-Aus | 0.126 | 0.101 | 1.24x |

Elevated DAF at selection peaks confirms that sweeps drove derived alleles toward fixation, validating both the selection scan procedures and the ancestral allele polarization.

### Supplementary Table S13: Fay & Wu's H Genome Scan — Top Sweep Windows

*GraphPop capability demonstrated: `graphpop.genome_scan` with Fay & Wu's H per window, computed from ancestral-allele-polarized data within the graph.*

**GJ-tmp** (7,468 windows, mean H = -0.049):

| Chr | Position | H | Top gene |
|-----|----------|---|----------|
| Chr8 | 17,701,028–17,801,027 | -0.558 | LOC_Os08g29030 |
| Chr4 | 18,401,004–18,501,003 | -0.528 | LOC_Os04g30900 |
| Chr9 | 8,287,521–8,387,520 | -0.502 | LOC_Os09g14120 |

**XI-1A** (7,468 windows, mean H = -0.082):

| Chr | Position | H | Top gene |
|-----|----------|---|----------|
| Chr7 | 14,151,017–14,251,016 | -1.314 | LOC_Os07g24910 |
| Chr12 | 3,651,048–3,751,047 | -0.644 | LOC_Os12g07530 |
| Chr12 | 5,301,048–5,401,047 | -0.622 | LOC_Os12g10110 |

**cA-Aus** (7,468 windows, mean H = -0.062):

| Chr | Position | H | Top gene |
|-----|----------|---|----------|
| Chr7 | 14,151,017–14,251,016 | -1.716 | LOC_Os07g24910 |
| Chr8 | 17,701,028–17,801,027 | -0.710 | LOC_Os08g29030 |
| Chr2 | 8,300,030–8,400,029 | -0.543 | LOC_Os02g14900 |

Note: Chr7:14.15 Mb shows the strongest sweep signal in both XI-1A and cA-Aus (H = -1.31 and -1.72), suggesting a shared ancient sweep predating the indica–aus divergence.

---

## Part C: Figure and Table Suggestions for Main Text

### Priority: figures that showcase GraphPop-unique capabilities

**Figure: Rice 3K case study (multi-panel)**

| Panel | Content | GraphPop-unique? | Message |
|-------|---------|-----------------|---------|
| **A** | piN/piS bar chart across 12 populations | **Yes** — annotation-conditioned diversity | "Cost of domestication" detected in single procedure calls |
| **B** | Missense vs synonymous Fst (paired bars, 3 pairs) | **Yes** — annotation-conditioned genome scan | Adaptive protein evolution detected; no classical tool equivalent |
| **C** | XP-EHH Manhattan with gene labels at peaks | **Partially** — gene labels via graph traversal | Domestication genes recovered via graph edges, not bedtools |
| **D** | Population tree (UPGMA on Fst) | No (validates correctness) | Correct recovery of known structure |
| **E** | Diversity–FROH scatter across 12 populations | No (validates correctness) | Coherent multi-statistic picture from one database |
| **F** | Workflow comparison diagram | N/A (conceptual) | Classical pipeline (6 tools) vs GraphPop (1 database) |

**Recommended priority for main text** (space is limited):
1. **Panels A + B** are the strongest — they show results that literally cannot be obtained from existing tools. These should be the centerpiece.
2. **Panel C** (XP-EHH + gene labels) is visually compelling and demonstrates graph traversal.
3. **Panel D** (tree) is necessary for orientation but can be small.
4. **Panel F** (workflow comparison) may work better as a figure in the Methods section of the main paper rather than in the rice results section.
5. **Panel E** is nice but secondary — could move to supplementary.

### Main Text Table

**Table: GraphPop-unique analyses enabled by the graph-native architecture**

| Analysis | GraphPop call | Classical equivalent | #Steps saved |
|----------|--------------|---------------------|-------------|
| piN for one pop, one chr | `graphpop.diversity(chr, s, e, pop, {consequence: 'missense_variant'})` | bcftools view + SnpSift filter + vcftools --site-pi | 3→1 |
| piN/piS for 12 pops × 12 chr | 288 × above | 288 × (VCF filter + stats + ratio script) | 864→288 |
| Missense-Fst genome scan | `graphpop.genome_scan(chr, pop, w, s, {pop2: p2, consequence: 'missense_variant'})` | No single-tool equivalent exists | N/A→1 |
| Gene at XP-EHH peak | Result already linked via graph edges | bedtools intersect + gene annotation GFF | 2→0 |
| Cross-statistic locus query | Cypher: `MATCH (v:Variant) WHERE v.ihs_pop > 3 AND v.xpehh_pop1_pop2 > 3` | Custom script merging outputs from 2+ tools | N→1 |

---

## Key References for Rice Section

1. Wang W et al. (2018). Genomic variation in 3,010 diverse accessions of Asian cultivated rice. *Nature* 557:43–49. doi:10.1038/s41586-018-0063-9
2. Huang X et al. (2012). A map of rice genome variation reveals the origin of cultivated rice. *Nature* 490:497–501. doi:10.1038/nature11532
3. Alam O et al. (2025). Evolutionary histories of functional mutations during the domestication and spread of japonica rice in Asia. *PNAS* 122(43):e2514614122. doi:10.1073/pnas.2514614122
4. Weir BS, Cockerham CC (1984). Estimating F-statistics for the analysis of population structure. *Evolution* 38:1358–1370.
5. Tajima F (1989). Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. *Genetics* 123:585–595.
6. Fay JC, Wu CI (2000). Hitchhiking under positive Darwinian selection. *Genetics* 155:1405–1413.
7. Garud NR et al. (2015). Recent selective sweeps in North American Drosophila melanogaster show signatures of soft sweeps. *PLoS Genet* 11:e1005004.
8. Yi X et al. (2010). Sequencing of 50 human exomes reveals adaptation to high altitude. *Science* 329:75–78.
9. Ferrer-Admetlla A et al. (2014). On detecting incomplete soft or hard selective sweeps using haplotype structure. *Mol Biol Evol* 31:1275–1291.
10. Narasimhan V et al. (2016). BCFtools/RoH: a hidden Markov model approach for detecting autozygosity from next-generation sequencing data. *Bioinformatics* 32:1749–1751.
11. Wigginton JE et al. (2005). A note on exact tests of Hardy-Weinberg equilibrium. *Am J Hum Genet* 76:887–893.
12. Lu J et al. (2006). The accumulation of deleterious mutations in rice genomes: a hypothesis on the cost of domestication. *Trends Genet* 22:126–131.
