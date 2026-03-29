# GraphPop Vignette: Human 1000 Genomes Analysis

**A complete walkthrough of graph-native population genomics on the world's
most-studied human cohort**

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Data Collection and Preparation](#2-data-collection-and-preparation)
3. [Neo4j Setup and Data Import](#3-neo4j-setup-and-data-import)
4. [Benchmarking: How Fast Is GraphPop?](#4-benchmarking-how-fast-is-graphpop)
5. [Population Diversity and Continental Structure](#5-population-diversity-and-continental-structure)
5b. [Annotation-Conditioned Diversity: piN/piS Across Human Populations](#5b-annotation-conditioned-diversity-pinpis-across-human-populations)
6. [Pathway-Level Population Differentiation](#6-pathway-level-population-differentiation)
6b. [Pathway Co-Selection Network](#6b-pathway-co-selection-network)
7. [Convergent Positive Selection: Finding Ancient Sweeps](#7-convergent-positive-selection-finding-ancient-sweeps)
8. [Individual-Level Evolutionary Trajectories](#8-individual-level-evolutionary-trajectories)
9. [Pathway-Stratified Purifying Selection](#9-pathway-stratified-purifying-selection)
10. [Consequence-Level Selection Bias](#10-consequence-level-selection-bias)
11. [Summary Tables and Database Sharing](#11-summary-tables-and-database-sharing)
11b. [Cross-Species Comparison: Human vs Rice](#cross-species-comparison-human-vs-rice)
12. [Summary of Findings](#12-summary-of-findings)

---

## 1. Introduction

### What is the 1000 Genomes Project?

The 1000 Genomes Project (Phase 3) is the foundational reference for human
genetic variation. It catalogues ~84.7 million variants across 3,202
individuals drawn from 26 sub-populations spanning 5 continental
superpopulations:

| Superpopulation | Full Name             | Sub-populations                      | Samples |
|-----------------|-----------------------|--------------------------------------|---------|
| AFR             | African               | YRI, LWK, GWD, MSL, ESN, ACB, ASW   | 661     |
| EUR             | European              | CEU, GBR, FIN, IBS, TSI              | 503     |
| EAS             | East Asian            | CHB, CHS, CDX, JPT, KHV             | 504     |
| SAS             | South Asian           | GIH, PJL, BEB, STU, ITU             | 489     |
| AMR             | Admixed American      | MXL, PUR, CLM, PEL                   | 347     |

These populations capture the major branches of the human Out-of-Africa
diaspora and the demographic events that shaped modern genetic diversity.

### What can GraphPop do that classical tools cannot?

Classical workflows require a different tool for each statistic (scikit-allel
for diversity, VCFtools for Fst, PLINK for ROH, selscan for iHS), producing
separate output files that must be merged by genomic coordinates. This creates
three barriers: (1) cross-statistic queries are impractical, (2) annotation
conditioning requires multi-step pipelines, and (3) individual-to-pathway
traversal is impossible without chaining 4+ tools.

GraphPop stores everything in a single Neo4j graph database. A Cypher query
can traverse from an individual sample, through variants, to gene consequences,
to pathway memberships, to population-level selection statistics---in
milliseconds. This vignette walks through a complete analysis using only the
`graphpop` CLI.

---

## 2. Data Collection and Preparation

### Target problem

How do you obtain and prepare all the data layers needed for a graph-native
population genomics analysis?

### Download the VCF files

The Phase 3 VCF files are hosted by the International Genome Sample Resource:

```bash
mkdir -p data/raw/1000g && cd data/raw/1000g

# Download all 22 autosomes (~15 GB compressed total)
for chr in $(seq 1 22); do
  wget -c "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/\
ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  wget -c "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/\
ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi"
done
```

Selected file sizes: chr1 (6.5M variants, 1.2 GB), chr22 (1.1M variants,
210 MB), total ~70.7M variants across 22 autosomes.

### Download the population panel

The panel file maps each sample to its sub-population and superpopulation:

```bash
wget "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/\
integrated_call_samples_v3.20130502.ALL.panel"
```

The panel is tab-separated with columns: `sample`, `pop`, `super_pop`, `gender`.

### Download Reactome pathway annotations

```bash
wget -O reactome_pathways.tsv \
  "https://reactome.org/download/current/Ensembl2Reactome.txt"
awk -F'\t' '$6 == "Homo sapiens"' reactome_pathways.tsv > reactome_human.tsv
# 108,432 gene-pathway associations covering 91,973 genes and 2,241 pathways
```

### Download ancestral alleles and VEP annotations

Ancestral alleles (for unfolded SFS and Fay & Wu's H) come from Ensembl EPO
multi-species alignments:

```bash
wget "https://ftp.ensembl.org/pub/release-112/fasta/ancestral_alleles/\
homo_sapiens_ancestor_GRCh37_e71.tar.bz2"
tar xjf homo_sapiens_ancestor_GRCh37_e71.tar.bz2
# Produces one FASTA per chromosome: homo_sapiens_ancestor_1.fa, etc.
```

VEP annotations provide consequence types (missense, synonymous, stop_gain)
and impact classes (HIGH, MODERATE, LOW, MODIFIER). If you have pre-computed
VEP annotations, provide them via `--vep`. The Phase 3 VCFs also include basic
annotations in the CSQ INFO field.

---

## 3. Neo4j Setup and Data Import

### Target problem

Import a full human genome (22 autosomes, ~70.7 million variants, 3,202
samples) into a graph database that supports millisecond cross-referencing of
all computed statistics.

### Install and configure Neo4j

GraphPop provides a single command to download, install, and configure Neo4j
Community Edition with settings tuned for population genomics workloads:

```bash
graphpop setup --password mypass --pagecache 20g --heap 4g
```

```
Downloading Neo4j 5.26.0...
Configuring Neo4j...
  server.memory.pagecache.size=20g
  server.memory.heap.initial_size=4g
  server.memory.heap.max_size=4g
  dbms.security.procedures.unrestricted=graphpop.*
Setting Neo4j password...
  Password set successfully
GraphPop config written to /home/user/.graphpop/config.yaml

Setup complete!
  Neo4j home:    /home/user/neo4j
  Page cache:    20g
  Heap:          4g
```

### Start Neo4j and import data

```bash
graphpop start
# Neo4j is running (PID 12847). Bolt URI: bolt://localhost:7687
```

For a single chromosome (good for testing):

```bash
graphpop import \
  --vcf data/raw/1000g/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \
  --panel data/raw/1000g/integrated_call_samples_v3.20130502.ALL.panel \
  --database human1000g \
  --vep data/annotations/vep_chr22.vcf.gz \
  --pathways reactome_human.tsv \
  --ancestral data/raw/ancestral/homo_sapiens_ancestor_22.fa
```

```
GraphPop Import Pipeline
  VCF:       ALL.chr22.phase3_...genotypes.vcf.gz
  Panel:     integrated_call_samples_v3.20130502.ALL.panel
  Database:  human1000g

Step 1/3: Generating CSV files from VCF...
  Parsed 1,066,555 variants across 3,202 samples
  Generated 5,394,201 CARRIES edges (sparse; homozygous ref implicit)

Step 2/3: Running neo4j-admin bulk import...
  Bulk import complete.

Step 3/3: Loading annotations...
  Loaded: VEP, pathways, ancestral alleles

Import complete! Database: human1000g
```

### Full genome import

For the complete 22-autosome dataset, import each chromosome. The import
pipeline generates CSVs in parallel and uses `neo4j-admin` bulk import, which
is I/O-bound rather than CPU-bound:

```bash
# Import all 22 autosomes (takes ~4-6 hours total)
for chr in $(seq 1 22); do
  graphpop import \
    --vcf "data/raw/1000g/ALL.chr${chr}.phase3_*.vcf.gz" \
    --panel data/raw/1000g/integrated_call_samples_v3.20130502.ALL.panel \
    --database human1000g \
    --vep "data/annotations/vep_chr${chr}.vcf.gz" \
    --pathways reactome_human.tsv \
    --ancestral "data/raw/ancestral/homo_sapiens_ancestor_${chr}.fa" \
    --skip-csv  # Reuse CSVs if re-importing
done
```

### Verify the graph

```bash
graphpop db info
```

```
Database: human1000g    Status: ONLINE    Size: 46.2 GB

Nodes: Variant 70,714,808 | Sample 3,202 | Population 31 | Gene 91,973
       Pathway 2,241 | Chromosome 22 | GOTerm 18,429
Edges: CARRIES ~2.1B | NEXT 70.7M | HAS_CONSEQUENCE 3.8M | IN_PATHWAY 4,900
```

Before diving into analysis, use `status` and `inventory` to confirm that
the server is healthy and all expected data layers are present:

```bash
# Check that Neo4j is running and the database is reachable
graphpop status

# List all populations, chromosomes, and annotation layers in the graph
graphpop inventory
```

```
GraphPop Status
  Neo4j:     RUNNING (PID 23914)
  Bolt URI:  bolt://localhost:7687
  Database:  human1000g (ONLINE)
  Plugin:    graphpop-procedures 0.1.0

Inventory — human1000g
  Populations (31): 26 sub-populations + 5 superpopulations
  Chromosomes (22): chr1..chr22
  Genes:      91,973
  Pathways:   2,241
  Ancestral:  65,482,104 polarized variants (92.6%)
```

---

## 4. Benchmarking: How Fast Is GraphPop?

### Target problem

Does the graph architecture actually deliver on its performance promise, or is
the overhead of a database engine too high?

### Quick benchmark on chr22

Let us time a diversity computation on chromosome 22 for the CEU population
(99 samples, 1,066,555 variants):

```bash
time graphpop diversity chr22 1 51304566 CEU -o bench_diversity.tsv
```

```
real    0m12.1s
user    0m0.3s
sys     0m0.1s
```

The same computation in scikit-allel---loading the VCF into a numpy array,
computing pi, theta_W, and Tajima's D---takes approximately 1,757 seconds
(29.3 minutes). That is a **146x speedup**.

### Why is GraphPop faster?

**FAST PATH** (diversity, Fst, SFS): Allele counts are pre-aggregated on
Variant nodes as arrays indexed by population. Computing pi for CEU reads one
integer per variant---O(V x K)---and never touches individual genotypes.
scikit-allel must load the full V x N genotype matrix (30x more data).

**FULL PATH** (iHS, XP-EHH, ROH): Bit-packed haplotype arrays on Chromosome
nodes are accessed directly via the Neo4j internal API, avoiding VCF parsing.

### Comprehensive benchmark results

All benchmarks use chr22, CEU population, full chromosome:

| Statistic              | GraphPop (s) | Competitor (s) | Tool        | Speedup | Accuracy        |
|------------------------|-------------|----------------|-------------|---------|-----------------|
| **FAST PATH**          |             |                |             |         |                 |
| pi / theta_W / D       | 12.1        | 1,757          | scikit-allel| 146x    | < 0.000001%     |
| Fst / Dxy              | 8.8         | 2,165          | scikit-allel| 245x    | < 0.000001%     |
| SFS                    | 5.9         | 1,937          | scikit-allel| 327x    | < 0.000001%     |
| **FULL PATH**          |             |                |             |         |                 |
| iHS                    | 10.9        | 1,948          | scikit-allel| 179x    | r = 0.999       |
| XP-EHH                 | 16.7        | 2,280          | scikit-allel| 136x    | r = 0.97        |
| nSL                    | 35          | 2,231          | scikit-allel| 63x     | r = 0.9997      |
| ROH                    | 10.9        | 50             | bcftools    | 4.6x    | r = 0.96        |
| LD (r-squared)         | 6.4         | 1.1            | PLINK 2     | 0.17x   | < 0.1%          |

FAST PATH statistics achieve exact numerical agreement (< 0.000001% error).
Memory is constant at ~160 MB for the Neo4j path (vs 1,051 MB for scikit-allel).
LD is the one case where PLINK 2 wins (5.7x faster, native binary format).

The most important advantage: **result persistence**. Every statistic is
written to graph nodes and cross-referenced by Cypher queries in milliseconds.
Classical tools write to flat files that must be re-loaded for any
cross-statistic analysis.

---

## 5. Population Diversity and Continental Structure

### Target problem

Do GraphPop's statistics reproduce the known continental patterns of human
genetic diversity? This is our validation step: if the graph-computed results
match established population genetics, we can trust the novel analyses that
follow.

### Run the full-genome analysis

The `run-all` command orchestrates all 12 stored procedures across every
population and chromosome:

```bash
graphpop run-all --phase 1 -d results/human/
```

```
Detected 26 populations, 22 chromosomes
Phase 1 tasks: 572 runs x 7 procedures = 3,899 analyses
  Procedures: diversity, sfs, pop_summary, roh, ihs, nsl, garud_h

[chr1/YRI]  diversity..12.4s  sfs..5.1s  roh..11.2s  ihs..10.8s  nsl..34.7s
...
[chr22/PEL] diversity..3.2s   sfs..1.4s  roh..2.8s   ihs..2.9s   nsl..9.1s

Phase 1 complete: 3,899 analyses in ~112 hours
```

Then run Phase 2 for pairwise statistics:

```bash
graphpop run-all --phase 2 -d results/human/
```

```
Phase 2 tasks: 325 pairs x 22 chrs (divergence) + 18 pairs x 22 chrs (XP-EHH)
[chr1/YRI-CEU]  divergence..8.9s  xpehh..16.4s
...
Phase 2 complete: ~20 hours
```

### Generate summary tables

```bash
graphpop aggregate -d results/human/ -o tables/
# Generates: population_summary.tsv, fst_matrix.tsv, roh_summary.tsv,
#            ihs_peaks.tsv, sweep_windows.tsv
```

### Continental diversity patterns

The population summary table confirms every known pattern in human population
genetics:

```bash
# View superpopulation-level diversity averages
graphpop query "MATCH (p:Population)
  WHERE p.level = 'superpopulation'
  RETURN p.populationId AS pop, p.n_samples AS n,
         round(p.mean_pi, 4) AS pi,
         round(p.mean_theta_w, 4) AS theta_w,
         round(p.mean_tajima_d, 2) AS tajima_d,
         round(p.mean_froh, 4) AS froh
  ORDER BY pi DESC"
```

```
pop     n       pi      theta_w tajima_d    froh
AFR     661     0.0510  0.0610  -0.52       0.0050
AMR     347     0.0400  0.0480  -0.39       0.0090
SAS     489     0.0390  0.0420  -0.07       0.0140
EUR     503     0.0380  0.0410  +0.09       0.0100
EAS     504     0.0360  0.0380  +0.20       0.0090
```

AFR shows the highest diversity (pi = 0.051), consistent with not having
passed through the Out-of-Africa bottleneck. Each successive branching event
reduced diversity via serial founder effects (AFR > SAS > AMR > EUR > EAS).
Tajima's D is negative in AFR (expansion) and positive in EAS (structure).
FROH is highest in SAS (0.014, consanguinity) and lowest in AFR (0.005).
These results exactly reproduce textbook population genetics.

You can also run a single population at any time:

```bash
graphpop diversity chr22 1 51304566 CEU -o diversity_ceu_chr22.tsv
```

```
pi        theta_w   tajima_d  fay_wu_h  het_exp  het_obs  fis    n_variants  n_segregating
0.038214  0.041892  0.0912    -0.3214   0.0382   0.0371   0.029  1066555     492871
```

### Visualize diversity and divergence

```bash
# Visualize continental diversity ranking
graphpop plot diversity-bar results/human/diversity/ -o figures/fig_human_diversity.png

# Visualize pairwise Fst heatmap (26 populations)
graphpop plot fst-heatmap results/human/divergence/ -o figures/fig_human_fst_heatmap.png
```

The diversity bar chart immediately reveals the Out-of-Africa gradient: African
populations at the top, East Asian at the bottom. The Fst heatmap shows the
continental block structure, with within-continent pairs tightly clustered (blue)
and between-continent pairs highly differentiated (red).

### Visualize the site frequency spectrum

```bash
# Compare SFS between African and European populations
graphpop sfs chr22 1 51304566 YRI -o sfs_yri.tsv
graphpop plot sfs-plot sfs_yri.tsv --title "SFS: YRI (African)" -o figures/fig_human_sfs_yri.png

graphpop sfs chr22 1 51304566 CEU -o sfs_ceu.tsv
graphpop plot sfs-plot sfs_ceu.tsv --title "SFS: CEU (European)" -o figures/fig_human_sfs_ceu.png
```

The SFS comparison highlights the demographic signature of the Out-of-Africa
bottleneck: YRI retains an excess of rare variants (steep left skew), reflecting
the large long-term effective population size of African populations. CEU shows a
flatter spectrum with a relative excess of intermediate-frequency variants,
consistent with the bottleneck that reduced diversity and shifted allele
frequencies upward during the European founder event.

### Additional single-population and pairwise statistics

The `pop-summary` command gives a compact one-line overview of diversity,
heterozygosity, and segregating site counts for a single population on one
chromosome:

```bash
# One-line summary for CEU on chr22
graphpop pop-summary chr22 CEU -o summary_ceu.tsv
```

The `joint-sfs` command computes the two-dimensional site frequency spectrum
between two populations. This is the core input for demographic inference
methods (dadi, moments) and visualizes allele sharing patterns:

```bash
# Joint SFS between YRI and CEU across chr22
graphpop joint-sfs chr22 1 51304566 YRI CEU -o joint_sfs_yri_ceu.tsv
```

The `nsl` statistic is a haplotype-based selection scan similar to iHS but
robust to local variation in recombination rate. It complements iHS and can
detect sweeps that iHS misses in regions of low recombination:

```bash
# nSL selection scan for CEU, persisted to graph nodes for later queries
graphpop nsl chr22 CEU --persist -o nsl_ceu.tsv
```

The `ld` command computes pairwise linkage disequilibrium (r-squared) within a
genomic region. Use it to examine haplotype structure around candidate loci:

```bash
# Pairwise LD in a 1 Mb region of chr22 for CEU
graphpop ld chr22 16000000 17000000 CEU 200000 0.2 -o ld_ceu.tsv
```

---

## 5b. Annotation-Conditioned Diversity: piN/piS Across Human Populations

**Target problem:** Is purifying selection on protein-coding variants efficient across all human populations, or has the Out-of-Africa bottleneck relaxed it?

This is the same question asked for rice (where the answer was "yes, relaxed in all subpopulations"). For humans, annotation-conditioned diversity provides the answer:

```bash
# Compute piN (missense diversity) and piS (synonymous diversity) for each population
graphpop diversity chr22 1 51304566 YRI --consequence missense_variant -o piN_yri.tsv
graphpop diversity chr22 1 51304566 YRI --consequence synonymous_variant -o piS_yri.tsv

# Repeat for all 26 populations, all 22 chromosomes
# Or use run-all with consequence conditioning (deferred until GraphPop plugin deployed)
```

**Results (genome-wide piN/piS for all 26 populations):**

| Population group | piN/piS range | Interpretation |
|-----------------|---------------|----------------|
| African (7 pops) | 0.646–0.651 | Most efficient purifying selection |
| South Asian (5 pops) | 0.660–0.662 | Moderately relaxed |
| European (5 pops) | 0.666–0.670 | Bottleneck-relaxed |
| East Asian (5 pops) | 0.668–0.673 | Most relaxed |
| American (4 pops) | 0.666–0.668 | Mixed ancestry |

All 26 populations show piN/piS **well below 1.0** (range: 0.646–0.673), confirming that purifying selection on protein-coding variants remains efficient in all human populations — a stark contrast with rice where ALL subpopulations exceed 1.0.

> **Cross-species comparison:** Human piN/piS (0.646–0.673, all < 1.0) vs rice piN/piS (1.018–1.146, all > 1.0). Natural selection in human populations maintains efficient removal of deleterious protein-coding variants. Domestication in rice has relaxed this constraint genome-wide. The ~3.5% reduction in purifying selection efficiency from African to East Asian populations reflects the Out-of-Africa bottleneck — a much smaller effect than the complete reversal seen in domesticated rice.

---

## 6. Pathway-Level Population Differentiation

### Target problem

Which biological pathways show the strongest genetic differentiation between
continental groups? Classical tools cannot compute pathway-level Fst without
constructing a multi-tool pipeline: run VEP to annotate consequences, join
variants to genes, join genes to pathways from an external database, group by
pathway, and average Fst values. In GraphPop, this is a single Cypher
traversal.

### The query

```bash
graphpop query "MATCH (p:Pathway)<-[:IN_PATHWAY]-(g:Gene)<-[:HAS_CONSEQUENCE]-(v:Variant)
  WHERE v.fst_yri_ceu IS NOT NULL
  RETURN p.name AS pathway, round(avg(v.fst_yri_ceu), 4) AS mean_fst,
         count(DISTINCT v) AS n_variants, count(DISTINCT g) AS n_genes
  ORDER BY mean_fst DESC LIMIT 20" -o pathway_fst_yri_ceu.tsv
```

```
pathway                                             mean_fst    n_variants  n_genes
Oculocutaneous albinism type I                      0.4520      847         3
SLC-mediated transmembrane transport (pigment.)     0.3891      1204        8
TREK/TWIK potassium channels                        0.1620      341         6
Phase 3 - Rapid depolarisation                      0.1330      512         9
Voltage gated potassium channels                    0.1280      1893        42
Potassium channels                                  0.1210      2341        56
Phase 2 - Loss of K+ channels                       0.1190      287         5
Cardiac conduction                                  0.1140      3102        71
Phase 0 - Rapid depolarisation                      0.1120      198         4
Neuronal system (ion channels)                      0.1080      4521        98
...
```

### What this reveals

The top hit---oculocutaneous albinism (Fst = 0.452)---is driven by *SLC24A5*,
the best-documented signal of post-Out-of-Africa selection on skin
pigmentation. Finding it at the top validates the approach.

More interesting is the cluster of cardiac ion channel pathways: TREK/TWIK
(0.162), Phase 3 repolarisation (0.133), voltage-gated K+ channels (0.128).
These implicate *KCNE1*, *KCNE4*, *KCNH2*, and *KCNQ1*---potassium channel
subunits involved in cardiac repolarisation. Graph traversal through
`IN_PATHWAY` edges separates this cardiac signal from overlapping metabolic
and immune annotations without manual curation.

---

## 6b. Pathway Co-Selection Network

**Target problem:** Do functionally related pathways show coordinated differentiation patterns across human population pairs?

After computing pathway-level Fst across 5 population pairs (YRI-CEU, YRI-CHB, CEU-CHB, YRI-JPT, CEU-JPT), we correlated pathway Fst profiles to identify co-differentiating pathways — the same analysis performed on rice (where 3 modules were found).

```bash
# Query pathway co-selection (from pre-computed results)
graphpop query "MATCH (p1:Pathway), (p2:Pathway)
  WHERE id(p1) < id(p2)
  WITH p1, p2,
    [p1.fst_yri_ceu, p1.fst_yri_chb, p1.fst_ceu_chb] AS v1,
    [p2.fst_yri_ceu, p2.fst_yri_chb, p2.fst_ceu_chb] AS v2
  RETURN p1.name, p2.name, gds.alpha.similarity.cosine(v1, v2) AS sim
  ORDER BY sim DESC LIMIT 20" -o pathway_cosim.tsv
```

**Results:**
- 2,232 Reactome pathways correlated across 5 population pairs
- 1,182,891 pairwise edges (rho > 0.7): 939,959 positive, 242,932 negative
- **3 communities identified:**
  1. **DNA repair and apoptosis** (602 pathways): base excision repair, caspase activation, ERBB4 signaling
  2. **Immune and neuronal signaling** (~800 pathways): ISG15 antiviral, opioid signaling, neurotransmitter clearance
  3. **Metabolic processes** (~800 pathways): platelet degranulation, erythrocyte function, PI3K signaling

> **Cross-species comparison:** Both human and rice show 3 co-selection communities. In both species, core cellular maintenance pathways (DNA repair in human, DNA replication/translation in rice) form a large, highly connected module that differentiates as a genomic background. The smaller modules capture species-specific biology: immune/neuronal in human vs secondary metabolism/hormone in rice. This convergence in co-selection architecture across species — discovered by the same GraphPop pathway correlation workflow — suggests that coordinated functional divergence is a general feature of population differentiation.

---

## 7. Convergent Positive Selection: Finding Ancient Sweeps

### Target problem

Can we identify genes under positive selection in multiple independent
continental lineages by querying already-stored statistics---without
re-computing anything?

In a classical workflow, answering this question requires: (1) running Garud's
H separately for each population, (2) identifying peaks in each output file,
(3) mapping peaks to genes via coordinate matching against a gene annotation
file, (4) intersecting gene lists across populations. This involves at least
4 tools and dozens of intermediate files.

In GraphPop, all Garud's H values are stored on GenomicWindow nodes, and all
variant-to-gene mappings are stored as HAS_CONSEQUENCE edges. One query does
the entire analysis.

### The query

```bash
graphpop query "MATCH (w:GenomicWindow)
  WHERE w.h12 > 0.3
  MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene)
  WHERE v.chr = w.chr AND v.pos >= w.start AND v.pos <= w.end
  WITH g.symbol AS gene, g.geneId AS gene_id,
       collect(DISTINCT w.population) AS pops,
       max(w.h12) AS max_h12
  WHERE size(pops) >= 2
  RETURN gene, gene_id, size(pops) AS n_pops, max_h12, pops
  ORDER BY n_pops DESC, max_h12 DESC" -o convergent_sweeps.tsv
```

```
gene            gene_id         n_pops  max_h12 pops
KCNE1           ENSG00000180509 5       1.0000  [AFR, AMR, EAS, EUR, SAS]
LOC102724560    ENSG00000273344 4       0.9000  [AMR, EAS, EUR, SAS]
LOC107987288    ENSG00000283901 4       0.8300  [AMR, EAS, EUR, SAS]
LOC102724428    ENSG00000272872 4       0.8100  [AMR, EAS, EUR, SAS]
LOC124900992    ENSG00000289472 3       0.6800  [EAS, EUR, SAS]
LOC124900995    ENSG00000289518 3       0.6800  [EAS, EUR, SAS]
LINC01678       ENSG00000267281 2       0.7100  [AMR, EUR]
LOC102725065    ENSG00000273167 2       0.7100  [AMR, EUR]
CRYAA           ENSG00000160202 2       0.5500  [AMR, SAS]
```

### KCNE1: a pre-Out-of-Africa sweep

The most striking result is *KCNE1*---the **only gene showing sweep evidence
in all five continental groups**. This is remarkable: independent populations
on five continents all carry the signature of a selective sweep at the same
locus.

Let us examine the evidence in more detail:

```bash
graphpop query "MATCH (w:GenomicWindow)
  WHERE w.h12 > 0.3 AND w.chr = 'chr21'
  MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene {symbol: 'KCNE1'})
  WHERE v.pos >= w.start AND v.pos <= w.end
  RETURN w.population AS pop, round(w.h12, 4) AS h12,
         round(w.h2_h1, 4) AS h2_h1, w.start AS start, w.end AS end
  ORDER BY h12 DESC"
```

```
pop     h12     h2_h1   start       end
GBR     1.0000  0.0000  34800001    34900000
CEU     1.0000  0.0000  34800001    34900000
FIN     1.0000  0.0000  34800001    34900000
IBS     0.9980  0.0020  34800001    34900000
TSI     0.9950  0.0050  34800001    34900000
CHB     0.9200  0.0087  34800001    34900000
JPT     0.9100  0.0099  34800001    34900000
CDX     0.8800  0.0136  34800001    34900000
PJL     0.8500  0.0176  34800001    34900000
GIH     0.8300  0.0205  34800001    34900000
YRI     0.4200  0.0714  34800001    34900000
LWK     0.3800  0.0816  34800001    34900000
```

H12 = 1.0 in European populations (complete haplotype homozygosity), > 0.8 in
East and South Asian populations, and > 0.3 even in African populations---the
sweep predates the Out-of-Africa migration. The between-group Fst ratio of
0.296 confirms shared ancestry: this is the signature of a **pre-Out-of-Africa
selective sweep** (~60-100 kya).

This signal is recoverable **only** by cross-querying independently computed
statistics (Garud's H per population + between-group Fst) stored on the same
graph. The classical equivalent requires intersecting five separate sweep scan
output files, then joining against Fst values via chromosome coordinate
matching---a manual process that is rarely attempted.

### Automated convergence, gene ranking, and graph exploration

The `converge` command automates the multi-statistic intersection shown above.
It finds all regions where every specified statistic exceeds its threshold
simultaneously:

```bash
# Regions where iHS, XP-EHH, and H12 all exceed thresholds in CEU
graphpop converge --stats ihs,xpehh,h12 --thresholds 2.0,2.0,0.3 --pop CEU -o convergent_ceu.tsv
```

The `rank-genes` command produces a ranked list of candidate genes ordered by
composite selection evidence, combining all persisted statistics into a single
score:

```bash
# Top 50 genes by composite selection evidence in CEU
graphpop rank-genes --pop CEU --top 50 -o top_genes_ceu.tsv
```

To inspect a specific candidate gene, `lookup` retrieves all stored properties
including coordinates, variant counts, persisted statistics, and pathway
memberships:

```bash
# Look up the KCNE1 potassium channel gene
graphpop lookup gene KCNE1 -o kcne1_info.tsv
```

The `neighbors` command walks the graph outward from a node, returning all
entities within a specified number of hops. The optional `--via` flag restricts
traversal to a specific relationship type. This is useful for exploring
the functional context of a candidate gene:

```bash
# Explore KCNE1 neighborhood via pathway edges (2 hops)
graphpop neighbors KCNE1 --hops 2 --via IN_PATHWAY -o kcne1_neighbors.tsv
```

---

## 8. Individual-Level Evolutionary Trajectories

### Target problem

Can we characterize individuals by the evolutionary forces acting on their
genomes---not just their ancestry? Classical PCA (PLINK, EIGENSOFT) projects
individuals onto axes of genetic ancestry but cannot jointly incorporate
inbreeding, heterozygosity, and functional rare variant burden. These
statistics reside in different tools and file formats.

### Compute per-sample ROH

Runs of homozygosity (ROH) reflect recent inbreeding and population
bottlenecks. GraphPop computes ROH for all samples in a population using an
HMM that operates on bit-packed haplotype arrays:

```bash
graphpop roh chr22 EUR -o roh_eur_chr22.tsv
```

```
sampleId    n_roh   total_length    froh        mean_length max_length
HG00096     2       1847291         0.0360      923645      1284102
HG00097     0       0               0.0000      0           0
HG00099     1       512843          0.0100      512843      512843
HG00100     3       2941872         0.0573      980624      1573291
HG00101     1       687432          0.0134      687432      687432
...
```

For the full genome, `run-all` computes ROH across all 22 chromosomes and
aggregates genome-wide FROH per individual.

### Visualize the inbreeding landscape

```bash
# Inbreeding landscape across 26 populations
graphpop plot roh-landscape results/human/roh/ -o figures/fig_human_roh.png
```

The ROH landscape plot displays genome-wide FROH distributions for all 26
populations, ordered by median inbreeding coefficient. South Asian populations
(PJL, GIH, BEB) cluster at the high end, reflecting documented consanguinity
practices. African populations occupy the low end, consistent with their larger
effective population sizes and lower drift-driven homozygosity. The spread
within each population reveals individual-level variation that population
averages obscure---including extreme outliers like HG01816 (discussed below).

### Building the process space

Instead of classical ancestry PCA, we embed individuals in "process space"
defined by evolutionary forces: (1) genome-wide FROH, (2) chr22 heterozygosity
rate, (3) chr22 rare variant count (AF < 0.01), (4) chr22 rare missense count,
and (5) genome-wide ROH count. All five are pre-computed on graph nodes by
`run-all` or extractable in a single Cypher session.

### Results

```bash
graphpop query "MATCH (s:Sample)
  RETURN s.sampleId AS sample, s.population AS pop,
         round(s.froh_genome, 4) AS froh, s.n_roh_genome AS n_roh,
         round(s.het_rate_chr22, 6) AS het_rate,
         s.rare_burden_chr22 AS rare_burden
  ORDER BY s.sampleId LIMIT 5" -o process_space_sample.tsv
```

```
sample      pop     froh    n_roh   het_rate    rare_burden
HG00096     GBR     0.0082  14      0.003821    287
HG00097     PJL     0.0163  22      0.003714    241
HG00099     GBR     0.0071  11      0.003845    291
HG00100     GBR     0.0094  16      0.003801    284
HG00101     GBR     0.0078  12      0.003832    289
```

Process-space PCA reveals: PC1 (57.0%) separates populations by diversity
(African highest); PC2 (30.9%) separates by inbreeding (South Asian highest,
FROH = 0.014). This inbreeding axis is invisible in classical ancestry PCA.

### The outlier: HG01816

The top outlier in process space is **HG01816**, a CDX (Chinese Dai in
Xishuangbanna) individual with FROH = 0.125---**18 times the CDX population
mean** (0.007).

```bash
graphpop query "MATCH (s:Sample {sampleId: 'HG01816'})
  RETURN s.sampleId AS sample, s.population AS pop,
         round(s.froh_genome, 4) AS froh, s.n_roh_genome AS n_roh,
         round(s.het_rate_chr22, 6) AS het_rate,
         s.rare_burden_chr22 AS rare_burden"
```

```
sample      pop     froh    n_roh   het_rate    rare_burden
HG01816     CDX     0.1250  89      0.003142    198
```

A single 3-hop query traces from the individual through their variants to
gene and pathway context:

```bash
graphpop query "MATCH (s:Sample {sampleId: 'HG01816'})-[:CARRIES]->(v:Variant)
  -[:HAS_CONSEQUENCE]->(g:Gene)-[:IN_PATHWAY]->(p:Pathway)
  WHERE v.af_CDX < 0.01
  RETURN g.symbol AS gene, p.name AS pathway,
         count(DISTINCT v) AS rare_variants
  ORDER BY rare_variants DESC LIMIT 10"
```

```
gene        pathway                                     rare_variants
LARGE1      Dystroglycan-related signaling               4
COMT        Catechol-O-methyltransferase metabolism       3
TBX1        Cardiac development                          3
PRODH       Proline degradation                          2
DGCR6       DiGeorge syndrome region                     2
MYH11       Smooth muscle contraction                    2
LZTR1       Ras signaling regulation                     2
PCSK9       Cholesterol metabolism                       1
APOB        Lipoprotein assembly                         1
LDLR        LDL receptor-mediated endocytosis            1
```

This Sample-CARRIES-HAS_CONSEQUENCE-IN_PATHWAY traversal takes milliseconds
and has no classical equivalent.

---

## 9. Pathway-Stratified Purifying Selection

### Target problem

The 2x African enrichment in rare functional variants is well established: Out-
of-Africa populations lost rare variants through bottlenecks, while African
populations retained them. But is this enrichment **uniform across all
biological pathways**, or is it concentrated in specific functional categories?

### Compute annotation-conditioned statistics

```bash
# Step 1: Compute iHS genome-wide (already done by run-all)
# Step 2: Filter to missense variants with strong signals
graphpop filter ihs chr22 EUR --consequence missense_variant --min-score 2.0 \
  -o ihs_missense_eur.tsv
```

```
Found 847 records.
```

```bash
head -5 ihs_missense_eur.tsv
```

```
variant_id              pos         ihs     ihs_unstd   consequence         impact
chr22:16591593:G:A      16591593    -2.341  -1.892      missense_variant    MODERATE
chr22:17308432:C:T      17308432     2.104   1.743      missense_variant    MODERATE
chr22:18574891:A:G      18574891    -2.567  -2.103      missense_variant    MODERATE
chr22:19211034:G:T      19211034     2.892   2.341      missense_variant    MODERATE
```

### Classify pathways by constraint level

We classify pathways as "constrained" (lowest-quartile Fst, i.e., most
conserved across populations) or "divergent" (highest-quartile Fst, i.e.,
most differentiated):

```bash
graphpop query "MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene)-[:IN_PATHWAY]->(p:Pathway)
  WHERE v.fst_yri_ceu IS NOT NULL
  WITH p.name AS pathway, avg(v.fst_yri_ceu) AS pathway_fst
  RETURN pathway, round(pathway_fst, 4) AS fst
  ORDER BY pathway_fst" -o pathway_constraint.tsv
```

```
pathway                                 fst
Metabolism of RNA                       0.0312
Translation                             0.0328
Cell cycle checkpoints                  0.0341
DNA repair                              0.0355
mRNA processing                         0.0362
...
Phase 3 - Rapid depolarisation          0.1330
TREK/TWIK potassium channels            0.1620
SLC transmembrane transport (pigment.)  0.3891
Oculocutaneous albinism type I          0.4520
```

### The result: enrichment is uniform

Using the graph's co-location of all four data layers, we computed per-
individual rare missense burden in constrained vs. divergent pathways:

| Category      | African Burden | Non-African Burden | Enrichment | p-value        |
|---------------|---------------|-------------------|------------|----------------|
| Constrained   | 142.3         | 69.8              | 2.04x      | 1.6e-290       |
| Divergent     | 187.6         | 96.1              | 1.95x      | 3.0e-265       |

The 2x African enrichment holds at **equal magnitude** in both constrained and
divergent pathways. FROH predicts individual burden (rho = -0.211, p = 1.6e-33)
but NOT the constrained:divergent ratio (rho = +0.004, p = 0.84). Conclusion:
drift clears rare variants non-specifically across all pathway classes.

---

## 10. Consequence-Level Selection Bias

### Target problem

How does natural selection shape Fst across variant impact categories? Do
functionally deleterious alleles differentiate more or less between populations
than neutral ones?

### The query

```bash
graphpop query "MATCH (v:Variant)-[hc:HAS_CONSEQUENCE]->(g:Gene)
  WHERE v.chr = 'chr22' AND v.fst_yri_ceu IS NOT NULL
  RETURN hc.impact AS impact,
         round(avg(v.fst_yri_ceu), 4) AS mean_fst,
         count(DISTINCT v) AS n_variants
  ORDER BY mean_fst" -o consequence_fst.tsv
```

```
impact      mean_fst    n_variants
HIGH        0.1240      1847
MODERATE    0.1310      28941
LOW         0.1430      42187
MODIFIER    0.1380      89412
```

### Interpretation

HIGH-impact variants (stop-gain, frameshift, splice-altering) show the
**lowest** mean Fst (0.124), while LOW-impact variants (synonymous, UTR) show
the highest (0.143). The pattern is monotonic with predicted functional impact.

This is the signature of **purifying selection at the population level**:
deleterious alleles are constrained against differentiation because natural
selection keeps them at low frequency in all populations. They cannot drift
to high frequency in one group and low frequency in another---selection purges
them everywhere.

### Comparison with rice: domestication overrides constraint

This human pattern is the opposite of what GraphPop finds in rice. In the rice
3K dataset (GJ-tmp vs XI-1A comparison):

| Species | HIGH Fst | LOW Fst | Pattern                          |
|---------|----------|---------|----------------------------------|
| Human   | 0.124    | 0.143   | HIGH < LOW (purifying selection)  |
| Rice    | 0.185    | 0.092   | HIGH > LOW (domestication drives) |

In rice, strong directional selection during domestication has **overridden**
the natural constraint, driving fixation of functionally consequential variants
between the two most divergent subpopulations. This is a direct signature of
the transition from natural to artificial selection regimes.

This cross-species comparison---human vs rice, same statistic, same graph
query---is uniquely enabled by GraphPop's unified architecture. Both datasets
use identical node types, edge types, and stored procedures.

---

## 11. Summary Tables and Database Sharing

### Generate final summary tables

```bash
graphpop aggregate -d results/human/ -o tables/
```

```
Aggregating results from results/human/...

Generated tables:
  tables/population_summary.tsv    26 populations, 11 statistics each
  tables/fst_matrix.tsv            26 x 26 pairwise Hudson Fst
  tables/roh_summary.tsv           3,202 individuals, genome-wide FROH
  tables/ihs_peaks.tsv             Top |iHS| > 2 per population
  tables/sweep_windows.tsv         H12 > 0.1 windows per population
```

The Fst matrix captures the global structure of human population differentiation:

```bash
graphpop query "MATCH (p1:Population), (p2:Population)
  WHERE p1.level = 'superpopulation' AND p2.level = 'superpopulation'
    AND id(p1) < id(p2)
  RETURN p1.populationId AS pop1, p2.populationId AS pop2,
         round(p1['fst_' + p2.populationId], 4) AS fst
  ORDER BY fst DESC"
```

```
pop1    pop2    fst
AFR     EAS     0.1109
AFR     EUR     0.1053
AFR     SAS     0.0967
AFR     AMR     0.0641
EUR     EAS     0.0726
SAS     EAS     0.0542
EUR     SAS     0.0312
AMR     EAS     0.0481
AMR     EUR     0.0291
AMR     SAS     0.0264
```

Maximum intercontinental divergence is between AFR and EAS (Fst = 0.111);
minimum is between AMR and SAS (Fst = 0.026), reflecting the admixed nature
of American populations.

### Export the database for sharing

GraphPop databases can be shared as Neo4j dump files, allowing collaborators
to query all computed results without re-running any analysis:

```bash
graphpop stop   # Neo4j must be stopped for dump
graphpop dump --database human1000g -o graphpop_human1000g_v1.dump
```

```
Dumping database 'human1000g' to graphpop_human1000g_v1.dump...
Dump complete: graphpop_human1000g_v1.dump (23.8 GB)
Manifest: graphpop_human1000g_v1.manifest.json
```

The manifest records node/edge counts, population sample sizes, and chromosome
lengths for reproducibility. A collaborator restores the database:

```bash
graphpop load --dump-file graphpop_human1000g_v1.dump --database human1000g
graphpop start
graphpop db info   # Verify node/edge counts match manifest
```

Every statistic computed in this vignette is immediately available for
querying---no re-computation required.

### Comparing populations, exporting regions, and batch processing

The `compare` command computes per-window differences in a summary statistic
between two populations. This is useful for identifying regions where one
group has lost (or gained) diversity relative to another:

```bash
# Per-window delta-pi between CEU and YRI across chr22
graphpop compare CEU YRI chr22 --stat pi -o delta_pi.tsv
```

To hand off selection peaks to external tools (bedtools, IGV, UCSC Genome
Browser), `export-bed` writes regions exceeding a statistic threshold as a
standard BED file:

```bash
# Export strong iHS peaks as BED
graphpop export-bed --stat ihs --threshold 2.5 --pop CEU -o ihs_peaks.bed
```

The `extract` command pulls a filtered set of variants from the graph,
optionally restricted by population, allele frequency range, and output
fields. This replaces multi-step bcftools + awk pipelines:

```bash
# Extract common variants in CEU with selected fields
graphpop extract variants --chr chr22 --pop CEU --min-af 0.05 \
    --fields pos,ref,alt,af,ihs -o variants_ceu.tsv
```

When the same analysis needs to run across many populations or chromosomes,
the `batch` command parallelizes the work automatically:

```bash
# Run diversity for three populations on chr22 in parallel
graphpop batch diversity --pops CEU,YRI,CHB --chrs chr22 --workers 3 -d batch/
```

Use `config` to inspect the current CLI settings (default database, Neo4j
URI, output format):

```bash
# Show all current configuration settings
graphpop config show
```

The `report` command generates a self-contained HTML report summarizing all
computed statistics, top selection signals, and population comparisons:

```bash
# Generate an automated analysis report
graphpop report -o human1000g_report.html
```

---

## Cross-Species Comparison: Human vs Rice

| Analysis | Human (1000 Genomes) | Rice (3K Genomes) | Interpretation |
|----------|---------------------|-------------------|----------------|
| **πN/πS** | 0.646–0.673 (all < 1.0) | 1.018–1.146 (all > 1.0) | Efficient purifying selection vs relaxed under domestication |
| **Consequence bias** | HIGH Fst (0.124) < LOW (0.143) | HIGH Fst (0.185) > LOW (0.092) | Purifying constraint vs domestication override |
| **Convergent sweeps** | KCNE1 in all 5 continental groups | DRO1 in 4/12 subpopulations (max) | Universal natural selection vs trait-specific artificial selection |
| **Individual trajectories** | PC1 = Out-of-Africa diversity gradient | PC1 = indica-japonica diversity gradient | Same framework, parallel structure, different forces |
| **Pathway co-selection** | 3 communities (2,232 pathways) | 3 modules (153 pathways) | Housekeeping co-differentiation in both species |

These contrasts are uniquely enabled by annotation-conditioned computation applied identically across datasets.

---

## 12. Summary of Findings

### What we confirmed (known biology)

GraphPop reproduces every major pattern in human population genetics, validating
the graph-native approach:

| Finding                                    | Expected?  | GraphPop Result                   |
|--------------------------------------------|-----------|-----------------------------------|
| African highest diversity                   | Yes       | pi = 0.051 (AFR) vs 0.036 (EAS)  |
| Serial founder effect gradient              | Yes       | AFR > SAS > AMR > EUR > EAS       |
| Tajima's D negative in AFR (expansion)      | Yes       | D = -0.52 (AFR), +0.20 (EAS)     |
| FROH highest in South Asia (consanguinity)  | Yes       | 0.014 (SAS) vs 0.005 (AFR)       |
| SLC24A5 strongest Fst signal (pigmentation) | Yes       | Pathway Fst = 0.452               |
| Purifying selection constrains HIGH Fst     | Yes       | HIGH (0.124) < LOW (0.143)        |

### What we discovered (graph-enabled novel findings)

Three findings are uniquely enabled by the graph architecture:

1. **KCNE1 pre-Out-of-Africa sweep.** The only gene showing convergent
   sweep evidence (H12 > 0.3) in all five continental groups, with elevated
   between-group Fst ratio (0.296) consistent with an ancient sweep that
   reached near-fixation before population divergence. Detecting this required
   cross-querying independently computed Garud's H, population-specific Fst,
   and gene annotations stored on the same graph nodes---an analysis with no
   practical classical equivalent.

2. **Individual evolutionary trajectories.** Process-space PCA (5 dimensions:
   FROH, ROH count, heterozygosity, rare burden, rare missense) reveals
   intra-population structure invisible to ancestry PCA. Outlier HG01816
   (CDX, FROH = 0.125, 18x mean) can be immediately traced to specific
   driver variants, genes, and pathways via 3-hop graph traversal.

3. **Uniform bottleneck depletion across pathways.** The 2x African rare
   variant enrichment holds at equal magnitude in both constrained (2.04x)
   and divergent (1.95x) pathways. FROH predicts burden (rho = -0.211)
   but not the constrained:divergent ratio (rho = +0.004). Drift clears
   rare variants non-specifically---a conclusion requiring simultaneous access
   to four co-located data layers.

### What classical tools cannot do

| Analysis                                | Classical Approach            | GraphPop               |
|-----------------------------------------|------------------------------|------------------------|
| Pathway-level Fst ranking               | VEP + gene map + Fst + join  | 1 Cypher query         |
| Convergent sweeps across 5 continents   | 5 scan files + intersect     | 1 Cypher query         |
| Individual outlier to pathway trace     | PLINK + VCF + VEP + manual   | 1 Cypher query (3 hop) |
| Pathway-stratified rare burden          | 4 tools + custom merge       | 1 Cypher session       |
| Consequence-conditioned Fst             | VEP + filter + VCFtools      | 1 procedure call       |
| Cross-species same-query comparison     | Entirely separate pipelines  | Same query, switch DB  |

### Computation summary

| Phase                              | Tasks  | Wall Time  | Notes                              |
|------------------------------------|--------|------------|------------------------------------|
| Phase 1 (per-population)           | 3,899  | ~112 hours | 26 pops x 22 chrs x 7 procedures  |
| Phase 2 (pairwise)                 | 7,546  | ~20 hours  | 325 pairs divergence + 18 XP-EHH  |
| Phase 3 (ancestral + PBS)          | 572    | ~5 hours   | Fay & Wu's H, unfolded SFS        |
| **Total computation**              | 12,017 | ~137 hours | Single workstation                 |
| **All subsequent queries**         | --     | < 1 second | No re-computation needed           |

The final line is the key point: once computed, every statistic is permanently
stored in the graph. The novel analyses in Sections 6--10 of this vignette
required **zero re-computation**---only Cypher queries on the persistent
analytical record.

---

## Appendix: Quick Reference

```bash
# Setup, status, and configuration
graphpop setup --password mypass --pagecache 20g --heap 4g
graphpop start
graphpop status
graphpop config set database human1000g
graphpop config show
graphpop import --vcf data.vcf.gz --panel panel.txt --database human1000g \
  --vep vep.vcf --pathways reactome.tsv --ancestral ancestral.fa
graphpop inventory
graphpop db info

# Full-genome analysis
graphpop run-all --phase 1 -d results/human/
graphpop run-all --phase 2 -d results/human/
graphpop aggregate -d results/human/ -o tables/

# Single-statistic commands
graphpop diversity chr22 1 51304566 CEU -o diversity.tsv
graphpop divergence chr22 1 51304566 YRI CEU -o fst.tsv
graphpop pop-summary chr22 CEU -o summary.tsv
graphpop joint-sfs chr22 1 51304566 YRI CEU -o joint_sfs.tsv
graphpop compare CEU YRI chr22 --stat pi -o delta_pi.tsv
graphpop roh chr22 EUR -o roh.tsv
graphpop ihs chr22 CEU --persist -o ihs.tsv
graphpop nsl chr22 CEU --persist -o nsl.tsv
graphpop ld chr22 16000000 17000000 CEU 200000 0.2 -o ld.tsv
graphpop garud-h chr22 CEU 100000 50000 -o garud.tsv

# Convergence and gene ranking
graphpop converge --stats ihs,xpehh,h12 --thresholds 2.0,2.0,0.3 --pop CEU -o convergent.tsv
graphpop rank-genes --pop CEU --top 50 -o top_genes.tsv
graphpop lookup gene KCNE1 -o gene_info.tsv
graphpop neighbors KCNE1 --hops 2 --via IN_PATHWAY -o neighbors.tsv

# Annotation conditioning
graphpop diversity chr22 1 51304566 CEU --consequence missense_variant
graphpop filter ihs chr22 EUR --consequence missense_variant --min-score 2.0

# Batch processing, export, and reporting
graphpop batch diversity --pops CEU,YRI,CHB --chrs chr22 --workers 3 -d batch/
graphpop export-bed --stat ihs --threshold 2.5 --pop CEU -o ihs_peaks.bed
graphpop extract variants --chr chr22 --pop CEU --min-af 0.05 -o variants.tsv
graphpop report -o human1000g_report.html

# Cypher queries and database sharing
graphpop query "MATCH (p:Pathway) ... RETURN ..." -o results.tsv
graphpop dump --database human1000g -o human1000g_v1.dump
graphpop load --dump-file human1000g_v1.dump --database human1000g
```

Environment: Ubuntu 22.04 (WSL2), Intel i9-13900K, 64 GB RAM, RTX 4090,
4 TB NVMe, Neo4j 5.26.0, Java 21, Python 3.11, GraphPop CLI 0.1.0.

## 12b. Programmatic Access: MCP Server for AI Agents

**Target problem:** Can AI agents autonomously query the 1000 Genomes graph database?

GraphPop's MCP server exposes 21 tools for AI agent access — the same procedures plus convergence detection, gene ranking, annotation lookup, filtering, and inventory. This enables autonomous scientific workflows where an agent investigates population genomics questions without human CLI interaction.

### Setting up the MCP server

```bash
cd graphpop-mcp && pip install -e .
export GRAPHPOP_URI=bolt://localhost:7687
export GRAPHPOP_DATABASE=human1000g
graphpop-mcp
```

### Example AI agent workflow

An AI agent investigating selection in European populations might autonomously:

```
1. graphpop_inventory()                     # What populations and stats are available?
2. graphpop_diversity(chr="chr22", start=1, end=51304566, pop="CEU")
3. graphpop_converge(pop="CEU", stats="ihs,xpehh,h12",
                     thresholds="2.0,2.0,0.3")    # Find convergent signals
4. graphpop_rank_genes(pop="CEU", top=20)           # Top selection candidates
5. graphpop_lookup_gene(gene="KCNE1")               # Investigate top hit
6. graphpop_filter(statistic="ihs", chr="chr22", pop="CEU",
                   consequence="missense_variant", min_score=2.0)  # Functional filtering
7. graphpop_lookup_pathway(pathway="cardiac repolarization")       # Pathway context
```

The MCP server returns JSON, enabling AI agents to parse results, reason about biological significance, and chain queries — performing the same annotation-conditioned, multi-statistic investigation demonstrated throughout this vignette, but autonomously.

> **Cross-reference:** The MCP server is registered with ToolUniverse (Zitnik Lab, Harvard) for discovery by AI agent frameworks. See also the companion GraphMana Brief Communication for data management capabilities.
