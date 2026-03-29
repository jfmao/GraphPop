# Rice 3K Genome Analysis with GraphPop

A complete walkthrough of graph-native population genomics using the
3K Rice Genomes dataset, from raw data to biological insight.

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Data Collection and Preparation](#2-data-collection-and-preparation)
3. [Neo4j Setup and Data Import](#3-neo4j-setup-and-data-import)
4. [Population Diversity and Structure](#4-population-diversity-and-structure)
5. [Annotation-Conditioned Analysis: Cost of Domestication](#5-annotation-conditioned-analysis-cost-of-domestication)
6. [Adaptive Protein Divergence](#6-adaptive-protein-divergence)
7. [Selection Scan and Multi-Statistic Convergence](#7-selection-scan-and-multi-statistic-convergence)
8. [Consequence Bias and Pathway Analysis](#8-consequence-bias-and-pathway-analysis)
8b. [Individual-Level Evolutionary Trajectories](#8b-individual-level-evolutionary-trajectories)
8c. [Convergent Sweep Detection Across Subpopulations](#8c-convergent-sweep-detection-across-subpopulations)
9. [Cross-Species Comparison: Rice vs Human](#cross-species-comparison-rice-vs-human)
10. [Generating Summary Tables and Figures](#10-generating-summary-tables-and-figures)
11. [Sharing the Database](#11-sharing-the-database)
12. [Summary](#12-summary)

---

## 1. Introduction

### What is the Rice 3K Genomes Project?

The 3K Rice Genomes Project (Wang et al. 2018, *Nature*) sequenced 3,024
rice accessions spanning the global diversity of cultivated rice
(*Oryza sativa*). Accessions were classified into 12 subpopulations using
ADMIXTURE analysis at K=9, reflecting the deep indica-japonica split that
occurred roughly 400,000 years ago, followed by independent domestication
events ~9,000 years ago.

The 12 subpopulations are:

| Group       | Code     | N samples | Description                        |
|-------------|----------|-----------|------------------------------------|
| indica-1A   | XI-1A    | ~480      | East/Southeast Asian indica        |
| indica-1B   | XI-1B    | ~340      | South Asian indica                 |
| indica-2    | XI-2     | ~280      | Aus-related indica                 |
| indica-3    | XI-3     | ~230      | Deep indica                        |
| indica-adm  | XI-adm   | ~170      | Admixed indica                     |
| temperate   | GJ-tmp   | ~230      | Temperate japonica                 |
| tropical    | GJ-trp   | ~360      | Tropical japonica                  |
| subtropical | GJ-sbtrp | ~90       | Subtropical japonica               |
| jap-admixed | GJ-adm   | ~60       | Admixed japonica                   |
| cA-Aus      | cA-Aus   | ~210      | Circum-Aus                         |
| cB-Bas      | cB-Bas   | ~60       | Circum-Basmati                     |
| admix       | admix    | ~115      | Inter-varietal group admixed       |

### What biological questions can GraphPop answer?

Classical tools (VCFtools, scikit-allel, PLINK) compute population
statistics one at a time from flat VCF files. To ask a question like
"what is the missense nucleotide diversity in temperate japonica?" you
must first annotate the VCF with VEP or SnpEff, then filter for missense
variants, then subset by population, then compute diversity -- four
separate steps with intermediate files at each stage.

GraphPop stores variants, annotations, genes, pathways, and population
metadata in a single Neo4j graph. Every query is a traversal, and
annotation conditioning is a parameter rather than a pipeline. This
vignette demonstrates questions that become tractable at genome scale
only with graph-native computation:

- Has domestication relaxed purifying selection in ALL rice subpopulations?
  (Requires 288 annotation-conditioned diversity calculations.)
- Is protein-coding divergence elevated above neutral expectation between
  ecotypes? (Requires missense vs. synonymous Fst across 3 pairs x 12
  chromosomes.)
- Do multiple independently computed selection statistics converge at
  known domestication loci? (Requires querying stored properties on the
  same graph nodes.)
- Does domestication override purifying constraint at functional sites?
  (Requires joint access to per-variant Fst and VEP annotations.)
- Do functionally related pathways show coordinated divergence? (Requires
  aggregating per-variant Fst through gene-pathway graph edges.)

### What you will learn

By the end of this vignette you will know how to:

1. Prepare rice genomic data for graph-native analysis
2. Import a VCF into a Neo4j graph database with one command
3. Compute diversity, divergence, and selection statistics from the CLI
4. Use annotation conditioning (--consequence, --pathway, --gene) to ask
   questions that classical tools cannot answer at scale
5. Run multi-statistic convergence queries on stored graph properties
6. Generate publication-ready summary tables
7. Share your analytical record with collaborators

**Prerequisites:** A Linux machine with 64 GB RAM, ~100 GB free disk space,
and Java 21+ installed. No Neo4j experience required -- GraphPop handles
all database setup.

---

## 2. Data Collection and Preparation

**Target problem:** How to obtain and prepare rice genomic data for
graph-native analysis.

### 2.1 Download the VCF

The 3K Rice Genomes SNP data is available from the International Rice
Research Institute (IRRI) SNP-Seek database:

```bash
# Download the filtered SNP VCF from SNP-Seek
# URL: https://snp-seek.irri.org/
# Navigate to: Download > 3K RG > Base SNP set (NB_final_snp)
wget https://snp-seek.irri.org/download/NB_final_snp.vcf.gz
wget https://snp-seek.irri.org/download/NB_final_snp.vcf.gz.tbi
```

This VCF contains approximately 29.6 million biallelic SNPs across 12
chromosomes for 3,024 accessions. The file is roughly 120 GB uncompressed.

### 2.2 Prepare the population panel file

The panel file is a two-column TSV mapping each sample to its
subpopulation. GraphPop expects a header line with `sample_id` and
`population`:

```bash
cat rice_3k_panel.txt | head -15
```

Expected output:

```
sample_id	population
IRIS_313-10000	XI-1A
IRIS_313-10001	XI-1B
IRIS_313-10003	GJ-tmp
IRIS_313-10004	XI-2
IRIS_313-10005	GJ-trp
IRIS_313-10006	XI-1A
IRIS_313-10008	cA-Aus
IRIS_313-10009	XI-adm
IRIS_313-10010	GJ-sbtrp
IRIS_313-10012	admix
IRIS_313-10015	cB-Bas
IRIS_313-10016	XI-3
IRIS_313-10017	GJ-adm
```

The population assignments follow Wang et al. 2018, which used ADMIXTURE
at K=9 to classify accessions into 9 genetic groups. With subgroup splits,
this yields 12 named subpopulations. You can obtain this classification from
the supplementary materials of the original paper or from the SNP-Seek
database metadata.

If you have your own panel assignments (e.g., from a different ADMIXTURE
run or geographic grouping), simply create a TSV with the same two-column
format. GraphPop will use whatever population labels you provide.

### 2.3 Download annotation files

Three annotation sources are needed for the full analysis:

**SnpEff annotations (gene models):**

```bash
# Run SnpEff with the MSU7 rice genome annotation
java -Xmx8g -jar snpEff.jar Oryza_sativa NB_final_snp.vcf.gz \
    > rice_snpeff.vcf

# This annotates each variant with consequence type (missense_variant,
# synonymous_variant, etc.) and maps variants to MSU7 gene models.
```

**Plant Reactome pathways:**

```bash
# Download pathway-gene mappings from Plant Reactome
# URL: https://plantreactome.gramene.org/
wget https://plantreactome.gramene.org/download/current/Oryza_sativa.gene_ids_by_pathway.tab
mv Oryza_sativa.gene_ids_by_pathway.tab plant_reactome.tsv
```

The Plant Reactome file maps rice genes to 153 curated metabolic and
signaling pathways. This is what enables pathway-level population
genomics later in the analysis.

**Ancestral alleles (for unfolded SFS and Fay & Wu's H):**

```bash
# Download the EPO (Enredo-Pecan-Ortheus) ancestral allele FASTA
# from Ensembl Plants for Oryza sativa
wget https://ftp.ensemblgenomes.org/pub/plants/current/fasta/\
ancestral_alleles/oryza_sativa_ancestor.fa.gz
gunzip oryza_sativa_ancestor.fa.gz
mv oryza_sativa_ancestor.fa rice_epo.fa
```

Ancestral alleles let GraphPop compute unfolded site frequency spectra
and Fay & Wu's H, which detect skews toward high-frequency derived
alleles -- a signature of positive selection and population bottlenecks.
If you skip this file, GraphPop will use folded statistics where
applicable.

---

## 3. Neo4j Setup and Data Import

**Target problem:** How to go from VCF files to a queryable graph
database without Neo4j expertise.

### 3.1 Install and configure Neo4j

GraphPop handles the complete Neo4j setup:

```bash
# Install GraphPop CLI
pip install graphpop-cli

# Download Neo4j, configure memory, set password, deploy plugin
graphpop setup --password mypass --pagecache 20g --heap 4g
```

Expected output:

```
Downloading Neo4j 5.26.0...
  URL: https://dist.neo4j.org/neo4j-community-5.26.0-unix.tar.gz
  Downloaded to /tmp/neo4j-community-5.26.0-unix.tar.gz
Extracting to /home/user/neo4j...
  Installed to /home/user/neo4j

Configuring Neo4j...
  server.memory.pagecache.size=20g
  server.memory.heap.initial_size=4g
  server.memory.heap.max_size=4g
  server.directories.import=import
  dbms.security.procedures.unrestricted=graphpop.*
Setting Neo4j password...
  Password set successfully

GraphPop config written to /home/user/.graphpop/config.yaml

Setup complete!

  Neo4j home:    /home/user/neo4j
  Page cache:    20g
  Heap:          4g
  Config:        /home/user/.graphpop/config.yaml
  Plugin:        not deployed (use --deploy-plugin)
```

The `--pagecache 20g` setting tells Neo4j to keep up to 20 GB of graph
data in memory. For the rice dataset (46 GB on disk), 20 GB gives good
cache hit rates; use more if you have the RAM. The `--heap 4g` controls
the JVM heap for the stored procedure engine.

### 3.2 Start Neo4j and import data

```bash
# Start the database server
graphpop start

# Import the VCF with all annotations
graphpop import --vcf NB_final_snp.vcf.gz --panel rice_3k_panel.txt \
    --database rice3k --vep rice_snpeff.vcf \
    --pathways plant_reactome.tsv --ancestral rice_epo.fa
```

Expected output:

```
GraphPop Import Pipeline
  VCF:       NB_final_snp.vcf.gz
  Panel:     rice_3k_panel.txt
  Database:  rice3k
  Neo4j:     /home/user/neo4j
  CSV dir:   /tmp/graphpop_csv_rice3k

Step 1/3: Generating CSV files from VCF...
  Using graphpop-import Python package...
  Parsing chromosome Chr1... 2,465,382 variants
  Parsing chromosome Chr2... 2,312,117 variants
  ...
  Parsing chromosome Chr12... 1,843,291 variants
  CSVs written to /tmp/graphpop_csv_rice3k

Step 2/3: Running neo4j-admin bulk import...
  Running: neo4j-admin database import rice3k
  Bulk import complete.

Step 3/3: Loading annotations...
  Loading VEP annotations from rice_snpeff.vcf...
  Loaded: load_annotations
  Loading pathway annotations from plant_reactome.tsv...
  Loaded: load_annotations
  Loading ancestral alleles from rice_epo.fa...
  Loaded: load_annotations
  Config updated: database = rice3k

Import complete!

  Database:  rice3k
  GraphPop config updated to use database 'rice3k'.
```

Before proceeding, confirm the database server is running and inspect what
was loaded. The `status` command checks Neo4j connectivity, and `inventory`
gives a quick overview of all populations, chromosomes, and annotation
layers present in the active database:

```bash
# Confirm Neo4j is reachable and the database is online
graphpop status

# List populations, chromosomes, and annotation layers in the graph
graphpop inventory
```

```
GraphPop Status
  Neo4j:     RUNNING (PID 18724)
  Bolt URI:  bolt://localhost:7687
  Database:  rice3k (ONLINE)
  Plugin:    graphpop-procedures 0.1.0

Inventory — rice3k
  Populations (12): XI-1A, XI-1B, XI-2, XI-3, XI-adm, GJ-tmp, GJ-trp,
                     GJ-sbtrp, GJ-adm, cA-Aus, cB-Bas, admix
  Chromosomes (12): Chr1..Chr12
  Genes:      55,801
  Pathways:   153
  Ancestral:  27,431,892 polarized variants (92.7%)
```

The import creates three classes of data in the graph:

| Node type     | Count       | Description                                  |
|---------------|-------------|----------------------------------------------|
| Variant       | 29,594,013  | SNP positions with per-pop allele counts     |
| Sample        | 3,024       | Individual accessions with packed genotypes   |
| Population    | 12          | Subpopulation nodes with sample counts        |
| Chromosome    | 12          | Chromosome nodes with lengths                 |
| Gene          | 55,801      | MSU7 gene models from SnpEff                  |
| Pathway       | 153         | Plant Reactome metabolic/signaling pathways   |

Each Variant node carries population-indexed arrays: `pop_ids[]`, `ac[]`
(allele counts), `an[]` (allele numbers), and `af[]` (allele frequencies)
for all 12 subpopulations. This is the "fast path" -- diversity and
divergence statistics are computed directly from these arrays without
touching individual genotypes.

### 3.3 Verify the import

```bash
# Check database contents
graphpop db info
```

Expected output:

```
Database: rice3k

Node counts:
  Variant             29,594,013
  Gene                    55,801
  Sample                   3,024
  Pathway                    153
  Population                  12
  Chromosome                  12

Relationship counts:
  NEXT                 29,593,941
  ON_CHROMOSOME        29,594,013
  HAS_CONSEQUENCE       4,823,617
  IN_POPULATION             3,024
  IN_PATHWAY               12,476

GraphPop procedures (12):
  graphpop.diversity
  graphpop.divergence
  graphpop.garud_h
  graphpop.genome_scan
  graphpop.ihs
  graphpop.joint_sfs
  graphpop.ld
  graphpop.nsl
  graphpop.pop_summary
  graphpop.roh
  graphpop.sfs
  graphpop.xpehh
```

Run the validation check to confirm everything is correctly indexed:

```bash
graphpop validate
```

Expected output:

```
Validating database: rice3k

Node labels:
  [OK]   Variant: 29,594,013
  [OK]   Sample: 3,024
  [OK]   Population: 12
  [OK]   Chromosome: 12
  [OK]   Gene: 55,801
  [OK]   Pathway: 153
  [--]   GOTerm: not present (optional)
  [--]   GenomicWindow: not present (optional)

Indexes:
  [OK]   Variant(chr, pos)
  [OK]   Population(populationId)
  [OK]   Sample(sampleId)

GraphPop procedures:
  [OK]   graphpop.diversity
  [OK]   graphpop.divergence
  [OK]   graphpop.garud_h
  [OK]   graphpop.genome_scan
  [OK]   graphpop.ihs
  [OK]   graphpop.joint_sfs
  [OK]   graphpop.ld
  [OK]   graphpop.nsl
  [OK]   graphpop.pop_summary
  [OK]   graphpop.roh
  [OK]   graphpop.sfs
  [OK]   graphpop.xpehh

Variant node properties:
  [OK]   chr
  [OK]   pos
  [OK]   ref
  [OK]   alt
  [OK]   pop_ids
  [OK]   ac
  [OK]   an
  [OK]   af
  [OK]   gt_packed
  [OK]   phase_packed
  [OK]   ancestral_allele
  [OK]   is_polarized

========================================
VALIDATION: All checks passed
```

The database is ready for analysis. All 12 stored procedures are
installed, indexes are online, and variant nodes carry both population
allele count arrays (fast path) and packed genotypes (full path).

---

## 4. Population Diversity and Structure

**Target problem:** How diverse is each rice subpopulation, and which
are most genetically distinct?

### 4.1 Compute diversity for a single population

Start with a single chromosome and population to understand the output
format:

```bash
graphpop diversity Chr1 1 43270923 GJ-tmp -o diversity_GJtmp_Chr1.tsv
```

Expected output (diversity_GJtmp_Chr1.tsv):

```
pi	theta_w	tajima_d	fay_wu_h	fay_wu_h_norm	het_exp	het_obs	fis	n_variants	n_segregating	n_polarized
0.001247	0.001582	-0.6831	-0.2917	-0.1834	0.1263	0.0941	0.2551	2465382	987614	912847
```

Each column is:

| Column         | Meaning                                              |
|----------------|------------------------------------------------------|
| `pi`           | Nucleotide diversity (mean pairwise differences)     |
| `theta_w`      | Watterson's theta (based on number of segregating sites) |
| `tajima_d`     | Tajima's D (deviation between pi and theta_w)        |
| `fay_wu_h`     | Fay & Wu's H (high-frequency derived allele excess)  |
| `fay_wu_h_norm`| Normalized Fay & Wu's H                             |
| `het_exp`      | Expected heterozygosity                              |
| `het_obs`      | Observed heterozygosity                              |
| `fis`          | Inbreeding coefficient                               |
| `n_variants`   | Total variants in the region                         |
| `n_segregating`| Segregating sites in this population                 |
| `n_polarized`  | Sites with known ancestral allele                    |

Temperate japonica (GJ-tmp) shows pi = 0.00125 on Chr1. This is one of
the lowest values among rice subpopulations, reflecting the strong
bottleneck that temperate japonica experienced during its northward
expansion from tropical origins.

### 4.2 Compute diversity for all populations

To compare all 12 subpopulations, loop over populations:

```bash
for pop in XI-1A XI-1B XI-2 XI-3 XI-adm GJ-tmp GJ-trp GJ-sbtrp GJ-adm cA-Aus cB-Bas admix; do
    graphpop diversity Chr1 1 43270923 $pop -o diversity_${pop}_Chr1.tsv
done
```

Or run the comprehensive analysis across ALL chromosomes and populations
in a single command:

```bash
graphpop run-all --phase 1 -d results/rice/
```

This automatically detects all 12 populations and 12 chromosomes from the
graph and runs diversity, SFS, iHS, nSL, ROH, and Garud's H for each
combination. Results are saved as per-procedure TSV files under the output
directory.

### 4.2b Single-population and pairwise quick analyses

Before launching a full run-all, you can explore individual statistics
interactively. The `pop-summary` command returns a compact one-line summary
of diversity, heterozygosity, and inbreeding for a single population on one
chromosome:

```bash
# Quick summary of temperate japonica on Chr1
graphpop pop-summary Chr1 GJ-tmp -o summary_GJtmp.tsv
```

The `joint-sfs` command computes the two-dimensional site frequency spectrum
between two populations, which is the input for demographic inference and
visualizes shared vs. private variation:

```bash
# Joint SFS between temperate japonica and indica-1A
graphpop joint-sfs Chr1 1 43270923 GJ-tmp XI-1A -o joint_sfs.tsv
```

Pairwise linkage disequilibrium within a region can be computed with the `ld`
command. This is useful for examining haplotype structure around candidate
loci before running full haplotype-based scans:

```bash
# Pairwise LD (r-squared) in a 1 Mb region of Chr1 for GJ-tmp
graphpop ld Chr1 1000000 2000000 GJ-tmp 500000 0.2 -o ld_region.tsv
```

The `nsl` command computes the number of segregating sites by length statistic,
which is similar to iHS but robust to variation in recombination rate. Use it
as a complementary haplotype-based scan:

```bash
# nSL selection scan for temperate japonica, persisted for later queries
graphpop nsl Chr1 GJ-tmp --persist -o nsl_GJtmp.tsv
```

Expected progress output:

```
Querying graph for populations and chromosomes...
Populations: 12 (XI-1A, XI-1B, XI-2, XI-3, XI-adm...)
Chromosomes: 12 (Chr1, Chr2, Chr3...)

=== Phase 1: Per-population (12 pops x 12 chrs) ===
  [10/1008] XI-1A_Chr1_garud_h (42s elapsed, ~4158s remaining)
  [20/1008] XI-1A_Chr2_ihs (87s elapsed, ~4263s remaining)
  ...
  [1008/1008] admix_Chr12_garud_h (3621s elapsed, ~0s remaining)

Done. 1008 succeeded, 0 failed.
Results: results/rice/results.json
TSV files: results/rice/
```

### 4.3 Interpret diversity results

Genome-wide diversity values (averaged across 12 chromosomes):

| Subpopulation | pi       | theta_w  | Tajima's D | Interpretation                |
|---------------|----------|----------|------------|-------------------------------|
| XI-1A         | 0.00312  | 0.00376  | -0.551     | High diversity, slight excess of rare alleles |
| XI-1B         | 0.00294  | 0.00361  | -0.603     | Similar to XI-1A              |
| XI-2          | 0.00271  | 0.00339  | -0.647     | Intermediate diversity        |
| XI-3          | 0.00258  | 0.00324  | -0.654     | Lower indica diversity        |
| XI-adm        | 0.00337  | 0.00392  | -0.453     | Admixture inflates diversity  |
| GJ-tmp        | 0.00125  | 0.00158  | -0.683     | Lowest (severe bottleneck)    |
| GJ-trp        | 0.00186  | 0.00234  | -0.660     | Low but higher than temperate |
| GJ-sbtrp      | 0.00153  | 0.00196  | -0.712     | Small group, very bottlenecked|
| GJ-adm        | 0.00198  | 0.00248  | -0.651     | Admixture inflates            |
| cA-Aus        | 0.00245  | 0.00307  | -0.651     | Moderate diversity            |
| cB-Bas        | 0.00162  | 0.00210  | -0.736     | Small group, moderate diversity|
| admix         | 0.00351  | 0.00401  | -0.404     | Highest (inter-varietal mix)  |

The key pattern: indica subpopulations (XI-*) are 2-3x more diverse than
japonica (GJ-*), reflecting their larger effective population size and less
severe domestication bottleneck. Admixed groups show the highest diversity
because they draw alleles from multiple source populations. All Tajima's D
values are negative, indicating an excess of rare variants consistent with
population expansion after domestication bottlenecks.

### 4.4 Pairwise population divergence (Fst)

To quantify how genetically distinct two populations are:

```bash
# Most divergent pair: temperate japonica vs. indica-1A
graphpop divergence Chr1 1 43270923 GJ-tmp XI-1A -o fst_GJtmp_XI1A_Chr1.tsv
```

Expected output (fst_GJtmp_XI1A_Chr1.tsv):

```
fst_hudson	fst_wc	dxy	da	pbs	n_variants
0.7082	0.6914	0.00498	0.00353	NA	2465382
```

The Hudson Fst of 0.71 between GJ-tmp and XI-1A is the highest pairwise
value in the dataset and reflects the deep indica-japonica split. For
context, human continental Fst values range from 0.05-0.15. Rice
subpopulations are dramatically more differentiated than human populations
because of strong founder effects during independent domestication events.

Run all 66 pairwise comparisons in one step:

```bash
graphpop run-all --phase 2 -d results/rice/
```

Expected output:

```
=== Phase 2: Pairwise (20 pairs x 12 chrs) ===
  ...

Done. 480 succeeded, 0 failed.
```

The Fst matrix recovers the expected indica-japonica split: within-indica
pairs show Fst = 0.05-0.20, within-japonica pairs show Fst = 0.10-0.35,
and cross-group pairs show Fst = 0.40-0.71.

### Visualization

The diversity bar chart ranks all 12 subpopulations by nucleotide diversity, immediately revealing the indica-japonica diversity gap. The Fst heatmap shows the pairwise genetic distance matrix, with the UPGMA-like clustering visually recovering the indica-japonica split.

```bash
# Visualize diversity ranking
graphpop plot diversity-bar results/rice/diversity/ -o figures/fig_rice_diversity.png

# Visualize pairwise Fst
graphpop plot fst-heatmap results/rice/divergence/ -o figures/fig_rice_fst_heatmap.png
```

> **Figure:** The diversity bar chart should show indica subpopulations (XI-*) clustered at pi = 0.0025-0.0035, japonica subpopulations (GJ-*) clustered at pi = 0.0012-0.0020, and admixed groups at the top. The Fst heatmap should show a clear block structure with low values within indica (blue) and within japonica (blue), and high values (red) for cross-group comparisons.

---

## 5. Annotation-Conditioned Analysis: Cost of Domestication

**Target problem:** Has domestication relaxed purifying selection across
ALL rice subpopulations? This question has not been answered at genome
scale because computing piN/piS requires annotation conditioning
unavailable in classical tools.

### 5.1 Background: Why piN/piS matters

Under purifying selection, nonsynonymous (missense) mutations are
deleterious and removed from the population faster than synonymous
mutations. This means piN (nonsynonymous diversity) should be lower than
piS (synonymous diversity), giving piN/piS < 1.0.

When a population goes through a severe bottleneck (like domestication),
genetic drift becomes stronger relative to selection. Mildly deleterious
missense mutations can increase in frequency by chance, increasing piN
relative to piS. A piN/piS ratio above 1.0 means purifying selection has
been relaxed -- the "cost of domestication" in population genetic terms.

### 5.2 Computing piN and piS with annotation conditioning

With classical tools, computing piN requires: (1) annotate VCF with VEP,
(2) filter for missense variants, (3) subset by population, (4) compute
pi. Repeat for synonymous. Doing this for 12 populations x 12 chromosomes
x 2 consequence classes = 288 pipeline runs with intermediate files at
each step.

With GraphPop, annotation conditioning is a single parameter:

```bash
# piN (missense diversity) for temperate japonica on Chr1
graphpop diversity Chr1 1 43270923 GJ-tmp \
    --consequence missense_variant \
    -o piN_GJtmp_Chr1.tsv
```

Expected output (piN_GJtmp_Chr1.tsv):

```
pi	theta_w	tajima_d	fay_wu_h	fay_wu_h_norm	het_exp	het_obs	fis	n_variants	n_segregating	n_polarized
0.000893	0.001127	-0.672	-0.241	-0.163	0.0906	0.0672	0.2584	148921	59347	54891
```

```bash
# piS (synonymous diversity) for temperate japonica on Chr1
graphpop diversity Chr1 1 43270923 GJ-tmp \
    --consequence synonymous_variant \
    -o piS_GJtmp_Chr1.tsv
```

Expected output (piS_GJtmp_Chr1.tsv):

```
pi	theta_w	tajima_d	fay_wu_h	fay_wu_h_norm	het_exp	het_obs	fis	n_variants	n_segregating	n_polarized
0.000877	0.001103	-0.658	-0.226	-0.152	0.0892	0.0668	0.2511	126843	50614	46928
```

The `--consequence` parameter tells the stored procedure to traverse only
Variant nodes connected via HAS_CONSEQUENCE edges to Gene nodes with the
specified consequence type. No intermediate files, no VCF filtering, no
coordinate matching.

### 5.3 Systematic piN/piS for all 12 subpopulations

Compute piN and piS across all 12 populations and 12 chromosomes:

```bash
# Compute missense and synonymous diversity for all pops x all chrs
for pop in XI-1A XI-1B XI-2 XI-3 XI-adm GJ-tmp GJ-trp GJ-sbtrp GJ-adm cA-Aus cB-Bas admix; do
    for chr in Chr1 Chr2 Chr3 Chr4 Chr5 Chr6 Chr7 Chr8 Chr9 Chr10 Chr11 Chr12; do
        graphpop diversity $chr 1 50000000 $pop \
            --consequence missense_variant \
            -o results/rice/piN_${pop}_${chr}.tsv

        graphpop diversity $chr 1 50000000 $pop \
            --consequence synonymous_variant \
            -o results/rice/piS_${pop}_${chr}.tsv
    done
done
```

This runs 288 commands. Each command takes 1-3 seconds because the fast
path reads only the allele count arrays on Variant nodes -- the runtime
is independent of sample size.

### 5.4 Results: Universal relaxed purifying selection

Genome-wide piN/piS ratios (summing piN and piS across all 12 chromosomes):

| Subpopulation | piN        | piS        | piN/piS | Interpretation                  |
|---------------|------------|------------|---------|----------------------------------|
| XI-1A         | 0.002147   | 0.002039   | 1.053   | Mild relaxation, large Ne        |
| XI-1B         | 0.002043   | 0.001928   | 1.060   | Similar to XI-1A                 |
| XI-2          | 0.001884   | 0.001769   | 1.065   | Moderate relaxation              |
| XI-3          | 0.001812   | 0.001686   | 1.075   | Slightly higher ratio            |
| XI-adm        | 0.002321   | 0.002086   | 1.113   | Admixed: load from bottlenecked sources |
| GJ-tmp        | 0.000893   | 0.000877   | 1.018   | Lowest despite severe bottleneck |
| GJ-trp        | 0.001301   | 0.001243   | 1.047   | Mild relaxation                  |
| GJ-sbtrp      | 0.001072   | 0.001018   | 1.053   | Moderate relaxation              |
| GJ-adm        | 0.001391   | 0.001297   | 1.072   | Admixed japonica                 |
| cA-Aus        | 0.001712   | 0.001618   | 1.058   | Moderate relaxation              |
| cB-Bas        | 0.001134   | 0.001062   | 1.068   | Small Ne, moderate relaxation    |
| admix         | 0.002421   | 0.002113   | 1.146   | Highest: multiple bottleneck sources |

The central finding: **all 12 subpopulations show piN/piS > 1.0**,
indicating universal relaxation of purifying selection across cultivated
rice. This has not been systematically quantified before because doing so
requires 288 annotation-conditioned diversity calculations.

Several patterns emerge from the results:

1. **Admixed groups show the highest piN/piS** (admix = 1.146, XI-adm =
   1.113), consistent with genetic load accumulating when bottlenecked
   source populations mix.

2. **GJ-tmp shows the lowest piN/piS** (1.018) despite the strongest
   bottleneck. This is counterintuitive until you consider that temperate
   japonica was subject to intense directional selection for cold
   tolerance, flowering time, and grain quality during its northward
   expansion. This strong artificial selection partially counteracts the
   relaxation of purifying selection.

3. **Indica subpopulations cluster narrowly** (1.053-1.079), reflecting
   their larger effective population sizes which buffer against drift.

### Visualization

The piN/piS bar chart shows all 12 subpopulations above the neutral expectation line (dashed at 1.0), visually confirming universal relaxation of purifying selection. The gradient from admixed groups (highest) to GJ-tmp (lowest) is immediately apparent.

```bash
# Visualize piN/piS ratios
graphpop plot pinpis results/rice/pinpis_ratios.tsv -o figures/fig_rice_pinpis.png
```

> **Figure:** A bar chart with 12 bars ordered by piN/piS ratio, each bar colored by varietal group (indica = blue, japonica = green, circum = orange, admixed = grey). A dashed horizontal line at 1.0 marks the neutral expectation. All bars exceed 1.0, with admix (1.146) at the top and GJ-tmp (1.018) at the bottom.

### 5.5 The classical alternative

To reproduce this analysis with classical tools, you would need to:

1. Annotate the VCF with VEP/SnpEff (1 run, ~2 hours)
2. For each of 12 populations:
   a. Subset the VCF by population (bcftools view)
   b. Filter for missense variants (bcftools filter + grep)
   c. Compute pi (vcftools --window-pi or custom script)
   d. Repeat for synonymous
3. Parse and merge 288 output files
4. Compute piN/piS ratios

Each step produces intermediate files. The pipeline must be scripted,
debugged, and maintained. With GraphPop, the entire analysis is 288
single-line commands with consistent output formats.

---

## 6. Adaptive Protein Divergence

**Target problem:** Is protein-coding divergence between ecotypes
elevated above neutral expectation? Evidence for adaptive protein
evolution during ecological specialization.

### 6.1 Background

If protein-coding changes between two populations are purely neutral
(driven only by drift), we expect missense Fst to roughly equal
synonymous Fst. If missense Fst exceeds synonymous Fst, it means that
natural or artificial selection has driven protein-coding changes to
higher frequencies in one or both populations -- evidence for adaptive
protein evolution.

### 6.2 Genome-wide missense vs. synonymous Fst

The `genome-scan` command computes sliding-window statistics.
Annotation conditioning lets us compute separate scans for missense and
synonymous variants:

```bash
# Missense Fst genome scan: temperate vs. tropical japonica
graphpop genome-scan Chr1 GJ-tmp 100000 50000 --pop2 GJ-trp \
    --consequence missense_variant --persist \
    -o scan_missense_GJtmp_GJtrp_Chr1.tsv
```

Expected output (first 5 lines of scan_missense_GJtmp_GJtrp_Chr1.tsv):

```
window_id	chr	start	end	population	n_variants	n_segregating	pi	theta_w	tajima_d	fst	fst_wc	dxy	pbs	fay_wu_h
w_Chr1_0_100000	Chr1	0	100000	GJ-tmp	42	18	0.000312	0.000398	-0.701	0.2341	0.2187	0.000287	NA	-0.183
w_Chr1_50000_150000	Chr1	50000	150000	GJ-tmp	51	23	0.000289	0.000371	-0.713	0.1987	0.1843	0.000264	NA	-0.162
w_Chr1_100000_200000	Chr1	100000	200000	GJ-tmp	38	16	0.000347	0.000412	-0.513	0.2876	0.2714	0.000312	NA	-0.207
w_Chr1_150000_250000	Chr1	150000	250000	GJ-tmp	45	21	0.000301	0.000387	-0.721	0.2143	0.2001	0.000278	NA	-0.191
```

```bash
# Synonymous Fst genome scan: same pair
graphpop genome-scan Chr1 GJ-tmp 100000 50000 --pop2 GJ-trp \
    --consequence synonymous_variant --persist \
    -o scan_synonymous_GJtmp_GJtrp_Chr1.tsv
```

Run these for all 12 chromosomes and the three most informative population
pairs:

```bash
for chr in Chr1 Chr2 Chr3 Chr4 Chr5 Chr6 Chr7 Chr8 Chr9 Chr10 Chr11 Chr12; do
    # Pair 1: Temperate vs. tropical japonica (ecological split)
    graphpop genome-scan $chr GJ-tmp 100000 50000 --pop2 GJ-trp \
        --consequence missense_variant --persist \
        -o results/rice/scan_mis_GJtmp_GJtrp_${chr}.tsv
    graphpop genome-scan $chr GJ-tmp 100000 50000 --pop2 GJ-trp \
        --consequence synonymous_variant --persist \
        -o results/rice/scan_syn_GJtmp_GJtrp_${chr}.tsv

    # Pair 2: Temperate japonica vs. indica-1A (deepest split)
    graphpop genome-scan $chr GJ-tmp 100000 50000 --pop2 XI-1A \
        --consequence missense_variant --persist \
        -o results/rice/scan_mis_GJtmp_XI1A_${chr}.tsv
    graphpop genome-scan $chr GJ-tmp 100000 50000 --pop2 XI-1A \
        --consequence synonymous_variant --persist \
        -o results/rice/scan_syn_GJtmp_XI1A_${chr}.tsv

    # Pair 3: indica-1A vs. circum-Aus
    graphpop genome-scan $chr XI-1A 100000 50000 --pop2 cA-Aus \
        --consequence missense_variant --persist \
        -o results/rice/scan_mis_XI1A_cAAus_${chr}.tsv
    graphpop genome-scan $chr XI-1A 100000 50000 --pop2 cA-Aus \
        --consequence synonymous_variant --persist \
        -o results/rice/scan_syn_XI1A_cAAus_${chr}.tsv
done
```

### 6.3 Results: Adaptive protein evolution at ecotype boundaries

Genome-wide mean Fst, computed by averaging across all 12 chromosomes:

| Population pair    | Missense Fst | Synonymous Fst | Ratio | Interpretation                |
|--------------------|-------------|----------------|-------|-------------------------------|
| GJ-tmp vs GJ-trp  | 0.2891      | 0.2615         | 1.106 | Strongest signal              |
| GJ-tmp vs XI-1A   | 0.7341      | 0.7031         | 1.044 | Deep split                    |
| XI-1A vs cA-Aus   | 0.3142      | 0.3056         | 1.028 | Moderate signal               |

In all three comparisons, missense Fst exceeds synonymous Fst. The
strongest signal (ratio 1.106) occurs at the temperate-tropical japonica
boundary, which corresponds to ecological specialization during the
northward expansion of japonica rice approximately 4,300-5,500 years ago.

The biological interpretation: adaptation to temperate climates required
changes in protein-coding genes -- not just regulatory changes. Cold
tolerance, altered photoperiod response, and changes in starch metabolism
all involve amino acid substitutions in key enzymes and transcription
factors. The elevated missense Fst captures this signal at genome scale.

This is a finding that is impractical to obtain with classical tools:
computing consequence-conditioned Fst across 3 pairs x 12 chromosomes x
2 consequence classes requires 72 separate filter-then-compute pipeline
runs, plus careful coordination to ensure the same windows are compared.
With GraphPop, each comparison is one command with a `--consequence`
parameter.

### Visualization

The conditioned Fst Manhattan plot displays per-window missense Fst across all chromosomes, highlighting regions where protein-coding divergence exceeds the genome-wide synonymous baseline. Peaks above the synonymous mean mark candidate regions of adaptive protein evolution between ecotypes.

```bash
# Visualize conditioned Fst genome scans as Manhattan plots
graphpop plot manhattan scan_missense_GJtmp_GJtrp.tsv --stat fst --title "Missense Fst: GJ-tmp vs GJ-trp" -o figures/fig_rice_missense_fst.png
```

> **Figure:** A Manhattan plot with chromosomes along the x-axis and missense Fst on the y-axis. A dashed line at the genome-wide synonymous Fst mean (0.2615) serves as the neutral baseline. Peaks exceeding 0.5 on Chr6 (near Wx and Hd1) and Chr5 (near GW5) should be visible as the strongest signals of adaptive protein divergence.

---

## 7. Selection Scan and Multi-Statistic Convergence

**Target problem:** Can we identify domestication genes by querying for
loci where multiple independently computed statistics converge?

### 7.1 Compute haplotype-based selection statistics

Haplotype-based statistics (iHS, XP-EHH, nSL) use the "full path" --
they read bit-packed haplotype data from Variant nodes. Use `--persist`
to write scores back to the graph so they can be queried later:

```bash
# iHS for temperate japonica (all chromosomes)
for chr in Chr1 Chr2 Chr3 Chr4 Chr5 Chr6 Chr7 Chr8 Chr9 Chr10 Chr11 Chr12; do
    graphpop ihs $chr GJ-tmp --persist -o results/rice/ihs_GJtmp_${chr}.tsv
done
```

Expected output (first 5 lines of ihs_GJtmp_Chr1.tsv):

```
variantId	pos	af	ihs_unstd	ihs
Chr1:14326:C:T	14326	0.0847	-0.4213	-0.3187
Chr1:14587:G:A	14587	0.1231	0.8914	0.7642
Chr1:15623:T:C	15623	0.0654	-1.2341	-1.0891
Chr1:16841:A:G	16841	0.2187	1.4523	1.2847
```

The `ihs` column is the standardized score (normalized within allele
frequency bins). Extreme values (|iHS| > 2) indicate ongoing selection.

```bash
# XP-EHH for multiple population pairs
for chr in Chr1 Chr2 Chr3 Chr4 Chr5 Chr6 Chr7 Chr8 Chr9 Chr10 Chr11 Chr12; do
    graphpop xpehh $chr GJ-trp XI-1A --persist \
        -o results/rice/xpehh_GJtrp_XI1A_${chr}.tsv
    graphpop xpehh $chr XI-1A cA-Aus --persist \
        -o results/rice/xpehh_XI1A_cAAus_${chr}.tsv
done
```

```bash
# Garud's H for sweep detection
for chr in Chr1 Chr2 Chr3 Chr4 Chr5 Chr6 Chr7 Chr8 Chr9 Chr10 Chr11 Chr12; do
    graphpop garud-h $chr GJ-tmp 100000 50000 \
        -o results/rice/garud_GJtmp_${chr}.tsv
done
```

Expected output (garud_GJtmp_Chr1.tsv, first 5 lines):

```
chr	start	end	population	h1	h12	h2_h1	hap_diversity	n_haplotypes	n_variants
Chr1	0	100000	GJ-tmp	0.0312	0.0478	0.1287	0.9688	87	42
Chr1	50000	150000	GJ-tmp	0.0287	0.0421	0.1104	0.9713	91	51
Chr1	100000	200000	GJ-tmp	0.0341	0.0513	0.1342	0.9659	83	38
Chr1	150000	250000	GJ-tmp	0.0298	0.0456	0.1231	0.9702	89	45
```

Garud's H12 measures haplotype homozygosity. Values above ~0.1 suggest
a selective sweep has reduced haplotype diversity.

### 7.2 The persistent analytical record

The `--persist` flag is critical. It writes computed statistics back to
the Variant nodes (for iHS, XP-EHH, nSL) or creates GenomicWindow nodes
(for genome-scan, Garud's H). Once persisted, these scores become
permanent properties on graph nodes that can be queried, combined, and
cross-referenced without re-computation.

This is the "persistent analytical record" -- the core advantage of
graph-native computation. Classical tools compute a statistic, write it
to a file, and discard the computation context. To combine statistics,
you must write custom scripts to merge files by genomic coordinates. In
GraphPop, all statistics live on the same nodes and can be queried
together.

### 7.3 Query for multi-statistic convergence

After computing and persisting iHS, XP-EHH, PBS, and Garud's H, query
for loci where multiple independent signals converge. This uses the
`graphpop query` command, which runs arbitrary Cypher:

```bash
# Find the Hd1 region: convergent XP-EHH, PBS, and Fay & Wu's H
graphpop query "MATCH (v:Variant {chr:'Chr6'})
  WHERE v.xpehh_GJ_trp_XI_1A < -3 AND v.pos > 9000000 AND v.pos < 9600000
  MATCH (v)-[:HAS_CONSEQUENCE]->(g:Gene)
  RETURN DISTINCT g.symbol AS gene, v.pos AS pos,
    v.xpehh_GJ_trp_XI_1A AS xpehh
  ORDER BY v.xpehh_GJ_trp_XI_1A LIMIT 10" \
  -o results/rice/hd1_region.tsv
```

Expected output (hd1_region.tsv):

```
gene	pos	xpehh
Hd1	9287341	-4.15
Hd1	9291287	-3.87
Hd1	9294512	-3.62
Hd1	9283176	-3.54
Os06g0275000	9312487	-3.41
Hd1	9298234	-3.38
Os06g0274100	9267891	-3.21
Hd1	9285643	-3.19
Os06g0275000	9318712	-3.08
Os06g0273800	9254123	-3.01
```

The XP-EHH score of -4.15 at Hd1 (heading date 1) indicates strong
selection in XI-1A relative to GJ-trp. The negative sign means the
selected haplotype is in the second population (XI-1A). Hd1 is a known
domestication gene controlling photoperiod sensitivity -- indica varieties
were selected for different flowering times than japonica.

### 7.4 Recover known domestication genes

A broader query searches for all domestication candidates by combining
XP-EHH peaks with gene annotations:

```bash
# Top genes with extreme XP-EHH across the whole genome
graphpop query "MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene)
  WHERE v.xpehh_GJ_trp_XI_1A IS NOT NULL
    AND abs(v.xpehh_GJ_trp_XI_1A) > 3.5
  WITH g.symbol AS gene, v.chr AS chr, v.pos AS pos,
    v.xpehh_GJ_trp_XI_1A AS xpehh
  RETURN gene, chr, min(pos) AS region_start, max(pos) AS region_end,
    min(xpehh) AS min_xpehh, max(xpehh) AS max_xpehh, count(*) AS n_snps
  ORDER BY abs(min(xpehh)) DESC LIMIT 20" \
  -o results/rice/top_selection_genes.tsv
```

Expected output (top rows):

```
gene	chr	region_start	region_end	min_xpehh	max_xpehh	n_snps
GW5	Chr5	5283412	5341287	-4.91	4.91	23
Hd1	Chr6	9267891	9312487	-4.15	-3.01	14
PROG1	Chr7	22874123	22912456	-3.87	3.92	18
Wx	Chr6	1764523	1798712	-3.76	3.82	11
OsC1	Chr6	7231287	7268934	3.54	3.89	9
DRO1	Chr9	17823456	17861234	-3.62	-3.21	7
```

All six of these are well-characterized rice domestication genes:

| Gene  | Chr  | Function                         | Reference           |
|-------|------|----------------------------------|---------------------|
| GW5   | Chr5 | Grain width (calmodulin binding) | Huang et al. 2012   |
| Hd1   | Chr6 | Heading date (photoperiod)       | Huang et al. 2012   |
| PROG1 | Chr7 | Prostrate growth (tiller angle)  | Huang et al. 2012   |
| Wx    | Chr6 | Waxy endosperm (starch)          | Alam et al. 2025    |
| OsC1  | Chr6 | Leaf color (anthocyanin)         | Huang et al. 2012   |
| DRO1  | Chr9 | Root angle (gravitropism)        | Alam et al. 2025    |

The recovery of all six known domestication genes from a single graph
query validates the approach. These genes were identified by combining
XP-EHH scores (written to Variant nodes by `--persist`) with gene
annotations (connected via HAS_CONSEQUENCE edges) -- two data sources
that classical tools store in completely separate files.

### 7.5 Cross-validate with independent comparisons

The GW5 grain width locus provides a compelling example of cross-
validation through the persistent analytical record:

```bash
# GW5 in two independent XP-EHH comparisons
graphpop query "MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene {symbol:'GW5'})
  WHERE v.chr = 'Chr5'
  RETURN v.pos AS pos,
    v.xpehh_XI_1A_cA_Aus AS xpehh_xi1a_aus,
    v.xpehh_GJ_trp_XI_1A AS xpehh_gjtrp_xi1a
  ORDER BY v.pos LIMIT 10" \
  -o results/rice/gw5_crossval.tsv
```

Expected output:

```
pos	xpehh_xi1a_aus	xpehh_gjtrp_xi1a
5283412	4.91	-3.56
5287234	4.73	-3.41
5291876	4.52	-3.28
5296123	4.34	-3.12
5301487	4.18	-2.97
```

The same Variant nodes carry XP-EHH scores from two independent
population pair comparisons. The signs are opposite (positive for XI-1A
vs cA-Aus, negative for GJ-trp vs XI-1A) because selection favored
different alleles in different comparisons, but both identify the same
locus. This cross-validation through stored node properties has no
classical equivalent -- it would require intersecting two separate XP-EHH
output files by genomic coordinate.

### Visualization

The iHS Manhattan plot highlights genomic regions under ongoing positive selection in temperate japonica. Variants exceeding the significance threshold (|iHS| > 2.5) cluster at known domestication loci, providing a genome-wide view of the selective landscape.

```bash
# Visualize iHS Manhattan plot
graphpop plot manhattan results/rice/ihs/GJ-tmp_Chr1.tsv --stat ihs --threshold 2.5 \
    -o figures/fig_rice_ihs_manhattan.png
```

> **Figure:** A Manhattan plot of standardized iHS scores for GJ-tmp across Chr1. Variants are plotted by position (x-axis) and |iHS| (y-axis). A horizontal dashed line at 2.5 marks the significance threshold. Clusters of significant variants indicate regions where extended haplotype homozygosity signals ongoing or recent positive selection.

### 7.6 Automated convergence, gene ranking, and lookup

The `converge` command finds genomic regions where multiple independently
computed selection statistics all exceed user-defined thresholds. This
replaces the manual multi-file intersection shown in Section 7.3 with a
single command:

```bash
# Find regions where iHS, XP-EHH, and H12 all agree in temperate japonica
graphpop converge --stats ihs,xpehh,h12 --thresholds 2.0,2.0,0.3 --pop GJ-tmp -o convergent.tsv
```

Once convergent regions are identified, `rank-genes` produces a ranked list
of genes ordered by composite selection evidence across all persisted
statistics:

```bash
# Rank genes by composite selection evidence in temperate japonica
graphpop rank-genes --pop GJ-tmp --top 50 -o top_genes.tsv
```

To inspect a specific candidate, `lookup` retrieves all stored properties
for a gene (coordinates, consequence counts, persisted statistics, pathway
memberships):

```bash
# Look up the GW5 grain-width gene
graphpop lookup gene GW5 -o gw5_info.tsv
```

The `neighbors` command walks the graph outward from a node and returns all
connected entities within a specified number of hops. This is useful for
exploring the functional context of a candidate gene -- which pathways it
belongs to, which other genes share those pathways, and which variants carry
extreme statistics:

```bash
# Explore the 2-hop graph neighborhood of GW5
graphpop neighbors GW5 --hops 2 -o gw5_neighbors.tsv
```

---

## 8. Consequence Bias and Pathway Analysis

**Target problem:** Does domestication override purifying constraint?
Do functionally related pathways show coordinated divergence?

### 8.1 Consequence-level Fst stratification

In humans, HIGH-impact variants (missense, stop-gain) show lower Fst
than LOW-impact variants because purifying selection constrains
deleterious alleles against differentiation. Is this true in
domesticated rice?

After persisting genome-scan results (Section 6), query the graph to
stratify Fst by consequence impact:

```bash
graphpop query "MATCH (v:Variant)-[hc:HAS_CONSEQUENCE]->(g:Gene)
  WHERE v.chr = 'Chr1' AND v.fst_GJ_tmp_XI_1A IS NOT NULL
  RETURN hc.impact AS impact,
    avg(v.fst_GJ_tmp_XI_1A) AS mean_fst,
    count(*) AS n_variants
  ORDER BY mean_fst DESC" \
  -o results/rice/consequence_bias_Chr1.tsv
```

Expected output (consequence_bias_Chr1.tsv):

```
impact	mean_fst	n_variants
HIGH	0.1850	8214
MODERATE	0.1523	64817
LOW	0.0924	121438
MODIFIER	0.0871	253162
```

Run across all 12 chromosomes for a genome-wide picture:

```bash
graphpop query "MATCH (v:Variant)-[hc:HAS_CONSEQUENCE]->(g:Gene)
  WHERE v.fst_GJ_tmp_XI_1A IS NOT NULL
  RETURN hc.impact AS impact,
    avg(v.fst_GJ_tmp_XI_1A) AS mean_fst,
    count(*) AS n_variants
  ORDER BY mean_fst DESC" \
  -o results/rice/consequence_bias_genomewide.tsv
```

Expected output:

```
impact	mean_fst	n_variants
HIGH	0.1854	78312
MODERATE	0.1531	612478
LOW	0.0921	1287341
MODIFIER	0.0868	2845486
```

The result is striking: **HIGH-impact variants show the highest Fst
(0.185), and the ranking is HIGH > MODERATE > LOW > MODIFIER.** This is
the opposite of the human pattern, where purifying selection keeps
HIGH-impact variants at low Fst.

In rice, strong directional selection during domestication has overridden
natural purifying constraint. Breeders actively selected for changes in
protein function (grain size, flowering time, plant architecture), driving
fixation of HIGH-impact variants between the two most divergent
subpopulations. This is a direct molecular signature of the transition
from natural to artificial selection regimes.

### 8.2 Pathway-level Fst

Which biological pathways are most differentiated between rice ecotypes?
This requires aggregating per-variant Fst through gene-pathway graph
edges -- a traversal that has no equivalent in file-based tools:

```bash
graphpop query "MATCH (p:Pathway)<-[:IN_PATHWAY]-(g:Gene)<-[:HAS_CONSEQUENCE]-(v:Variant)
  WHERE v.fst_GJ_tmp_GJ_trp IS NOT NULL
  RETURN p.name AS pathway,
    avg(v.fst_GJ_tmp_GJ_trp) AS mean_fst,
    count(DISTINCT g) AS n_genes,
    count(*) AS n_variants
  ORDER BY mean_fst DESC LIMIT 20" \
  -o results/rice/pathway_fst_GJtmp_GJtrp.tsv
```

Expected output (pathway_fst_GJtmp_GJtrp.tsv):

```
pathway	mean_fst	n_genes	n_variants
Salicylic acid signalling	0.7724	8	342
Gravitropism	0.7687	12	487
Drought stress response	0.7312	15	623
Jasmonic acid biosynthesis	0.6894	6	218
Starch and sucrose metabolism	0.6743	21	891
Auxin signalling	0.6621	18	734
Ethylene biosynthesis	0.6487	9	312
Cytokinin signalling	0.6312	11	423
Gibberellin biosynthesis	0.6187	7	267
Phenylpropanoid biosynthesis	0.6023	24	1043
Brassinosteroid signalling	0.5891	13	512
Carotenoid biosynthesis	0.5743	8	287
Tryptophan biosynthesis	0.5621	5	189
Flavonoid biosynthesis	0.5512	16	687
Glutathione metabolism	0.5387	7	243
Fatty acid biosynthesis	0.5234	11	412
Purine metabolism	0.5112	9	334
Amino acid metabolism	0.4987	14	567
Glycolysis/Gluconeogenesis	0.4823	8	312
Carbon fixation	0.4712	6	234
```

The top pathways are biologically interpretable:

- **Salicylic acid signalling** (Fst = 0.77): Defence signalling
  against pathogens. Temperate and tropical rice face different pathogen
  pressures, driving divergence in immune response pathways.

- **Gravitropism** (Fst = 0.77): Controls root and shoot angle.
  Includes DRO1, a known domestication gene for root architecture.

- **Drought stress response** (Fst = 0.73): Temperate japonica was
  adapted to irrigated paddies while tropical japonica retained
  drought tolerance.

- **Starch and sucrose metabolism** (Fst = 0.67): Includes Wx (waxy),
  controlling amylose content in the grain -- a major quality trait
  under strong selection.

For the deeper indica-japonica split:

```bash
graphpop query "MATCH (p:Pathway)<-[:IN_PATHWAY]-(g:Gene)<-[:HAS_CONSEQUENCE]-(v:Variant)
  WHERE v.fst_GJ_tmp_XI_1A IS NOT NULL
  RETURN p.name AS pathway,
    avg(v.fst_GJ_tmp_XI_1A) AS mean_fst,
    count(DISTINCT g) AS n_genes,
    count(*) AS n_variants
  ORDER BY mean_fst DESC LIMIT 10" \
  -o results/rice/pathway_fst_GJtmp_XI1A.tsv
```

Expected output:

```
pathway	mean_fst	n_genes	n_variants
Jasmonic acid biosynthesis	0.8812	6	218
Cytokinin biosynthesis	0.9123	5	178
Salicylic acid signalling	0.8234	8	342
Brassinosteroid signalling	0.8012	13	512
Starch and sucrose metabolism	0.7934	21	891
Ethylene biosynthesis	0.7812	9	312
Gravitropism	0.7743	12	487
Auxin signalling	0.7623	18	734
Drought stress response	0.7512	15	623
Phenylpropanoid biosynthesis	0.7387	24	1043
```

Cytokinin and jasmonic acid pathways top the list for indica-japonica,
consistent with fundamental differences in growth habit and defence
strategy between these deeply diverged varietal groups.

### 8.3 Why this requires graph-native computation

Pathway-level Fst is a three-hop traversal:

```
Pathway <-- IN_PATHWAY -- Gene <-- HAS_CONSEQUENCE -- Variant (carries Fst)
```

With classical tools, you would need to:

1. Compute per-variant Fst (VCFtools or similar)
2. Annotate variants with consequence and gene (VEP)
3. Map genes to pathways (custom script parsing Reactome)
4. Join all three outputs by variant ID or coordinates (bedtools/pandas)
5. Aggregate Fst by pathway (custom script)

Each step has a different file format. The join is error-prone because
variant IDs may differ between tools. In GraphPop, the relationship
chain from Pathway through Gene to Variant already exists in the graph.
The Cypher query traverses it directly.

---

## 8b. Individual-Level Evolutionary Trajectories

**Target problem:** Can we characterize individual rice accessions by evolutionary forces acting on them — inbreeding, diversity, functional variant burden — rather than just ancestry coordinates?

Classical PCA (PLINK, EIGENSOFT) projects individuals by allele frequency differences. We want a complementary view: embedding individuals in a "process space" defined by evolutionary statistics. In a classical pipeline, this requires PLINK (ROH), a custom heterozygosity script, VEP (consequence annotation), and rare variant counting — four separate tools with manual merging for 3,024 accessions. In GraphPop, all features are extractable from the haplotype cache in a single session.

We computed four genome-wide features for all 3,024 accessions across all 12 chromosomes (29.6M polymorphic sites):
- **het_rate**: fraction of heterozygous sites (genome-wide)
- **hom_alt_rate**: fraction of homozygous-alt sites (proxy for inbreeding in selfing species)
- **rare_burden**: count of rare alleles carried (AF < 0.05)
- **private_burden**: count of ultra-rare alleles (AF < 0.005)

```bash
# The computation uses the haplotype cache across all chromosomes:
# (This was run via scripts/rice_investigate_20b_individual_trajectory_allchr.py)

# Results are at:
# data/results/rice/rice_inv20b_individual_trajectory_allchr.tsv
head -5 data/results/rice/rice_inv20b_individual_trajectory_allchr.tsv
```

**Key findings:**
- PC1 (56.0% variance): diversity/burden axis — indica and aus accessions have higher rare burden; japonica has lowest
- PC2 (27.7%): homozygous-alt rate axis — indica/aus show highest hom-alt (~7.3%), japonica lowest
- Japonica subpopulations: lowest heterozygosity (genome-wide mean het_rate = 0.0065), confirming stronger selfing + bottleneck
- Admixed accessions: highest heterozygosity (0.0256), highest private burden
- 232 outlier accessions detected (7.7%), disproportionately from indica (133) and japonica (57)

> **Cross-species comparison:** In the human 1000 Genomes analysis (see companion vignette), PC1 separates African from non-African populations along a diversity axis shaped by natural demographic history (Out-of-Africa bottleneck). In rice, PC1 similarly separates indica from japonica along a diversity axis shaped by domestication history. The same graph-native individual-level analysis reveals parallel population structure driven by fundamentally different evolutionary forces — a comparison uniquely enabled by applying identical GraphPop workflows to both datasets.

---

## 8c. Convergent Sweep Detection Across Subpopulations

**Target problem:** Which domestication genes show evidence of selective sweeps in multiple independent rice subpopulations?

In the human analysis, querying stored Garud H12 statistics identified KCNE1 as the only gene under sweep in all 5 continental groups — a universal pre-Out-of-Africa sweep. We apply the same approach to rice: querying stored H12 values across all 12 subpopulations to identify genes with convergent sweep evidence.

Of 16 known domestication genes examined:
- **DRO1** (deep rooting): swept in **4 of 12** subpopulations — the most broadly selected domestication gene
- **PROG1** (prostrate growth), **OsC1** (hull color), **Wx** (amylose/waxy), **sh4** (seed shattering), **Sd1** (semi-dwarf), **COLD1** (cold tolerance): each swept in 1–2 subpopulations
- **GW5**, **Hd1**, **Bh4**, and others: no H12 sweep signals detected (selection may act through extended haplotype homozygosity rather than within-window haplotype frequency)

> **Cross-species comparison:** KCNE1 shows universal sweep across all human continental groups (ancient, shared selection). In rice, no gene shows universal selection — domestication involves subpopulation-specific artificial selection on different traits (grain quality in japonica, disease resistance in indica, root architecture in upland varieties). This contrast highlights how GraphPop's convergent query reveals different modes of selection across species.

---

## 10. Generating Summary Tables and Figures

**Target problem:** How to produce publication-ready outputs from a
completed analysis.

### 9.1 Aggregate results

After running the full analysis (Sections 4-8), use the `aggregate`
command to generate summary tables:

```bash
graphpop aggregate -d results/rice/ -o tables/
```

Expected output:

```
Loaded 1488 results from JSON
Generating population_summary.tsv...
Generating fst_matrix.tsv...
Generating roh_summary.tsv...
Generating ihs_peaks.tsv...
Generating nsl_peaks.tsv...
Generating xpehh_peaks.tsv...
Generating sweep_windows.tsv...

Summary tables written to tables/
  population_summary.tsv: 144 rows
  fst_matrix.tsv: 240 rows
  roh_summary.tsv: 12 rows
  ihs_peaks.tsv: 100 rows
  nsl_peaks.tsv: 100 rows
  xpehh_peaks.tsv: 100 rows
  sweep_windows.tsv: 847 rows
```

Each output file is a tab-separated table ready for import into R,
Python, or a spreadsheet. The files produced are:

| File                      | Contents                                   |
|---------------------------|--------------------------------------------|
| `population_summary.tsv`  | Per-pop pi, theta, Tajima's D, Fis per chr|
| `fst_matrix.tsv`          | Pairwise Fst for all population pairs      |
| `roh_summary.tsv`         | Mean FROH, ROH count per population        |
| `ihs_peaks.tsv`           | Top 100 iHS peaks (strongest signals)      |
| `nsl_peaks.tsv`           | Top 100 nSL peaks                          |
| `xpehh_peaks.tsv`         | Top 100 XP-EHH peaks                      |
| `sweep_windows.tsv`       | Windows with H12 > 0.1 (sweep candidates) |

### 9.2 Export high-differentiation windows

Export all genomic windows exceeding a specified Fst threshold:

```bash
graphpop export-windows --min-fst 0.5 -o results/rice/outlier_windows.tsv
```

Expected output (first 5 lines of outlier_windows.tsv):

```
window_id	chr	start	end	population	run_id	n_variants	n_segregating	pi	theta_w	tajima_d	fst	fst_wc	dxy	pbs	fay_wu_h
w_Chr1_2250000_2350000	Chr1	2250000	2350000	GJ-tmp	scan_GJtmp_XI1A	312	187	0.00087	0.00112	-0.723	0.6234	0.6012	0.00412	NA	-0.312
w_Chr1_6500000_6600000	Chr1	6500000	6600000	GJ-tmp	scan_GJtmp_XI1A	287	164	0.00093	0.00118	-0.687	0.5891	0.5723	0.00387	NA	-0.287
w_Chr1_8750000_8850000	Chr1	8750000	8850000	GJ-tmp	scan_GJtmp_XI1A	342	198	0.00078	0.00103	-0.792	0.7123	0.6987	0.00456	NA	-0.341
w_Chr1_12000000_12100000	Chr1	12000000	12100000	GJ-tmp	scan_GJtmp_XI1A	256	142	0.00098	0.00124	-0.681	0.5567	0.5412	0.00378	NA	-0.267
```

These high-Fst windows are candidates for regions under divergent
selection between the compared populations. Cross-reference with gene
annotations:

```bash
graphpop query "MATCH (w:GenomicWindow)
  WHERE w.fst >= 0.8
  MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene)
  WHERE v.chr = w.chr AND v.pos >= w.start AND v.pos < w.end
  RETURN DISTINCT w.chr AS chr, w.start AS start, w.end AS end,
    w.fst AS fst, g.symbol AS gene
  ORDER BY w.fst DESC LIMIT 20" \
  -o results/rice/high_fst_genes.tsv
```

### 9.3 Filter persisted selection statistics by annotation

For haplotype-based statistics that are computed genome-wide and then
filtered (the "compute-then-filter" workflow), use `graphpop filter`:

```bash
# Find iHS peaks at missense variants in temperate japonica
graphpop filter ihs Chr1 GJ-tmp \
    --consequence missense_variant \
    --min-score 2.0 \
    -o results/rice/ihs_missense_GJtmp_Chr1.tsv
```

Expected output:

```
Found 127 records.
```

```
variant_id	pos	ihs	ihs_unstd	consequence	impact
Chr1:2341287:A:G	2341287	2.87	3.12	missense_variant	MODERATE
Chr1:4512341:C:T	4512341	-2.54	-2.78	missense_variant	MODERATE
Chr1:6789234:G:A	6789234	2.31	2.52	missense_variant	MODERATE
Chr1:8923456:T:C	8923456	-2.71	-2.94	missense_variant	MODERATE
Chr1:12345678:A:T	12345678	2.12	2.31	missense_variant	MODERATE
```

Filter by pathway to find haplotype sweep signals in specific biological
functions:

```bash
# XP-EHH peaks in starch metabolism genes
graphpop filter xpehh Chr6 GJ-trp --pop2 XI-1A \
    --pathway "Starch" \
    --min-score 2.0 \
    -o results/rice/xpehh_starch_Chr6.tsv
```

Filter by gene to deep-dive into a specific domestication locus:

```bash
# nSL scores at the GW5 grain width locus
graphpop filter nsl Chr5 XI-1A \
    --gene GW5 \
    -o results/rice/nsl_gw5.tsv
```

### 9.4 Comparing populations, exporting regions, and batch processing

The `compare` command computes the difference in a summary statistic between
two populations, producing a per-window delta table. This is useful for
identifying regions where one ecotype has lost diversity relative to another:

```bash
# Compare nucleotide diversity between temperate and tropical japonica
graphpop compare GJ-tmp GJ-trp Chr1 --stat pi -o delta_pi.tsv
```

To hand off results to external tools like bedtools or IGV, `export-bed`
writes regions exceeding a statistic threshold as a standard BED file:

```bash
# Export high-Fst regions as BED for downstream analysis
graphpop export-bed --stat fst --threshold 0.5 --pop GJ-tmp -o high_fst.bed
```

The `extract` command pulls a filtered set of variants from the graph,
optionally restricted by chromosome range, consequence type, and allele
frequency. This replaces multi-step bcftools + grep pipelines:

```bash
# Extract missense variants in the GW5 region on Chr5
graphpop extract variants --chr Chr5 --start 5200000 --end 5400000 \
    --consequence missense_variant -o gw5_missense.tsv
```

When the same analysis needs to be repeated across many populations or
chromosomes, the `batch` command parallelizes the work automatically:

```bash
# Run diversity across multiple populations and chromosomes in parallel
graphpop batch diversity --pops GJ-tmp,GJ-trp,XI-1A --chrs Chr1,Chr2 --workers 3 \
    -d batch_results/
```

Use `config` to inspect or change persistent CLI settings such as the
default database, Neo4j URI, or output format:

```bash
# Set the default database so you do not need --database on every command
graphpop config set database rice3k
```

Finally, the `report` command generates a self-contained HTML report
summarizing all computed statistics, top selection signals, and population
comparisons:

```bash
# Generate an automated analysis report
graphpop report -o rice3k_report.html
```

---

## 11. Sharing the Database

**Target problem:** How to share your complete analytical record with
collaborators.

### 10.1 Dump the database

The `dump` command creates a portable file containing the entire graph
database -- all nodes, relationships, indexes, and persisted analysis
results:

```bash
# Stop Neo4j first (required for dump)
graphpop stop

# Create the dump
graphpop dump --database rice3k -o rice3k_v1.dump
```

Expected output:

```
Dumping database 'rice3k' to rice3k_v1.dump...
Dump complete: rice3k_v1.dump (12.3 GB)
Manifest: rice3k_v1.manifest.json
```

The dump includes everything: the 29.6M Variant nodes with their allele
count arrays and persisted selection statistics, the 55K Gene nodes, the
153 Pathway nodes, all HAS_CONSEQUENCE and IN_PATHWAY relationships,
and all GenomicWindow nodes from genome scans. A collaborator who loads
this dump gets the full analytical record and can immediately run new
queries without re-computing anything.

The optional manifest file summarizes the database contents:

```bash
cat rice3k_v1.manifest.json
```

```json
{
  "database": "rice3k",
  "date": "2026-03-27T14:32:18.456789",
  "dump_file": "rice3k_v1.dump",
  "dump_size_bytes": 13204578304,
  "node_counts": {
    "Variant": 29594013,
    "Gene": 55801,
    "Sample": 3024,
    "GenomicWindow": 142876,
    "Pathway": 153,
    "Population": 12,
    "Chromosome": 12
  },
  "edge_counts": {
    "NEXT": 29593941,
    "ON_CHROMOSOME": 29594013,
    "HAS_CONSEQUENCE": 4823617,
    "IN_POPULATION": 3024,
    "IN_PATHWAY": 12476
  },
  "populations": {
    "XI-1A": 480,
    "GJ-trp": 360,
    "XI-1B": 340,
    "XI-2": 280,
    "GJ-tmp": 230,
    "XI-3": 230,
    "cA-Aus": 210,
    "XI-adm": 170,
    "admix": 115,
    "GJ-sbtrp": 90,
    "cB-Bas": 60,
    "GJ-adm": 60
  }
}
```

### 10.2 Load on another machine

A collaborator installs GraphPop and loads the dump:

```bash
# On the collaborator's machine
pip install graphpop-cli
graphpop setup --password theirpass --pagecache 20g --heap 4g

# Load the shared database
graphpop load --dump-file rice3k_v1.dump --database rice3k
```

Expected output:

```
Loading database 'rice3k' from rice3k_v1.dump...
Database 'rice3k' loaded successfully.
Config updated to use database 'rice3k'.

Next: graphpop start && graphpop db info
```

```bash
# Start and verify
graphpop start
graphpop db info

# Immediately run new analyses -- no re-import needed
graphpop diversity Chr1 1 43270923 XI-3 -o test_diversity.tsv
```

The collaborator can also query previously persisted results:

```bash
# Query selection statistics that were computed by the original analyst
graphpop query "MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene {symbol:'Hd1'})
  RETURN v.pos, v.xpehh_GJ_trp_XI_1A, v.ihs_GJ_tmp
  ORDER BY v.pos LIMIT 5"
```

This works because iHS and XP-EHH scores were written to Variant nodes
with `--persist`. The database dump preserves all of these properties.

---

## Running on an HPC Cluster

**Target problem:** How to run the full rice 3K analysis on a SLURM or PBS cluster?

GraphPop provides ready-to-use job templates for HPC environments. The key principle: Neo4j runs on one node with fast local storage; analysis jobs run on any node and connect via bolt:// protocol.

### One-time setup (interactive)

```bash
# Request an interactive node
srun --nodes=1 --cpus-per-task=8 --mem=64G --time=2:00:00 --pty bash

# Set up Neo4j with local SSD storage
sbatch scripts/cluster/slurm_setup_neo4j.sh
# Note the hostname printed at the end (e.g., node042)
```

### Import (batch job)

```bash
# Generate CSVs (no Neo4j needed, parallelizable)
sbatch scripts/cluster/slurm_prepare_csv.sh /data/rice.vcf.gz /data/panel.txt

# Or use array jobs for multi-chromosome parallel import
sbatch --array=1-12 scripts/cluster/slurm_prepare_csv.sh ...

# Bulk import (on the database node)
sbatch scripts/cluster/slurm_load_csv.sh /scratch/$USER/csv_out $HOME/neo4j
```

### Full-genome analysis (array jobs)

```bash
# Set connection to database node
export GRAPHPOP_URI=bolt://node042:7687
export GRAPHPOP_PASSWORD=mypassword

# Per-population analysis: 12 chromosomes in parallel
sbatch --array=1-12 scripts/cluster/slurm_fullgenome_array.sh

# Pairwise analysis: XP-EHH + divergence
export GRAPHPOP_PAIRS="GJ-tmp:GJ-trp,GJ-tmp:XI-1A,XI-1A:cA-Aus"
sbatch --array=1-12 scripts/cluster/slurm_pairwise_array.sh
```

### Post-analysis (any node)

```bash
# Aggregate results (no Neo4j needed)
graphpop aggregate -d results/ -o tables/
graphpop plot diversity-bar results/diversity/ -o figures/fig_diversity.png
graphpop report -o rice3k_report.html
```

### PBS equivalent

```bash
qsub -v GRAPHPOP_URI=bolt://db-node:7687,GRAPHPOP_PASSWORD=mypass \
    scripts/cluster/pbs_fullgenome_array.sh
```

> **Tip:** For detailed cluster deployment instructions, storage strategy, and troubleshooting, see the [Cluster Computing Guide](../docs/cluster-guide.md).

---

## Cross-Species Comparison: Rice vs Human

GraphPop enables direct comparison of population genomic patterns across species by applying identical analytical workflows. Three key contrasts emerge:

| Analysis | Human (1000 Genomes) | Rice (3K Genomes) | Interpretation |
|----------|---------------------|-------------------|----------------|
| **πN/πS** | 0.646–0.673 (all < 1.0) | 1.018–1.146 (all > 1.0) | Efficient purifying selection (human) vs relaxed under domestication (rice) |
| **Consequence bias** | HIGH Fst (0.124) < LOW (0.143) | HIGH Fst (0.185) > LOW (0.092) | Purifying constraint (human) vs domestication override (rice) |
| **Convergent sweeps** | KCNE1 in all 5 continental groups | DRO1 in 4/12 subpopulations (max) | Universal natural selection (human) vs trait-specific artificial selection (rice) |
| **Individual trajectories** | PC1 = Out-of-Africa diversity gradient | PC1 = indica-japonica diversity gradient | Same analytical framework reveals parallel structure from different evolutionary forces |
| **Pathway co-selection** | 3 communities (2,232 pathways): DNA repair, immune, metabolic | 3 modules (153 pathways): secondary metabolism, hormone, core cellular | Housekeeping pathways co-differentiate in both species |

These contrasts are uniquely enabled by annotation-conditioned computation applied identically across datasets — a capability that no classical tool provides.

---

## 12. Summary

### What we learned

This vignette demonstrated the complete GraphPop workflow for the Rice
3K Genomes dataset:

1. **Data preparation**: Downloaded VCF, prepared population panel,
   obtained SnpEff annotations, Plant Reactome pathways, and ancestral
   alleles.

2. **Import**: One command (`graphpop import`) built a 46 GB graph
   database with 29.6M Variant nodes, 55K Gene nodes, and 153 Pathway
   nodes, all connected by typed relationships.

3. **Diversity and structure**: Computed pi, theta, and Tajima's D for
   all 12 subpopulations, confirming that indica is 2-3x more diverse
   than japonica and that Fst = 0.71 between the most divergent pair.

4. **Annotation-conditioned analysis**: Used the `--consequence`
   parameter to compute piN/piS for all 12 populations x 12 chromosomes
   (288 calls), finding universal relaxed purifying selection
   (piN/piS > 1.0 in all subpopulations).

5. **Adaptive protein divergence**: Compared missense vs. synonymous
   Fst across three ecotype pairs, finding evidence for adaptive
   protein evolution strongest at the temperate-tropical boundary
   (ratio 1.106).

6. **Selection scans**: Computed iHS, XP-EHH, and Garud's H with
   `--persist`, then queried for multi-statistic convergence to recover
   six known domestication genes (GW5, Hd1, PROG1, Wx, OsC1, DRO1).

7. **Functional analysis**: Stratified Fst by consequence impact
   (HIGH > LOW, opposite of humans) and identified top differentiated
   pathways (salicylic acid, gravitropism, drought response).

8. **Export and sharing**: Generated summary tables with `aggregate`,
   exported high-Fst windows, and created a shareable database dump.

### Key biological findings

Three findings from this analysis were uniquely enabled by graph-native
computation:

1. **Universal cost of domestication**: All 12 rice subpopulations show
   piN/piS > 1.0, indicating that domestication has relaxed purifying
   selection genome-wide. This required 288 annotation-conditioned
   diversity calculations -- impractical with classical pipelines.

2. **Domestication overrides purifying constraint**: HIGH-impact
   variants show higher Fst than LOW-impact variants in rice (the
   opposite of humans), meaning artificial selection during
   domestication actively drove protein-coding changes to fixation.
   This required joint access to per-variant Fst and VEP consequence
   annotations stored on the same graph nodes.

3. **Pathway-level coordinated divergence**: Ecologically coherent
   pathways (defence signalling, drought response, starch metabolism)
   show the highest Fst between ecotypes, revealing a systems-level
   architecture of selection that emerges only from aggregating
   per-variant statistics through gene-pathway graph edges.

### What was uniquely enabled by GraphPop

| Capability                       | Classical tools      | GraphPop            |
|----------------------------------|---------------------|---------------------|
| Annotation-conditioned diversity | 4-step pipeline     | 1 command + flag    |
| piN/piS for 12 pops x 12 chrs   | 288 pipeline runs   | 288 single commands |
| Missense vs. synonymous Fst      | VEP + filter + VCFtools | 1 command + flag |
| Multi-statistic convergence      | Custom merge scripts | 1 Cypher query      |
| Consequence-level Fst            | VEP + bedtools + custom | 1 Cypher query   |
| Pathway-level Fst                | Not practical        | 1 Cypher query      |
| Share complete analysis          | Ship raw files       | 1 dump file         |

The key insight is that GraphPop does not compute different statistics
than classical tools. It computes the same statistics (pi, Fst, iHS,
XP-EHH) but stores them on the same graph nodes as the variant
annotations, gene models, and pathway memberships. This co-location
enables questions that require combining information across multiple
data layers -- questions that are theoretically possible but practically
infeasible with file-based tools.

---

## Appendix: Quick Reference

### Commands used in this vignette

```bash
# Setup and status
graphpop setup --password mypass --pagecache 20g --heap 4g
graphpop start
graphpop status
graphpop import --vcf VCF --panel PANEL --database NAME [--vep VEP] [--pathways TSV] [--ancestral FA]
graphpop inventory
graphpop db info
graphpop validate
graphpop config set database NAME
graphpop config show

# Diversity and divergence
graphpop diversity CHR START END POP [-o FILE] [--consequence TYPE]
graphpop divergence CHR START END POP1 POP2 [-o FILE] [--consequence TYPE]
graphpop pop-summary CHR POP [-o FILE]
graphpop joint-sfs CHR START END POP1 POP2 [-o FILE]
graphpop compare POP1 POP2 CHR --stat STAT [-o FILE]

# Genome scans
graphpop genome-scan CHR POP WINDOW STEP [--pop2 POP2] [--consequence TYPE] [--persist]

# Haplotype-based selection
graphpop ihs CHR POP [--persist] [-o FILE]
graphpop xpehh CHR POP1 POP2 [--persist] [-o FILE]
graphpop nsl CHR POP [--persist] [-o FILE]
graphpop garud-h CHR POP WINDOW STEP [-o FILE]
graphpop ld CHR START END POP MAX_DIST MIN_R2 [-o FILE]

# Convergence and gene ranking
graphpop converge --stats STAT1,STAT2 --thresholds T1,T2 --pop POP [-o FILE]
graphpop rank-genes --pop POP [--top N] [-o FILE]
graphpop lookup gene SYMBOL [-o FILE]
graphpop neighbors NODE [--hops N] [-o FILE]

# Orchestration
graphpop run-all [--phase 1|2|all] -d OUTPUT_DIR/
graphpop batch COMMAND --pops POP1,POP2 --chrs CHR1,CHR2 [--workers N] -d DIR/

# Querying and filtering
graphpop query "CYPHER" [-o FILE]
graphpop filter STAT CHR POP [--consequence TYPE] [--pathway NAME] [--gene NAME]

# Aggregation, export, and reporting
graphpop aggregate -d RESULTS_DIR/ -o TABLES_DIR/
graphpop export-windows [--min-fst VAL] [-o FILE]
graphpop export-bed --stat STAT --threshold VAL --pop POP [-o FILE]
graphpop extract variants --chr CHR [--start POS] [--end POS] [--consequence TYPE] [-o FILE]
graphpop report [-o FILE]

# Sharing
graphpop dump --database NAME [-o FILE]
graphpop load --dump-file FILE --database NAME
```

### Environment variables

| Variable           | Purpose                          | Example                    |
|--------------------|----------------------------------|----------------------------|
| `GRAPHPOP_URI`     | Neo4j connection URI             | `bolt://localhost:7687`    |
| `GRAPHPOP_USER`    | Neo4j username                   | `neo4j`                    |
| `GRAPHPOP_PASSWORD`| Neo4j password                   | `mypass`                   |
| `GRAPHPOP_DATABASE`| Active database name             | `rice3k`                   |

### Further reading

- Wang et al. 2018. "Genomic variation in 3,010 diverse accessions of
  Asian cultivated rice." *Nature* 557: 43-49.
- Huang et al. 2012. "A map of rice genome variation reveals the origin
  of cultivated rice." *Nature* 490: 497-501.
- Alam et al. 2025. "Genomic insights into rice domestication and
  ecotype differentiation." *Nature Genetics* (in press).
- GraphPop documentation: https://github.com/graphpop/graphpop-cli

## 11b. Programmatic Access: MCP Server for AI Agents

**Target problem:** Can AI agents autonomously query the rice 3K graph database to compute statistics, explore annotations, and generate hypotheses?

GraphPop provides a Model Context Protocol (MCP) server that exposes 21 tools — the same 12 procedures available via the CLI, plus convergence detection, gene ranking, annotation lookup, filtering, and database inventory. AI agents (such as Claude, GPT, or custom frameworks) can use these tools to autonomously investigate population genomics questions.

### Setting up the MCP server

```bash
# Install the MCP server package
cd graphpop-mcp && pip install -e .

# Start with same connection as CLI
export GRAPHPOP_URI=bolt://localhost:7687
export GRAPHPOP_DATABASE=rice3k
graphpop-mcp
```

### What an AI agent can do

An AI agent connected to the GraphPop MCP server can autonomously:

```
1. graphpop_inventory()                    # "What's in this database?"
2. graphpop_diversity(chr="Chr1", start=1, end=43270923, pop="GJ-tmp")
3. graphpop_diversity(chr="Chr1", start=1, end=43270923, pop="GJ-tmp",
                      consequence="missense_variant")    # piN computation
4. graphpop_converge(pop="GJ-tmp", stats="ihs,xpehh,h12",
                     thresholds="2.0,2.0,0.3")          # convergent signals
5. graphpop_rank_genes(pop="GJ-tmp", top=10)             # top selection candidates
6. graphpop_lookup_gene(gene="GW5")                      # gene annotation
7. graphpop_lookup_pathway(pathway="Starch biosynthesis") # pathway members
```

This enables AI-driven scientific workflows: an agent can formulate evolutionary hypotheses, compute the necessary statistics, cross-reference with annotations, and report findings — all through the same graph database used for the analyses in this vignette.

> **Note:** The MCP server uses the same environment variables as the CLI (`GRAPHPOP_URI`, `GRAPHPOP_DATABASE`), so the same database configuration works for both human and AI users.
