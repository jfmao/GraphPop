# GraphPop Data Import Strategy

## Overview

GraphPop transforms VCF (Variant Call Format) files into a Neo4j property graph
for population genomics analysis. The import pipeline is a two-stage process:

1. **Python stage** — Parse VCF, compute per-population statistics, emit CSV files
2. **Neo4j stage** — Bulk-load CSVs via `neo4j-admin database import full`

The pipeline is species-agnostic (tested on human 1000 Genomes and rice 3K) and
handles diploid, haploid, and mixed-ploidy chromosomes.

---

## Graph Schema

### Nodes

| Node Label   | Key Properties | Description |
|-------------|----------------|-------------|
| **Variant** | `variantId` (chr:pos:ref:alt), `chr`, `pos`, `ref`, `alt`, `variant_type`, per-pop arrays (`ac`, `an`, `af`, `het_count`, `hom_alt_count`, `het_exp`), global summaries (`ac_total`, `an_total`, `af_total`, `call_rate`), ancestral annotation (`ancestral_allele`, `is_polarized`), **packed genotype arrays** (`gt_packed`, `phase_packed`, `ploidy_packed`) | One node per biallelic variant site |
| **Sample**  | `sampleId`, `population`, `packed_index`, `sex` | One node per individual in the panel |
| **Population** | `populationId`, `name`, `n_samples`, `a_n`, `a_n2` | One node per population group |
| **Chromosome** | `chromosomeId`, `length` | One node per chromosome/contig |
| **Gene** *(optional)* | `geneId`, `symbol` | From VEP/SnpEff annotation |

### Relationships

| Relationship | Pattern | Description |
|-------------|---------|-------------|
| **NEXT** | `(:Variant)-[:NEXT {distance_bp}]->(:Variant)` | Linked list of variants along each chromosome, sorted by position |
| **ON_CHROMOSOME** | `(:Variant)-[:ON_CHROMOSOME]->(:Chromosome)` | Which chromosome a variant belongs to |
| **IN_POPULATION** | `(:Sample)-[:IN_POPULATION]->(:Population)` | Population membership |
| **HAS_CONSEQUENCE** *(optional)* | `(:Variant)-[:HAS_CONSEQUENCE {consequence}]->(:Gene)` | Functional annotation from VEP/SnpEff |

---

## Why Not CARRIES Edges?

The original design used per-genotype relationships:

```
(:Sample)-[:CARRIES {dosage, phase}]->(:Variant)
```

One CARRIES edge was created for every non-reference genotype call. This was
abandoned in Milestone 3.1 for the following reasons:

### 1. Edge explosion

For 1,070,401 variants × 3,202 samples (1000 Genomes chr22), the majority of
genotype calls are non-reference in at least one population. This produces
**billions of CARRIES edges** — roughly 6.3 billion for the full chr22 dataset.

At ~60-80 bytes per relationship in Neo4j's store (relationship record + properties),
this translates to **~200+ GB** of on-disk storage for a single chromosome. A
whole-genome dataset (80-90M variants) would require multiple terabytes just for
genotype edges.

### 2. Traversal cost

Every population genetics computation — allele frequency, heterozygosity,
haplotype reconstruction — requires iterating over all genotypes for a set of
variants. With CARRIES edges, this means:

- For each variant, traverse all incoming CARRIES edges
- Filter by population (join to Sample → Population)
- Extract dosage and phase from edge properties

This is O(V × S) random-access graph traversal, where V = variants and
S = samples. Even with Neo4j's traversal engine, this is orders of magnitude
slower than sequential array reads.

### 3. Haplotype matrix construction

EHH-based statistics (iHS, XP-EHH, nSL) and LD computation require building
a dense haplotype matrix: `haplotypes[variant][sample*2 + phase]`. With CARRIES
edges, constructing this matrix requires one graph traversal per variant,
with random-access pattern across the relationship store. The packed-array
approach reads a contiguous byte sequence per variant — one property access
instead of hundreds of edge traversals.

### 4. Import time

Creating billions of relationships during `neo4j-admin import` is slow — both
in CSV generation (writing the carries_edges.csv file) and in the import itself
(building relationship store + indexes). Removing CARRIES edges cuts import time
roughly in half.

---

## Packed Genotype Arrays: The Replacement

Instead of per-sample edges, genotype data is stored as **compact byte arrays
directly on Variant nodes**:

### gt_packed (2 bits per sample)

Four samples are packed into each byte, LSB-first:

```
Byte layout:  [s0:2bits][s1:2bits][s2:2bits][s3:2bits]
Bit positions: [1:0]     [3:2]     [5:4]     [7:6]

Encoding:
  00 = HomRef
  01 = Het
  10 = HomAlt
  11 = Missing
```

For 3,202 samples: `ceil(3202/4) = 801 bytes` per variant.

The 2-bit encoding is remapped from cyvcf2's convention during import:
- cyvcf2: 0=HomRef, 1=Het, **2=Missing, 3=HomAlt**
- packed: 0=HomRef, 1=Het, **2=HomAlt, 3=Missing**

This remap makes the ALT allele count computable from the 2-bit value directly
(`gt_value` for HomRef=0, Het=1, HomAlt=2, with Missing=3 as sentinel).

### phase_packed (1 bit per sample)

Eight samples packed per byte, LSB-first. The bit indicates which haplotype
carries the ALT allele for heterozygous sites (0 or 1). Only meaningful when
`gt_packed` indicates Het.

For 3,202 samples: `ceil(3202/8) = 401 bytes` per variant.

### ploidy_packed (1 bit per sample)

Same layout as phase_packed. Bit = 1 means the sample is haploid at this
variant (e.g., males on chrX). **Null or empty** means all samples are diploid
(the common case for autosomes), adding zero storage overhead.

### Storage comparison

| Approach | Per-variant storage (3,202 samples) | chr22 (1.07M variants) | Whole genome (~85M variants) |
|----------|-------------------------------------|------------------------|------------------------------|
| CARRIES edges | ~150 KB (avg ~2,400 edges × 64 B) | ~160 GB | ~12 TB |
| Packed arrays | **1,202 bytes** (801 + 401) | **~1.3 GB** | **~100 GB** |
| Reduction | **125×** | **125×** | **125×** |

---

## Import Pipeline Workflow

### Step 1: Python — VCF Parsing and CSV Generation

```
VCF file + Panel file
        │
        ▼
   VCFParser (cyvcf2, streaming)
        │
        ├── Intersects VCF samples with panel → PopulationMap
        ├── Detects ploidy per chromosome (auto/diploid mode)
        ├── Reads ##contig headers for chromosome lengths
        ├── Optionally loads ancestral allele FASTA
        │
        ▼
   Per-variant streaming loop
        │
        ├── Filter: PASS-only, biallelic-only, contig whitelist
        ├── Compute per-population: AC, AN, AF, het_count, hom_alt_count, het_exp
        ├── Pack genotypes: gt_packed (2 bit/sample, vectorized numpy)
        ├── Pack phase: phase_packed (1 bit/sample, het sites only)
        ├── Pack ploidy: ploidy_packed (1 bit/sample, non-diploid chr only)
        ├── Annotate ancestral allele (if FASTA provided)
        │
        ▼
   CSVEmitter (streaming, 100K-variant chunks)
        │
        ├── variant_nodes.csv     — Variant nodes with all properties + packed arrays
        ├── sample_nodes.csv      — Sample nodes with packed_index + sex
        ├── population_nodes.csv  — Population nodes with harmonic numbers
        ├── chromosome_nodes.csv  — Chromosome nodes with lengths
        ├── next_edges.csv        — NEXT linked list (built during streaming, sorted input assumed)
        ├── on_chromosome_edges.csv — ON_CHROMOSOME relationships
        └── in_population_edges.csv — IN_POPULATION relationships
```

**7 CSV files total.** No carries_edges.csv.

The packed byte arrays are serialized in CSV as semicolon-delimited signed
Java bytes (range -128 to 127), matching Neo4j's `byte[]` import type.

### Step 2: Optional — VEP/SnpEff Annotation

If functional annotation is available, a separate step produces:
- `gene_nodes.csv` — Gene nodes (geneId, symbol)
- `has_consequence_edges.csv` — HAS_CONSEQUENCE relationships (variant → gene with consequence type)

These are loaded alongside the core 7 CSVs.

### Step 3: neo4j-admin Bulk Import

```bash
neo4j-admin database import full \
    --array-delimiter=";" \
    --nodes=Variant=variant_nodes.csv \
    --nodes=Sample=sample_nodes.csv \
    --nodes=Population=population_nodes.csv \
    --nodes=Chromosome=chromosome_nodes.csv \
    --nodes=Gene=gene_nodes.csv \          # optional
    --relationships=NEXT=next_edges.csv \
    --relationships=ON_CHROMOSOME=on_chromosome_edges.csv \
    --relationships=IN_POPULATION=in_population_edges.csv \
    --relationships=HAS_CONSEQUENCE=has_consequence_edges.csv \  # optional
    neo4j
```

### Step 4: Post-Import Indexing

After starting Neo4j, create indexes:
- Uniqueness constraints on all node ID properties
- Composite index on `(Variant.chr, Variant.pos)` for range queries
- Index on `Variant.variant_type` for SNP/INDEL filtering
- Index on `Gene.symbol` for gene-level queries

---

## How Procedures Access Packed Genotypes

### FAST PATH (allele-count statistics)

Procedures like `graphpop.diversity`, `graphpop.divergence`, `graphpop.sfs`,
`graphpop.joint_sfs`, and `graphpop.pop_summary` primarily use the **pre-computed
per-population arrays** stored on Variant nodes (`ac[]`, `an[]`, `af[]`, etc.).
This requires zero genotype unpacking — just array indexing by population index.

When a user provides a **custom sample subset** via the `samples` option,
`SampleSubsetComputer` reads the packed genotype arrays and recomputes AC/AN/AF
on the fly for only the specified samples, using their `packed_index` values to
address the correct bits.

### FULL PATH (haplotype-level statistics)

Procedures like `graphpop.ld`, `graphpop.ihs`, `graphpop.xpehh`, `graphpop.nsl`,
`graphpop.roh`, and `graphpop.garud_h` need per-sample haplotype data. These use
`HaplotypeMatrix`, which:

1. Queries all Variant nodes in the genomic region
2. For each variant, reads `gt_packed` and `phase_packed` from the node
3. Unpacks into a dense `byte[variant][2*sample+phase]` matrix in memory
4. All subsequent computation runs on this in-memory matrix with zero further DB access

### PackedGenotypeReader

A shared utility class provides bit-extraction operations:
- `genotype(gtPacked, sampleIdx)` → 2-bit genotype (0/1/2/3)
- `phase(phasePacked, sampleIdx)` → 1-bit phase (0/1)
- `ploidy(ploidyPacked, sampleIdx)` → 1-bit ploidy (0=diploid, 1=haploid)

All operations are branchless bit shifts — `O(1)` per sample, no conditionals.

---

## Design Decisions

### Pre-computed population arrays

Per-population AC/AN/AF arrays are computed during import and stored on Variant
nodes. This means FAST PATH procedures do zero genotype unpacking for standard
population queries — they read the array, index by population position, and
compute. This is why diversity/Fst/SFS queries return in ~1 second regardless of
sample count.

### NEXT linked list

Variants on each chromosome are connected by NEXT edges in position order. This
enables:
- Efficient sliding-window traversal (walk the chain, no sort needed)
- Distance-aware computations (EHH decay uses physical distance from NEXT edges)
- Sequential scan without index lookups

The NEXT chain is built during streaming import — the parser assumes position-sorted
VCF input (standard for bcftools/GATK output).

### Semicolon array delimiter

Neo4j's `neo4j-admin import` uses a configurable array delimiter. We use `;`
because commas appear in CSV field separators and colons appear in variant IDs.

### Variant ID format

`chr:pos:ref:alt` (e.g., `chr22:16050075:A:G`). Globally unique, human-readable,
and encodes location directly in the ID.

### Ploidy handling

Ploidy is detected once per chromosome from the first variant's genotype field
width (cyvcf2). Three modes:
- **all_diploid**: no ploidy_packed emitted (null = diploid, zero overhead)
- **all_haploid**: all samples haploid (mitochondria, chloroplast)
- **mixed**: per-sample ploidy recorded (chrX males haploid, females diploid)

This avoids storing ploidy metadata for autosomes where it adds no information.

---

## Pipeline Usage

```bash
# Generate CSVs from VCF
conda run -n graphevo python -m graphpop_import.bulk_import \
    data/raw/1000g/chr22.vcf.gz \
    data/raw/1000g/integrated_call_samples_v3.20130502.ALL.panel \
    data/raw/1000g/csv_out \
    --stratify-by superpopulation

# Import into Neo4j
sudo bash scripts/neo4j-import.sh data/raw/1000g/csv_out neo4j

# Optional: load VEP annotations
bash scripts/load-annotations.sh data/raw/1000g/csv_out

# Start Neo4j and apply indexes
sudo neo4j start
cat data/raw/1000g/csv_out/post_import_indexes.cypher | cypher-shell -d neo4j
```
