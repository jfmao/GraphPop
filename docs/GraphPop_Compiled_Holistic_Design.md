# GraphPop: Compiled Holistic Design for a Graph-Native Population Genomics Platform

## Integrative Evolutionary Inference from Population Variation Graphs — Version 3.0 (Unified)

*Compiled from Design v1 (original architecture) and Revised Design v2 (All of Us scale redesign). All naming unified to GraphPop.*

---

## 1. Project Identity and Thesis

**Full Title:** GraphPop — Integrative Evolutionary Inference from Population Variation Graphs

**Core Thesis:** Every computational task in population genomics — from computing nucleotide diversity in a 1 Mb window, to inferring demographic history, to detecting polygenic selection, to conditioning forward simulations on empirical data — is a traversal, aggregation, diffusion, or augmentation of a single relational variation graph. By making this graph the persistent, queryable, evolvable substrate for all analysis, we eliminate the fragmented VCF-centric workflow that currently dominates the field, and we open the door to relational questions that flat representations cannot express.

**Two Foundational Requirements (from Revised Design v2):**

1. **Target All of Us scale** — 1+ billion variants, 1+ million samples, 100+ TB raw genotype data
2. **All computation must be graph-native** — no external columnar store; the graph database IS the computational substrate

**What GraphPop IS:**

- A graph-native compute engine with a query interface for population genomic inference
- A population variation graph where variants are nodes and their relationships encode genomic, statistical, and functional structure
- A scalable platform from single-chromosome prototyping to All of Us scale (1B+ variants, 1M+ samples)

**What GraphPop is NOT:**

- It is not a wrapper around existing tools — it replaces their core computations with graph-native implementations
- It is not a visualization platform — it is a compute engine with a query interface
- It is not a pan-genome graph (which represents reference sequences as nodes) — it is a population variation graph where variants are nodes

---

## 2. What Has Been Established

### A. Graph Schema Evolution: v2 → v3 (Stratified Aggregation Model)

The original v2 schema stored dosage vectors on Variant nodes. The revised v3 schema (from Revised Design v2) inverts the representation: Variant nodes carry pre-aggregated population-level statistics, while individual genotypes are stored as sparse CARRIES relationships. This redesign was forced by the impossibility of storing million-element dosage vectors at All of Us scale (1B variants × 1M samples = 1 petabyte).

### B. Compute Engine: Dual-Path Architecture

A single Java plugin (`graphpop-procedures.jar`) containing SIMD-accelerated VectorOps primitives (Java Vector API, JDK 21+) and domain-specific procedure classes. Cypher serves as the orchestration and integration layer, not the compute language. The engine now operates in two modes:

- **Fast path:** Population-level statistics from allele count arrays on Variant nodes — complexity O(V × K) where K = number of populations, independent of sample size N
- **Full path:** Individual-level statistics via CARRIES relationship traversal — sparse, efficient for rare variants (>90% of all variants at All of Us scale)

### C. Development Sequence (Phases 1–7)

1. Summary statistics (SFS-based, diversity, divergence, LD, haplotype)
2. Pairwise sample statistics (kinship, IBS/IBD distances)
3. Structure and dimensionality reduction (PCA, community detection, graph embeddings)
4. Selection inference beyond site statistics (composite likelihood, ML classifiers, embedding anomaly detection)
5. Demographic inference (SFS-based: ∂a∂i/moments integration; ABC)
6. Simulation integration (msprime/SLiM conditioning and feedback)
7. AI-augmented inference (GraphRAG + agentic loops)

### D. Hardware Platform

**Development (Tier 1):** Intel i9-13900K, 64 GB RAM, RTX 4090 (24 GB VRAM), 1 TB + planned 4 TB NVMe, running WSL2 under Windows with Claude Desktop for development assistance.

**Production (Tiers 2–3):** Three-tier deployment from laptop to Neo4j Infinigraph cluster (see Section 8).

### E. Benchmark Datasets

1. **1000 Genomes chr22** — phased, fully annotated, ancestral alleles from EPO alignment (human, small-scale validation)
2. **Rice 3K chr1** — unphased but highly homozygous, partial ancestral alleles from O. rufipogon (crop, medium-scale)
3. **msprime-simulated data** — exact ground truth for all statistics (validation)
4. **All of Us chr22 subset** (if accessible; or simulated equivalent with 50,000 samples) — large-scale scaling validation

### F. Key Architectural Decisions

| Aspect | Original Design (v1) | Revised Design (v2/v3) |
|--------|---------------------|----------------------|
| Genotype storage | Dosage vectors on Variant nodes | Sparse CARRIES relationships (Sample→Variant) |
| Allele frequency data | Computed from dosage vectors | Pre-aggregated population arrays on Variant nodes |
| Summary statistics | Read dosage vectors, compute | Fast path: read allele count arrays (O(V×K)) |
| LD computation | Read dosage vectors for pair | Traverse CARRIES for pair (sparse join) |
| Max sample size | ~10,000 (dosage vector size limit) | 1,000,000+ (sparse relationships scale) |
| Max variant count | ~50M (vector storage limit) | 1.5B+ (node count scales with Infinigraph) |
| External store | Required above 10k samples | None — all in Neo4j |
| Adding new samples | Extend all dosage vectors (touch all nodes) | Add CARRIES edges + update allele counts |
| Per-population vectors | Pre-sliced during import | Population arrays on Variant nodes (K×40 bytes) |
| Sparse LD edges | Above threshold, graph-native | Unchanged |
| Analysis versioning | run_id on all materialized outputs | Unchanged |
| Deployment | Single machine | Tier 1: laptop, Tier 2: small cluster, Tier 3: Infinigraph cluster |
| Neo4j edition | Community (development) | Community (Tier 1), Enterprise + Infinigraph (Tier 2–3) |

---

## 3. The Scale Problem — Honest Arithmetic

### All of Us Dataset (Current and Projected)

From the Nature 2024 paper and the February 2025 data release:

- 245,388 whole genomes (2024 paper), 414,830 srWGS samples (Feb 2025 release)
- 535,000 short-read genomes planned for spring 2026 release
- Target: 1,000,000+ participants with WGS
- More than 1 billion genetic variants identified from 245k genomes
- With 1 million genomes, variant count will likely reach 1.5–2 billion (rare variants scale roughly linearly with sample size)

### Why the Original Dosage-Vector Design Breaks

- 1B Variant nodes × 1M samples × 1 byte/sample = **1 petabyte** — impossible
- Even 500k samples: 1B × 500k × 1 byte = **500 TB** — far exceeds any node-property model

### The Stratified Aggregation Model Solves This

Instead of storing the full genotype matrix, store:

- **Level 1 — Population summaries on Variant nodes:** ~200 bytes base + K × 40 bytes per population. With 50 strata: ~2,200 bytes/variant. For 1.5B variants: **~3.3 TB** — feasible with Infinigraph
- **Level 2 — Sparse CARRIES relationships:** Only non-reference genotypes. ~4.5M CARRIES per sample. For 1M samples: **~4.5 trillion edges, ~135 TB** — at edge of Infinigraph capability
- **Level 3 — Sparse pairwise sample statistics:** Kinship above threshold, IBD segments, k-NN IBS edges — **~160 GB** total

---

## 4. The Holistic Architecture — Seven Layers

```
┌─────────────────────────────────────────────────────────────────────────┐
│  LAYER 7: AI-Augmented Inference                                        │
│  GraphRAG retrieval · LangGraph/AutoGen agents · hypothesis generation  │
│  simulation design · iterative refinement · natural language interface   │
├─────────────────────────────────────────────────────────────────────────┤
│  LAYER 6: Simulation Integration                                        │
│  msprime/SLiM export · simulated data re-import · ABC framework         │
│  parameter estimation · model comparison · simulation-conditioned stats  │
├─────────────────────────────────────────────────────────────────────────┤
│  LAYER 5: Inference Engines                                             │
│  Demographic inference (SFS fitting, ∂a∂i/moments bridge)               │
│  Selection inference (CLR, ML classifiers, embedding anomaly detection) │
│  Recombination inference (LD-decay fitting, hotspot detection)           │
├─────────────────────────────────────────────────────────────────────────┤
│  LAYER 4: Graph Algorithms & Embeddings                                 │
│  Community detection (Leiden/Louvain) · PCA via GRM                     │
│  node2vec / GraphSAGE / GAT · spectral embeddings                       │
│  Label propagation · diffusion-based ancestry · manifold analysis       │
│  Bipartite (Sample–Variant) community detection [v3 new capability]     │
├─────────────────────────────────────────────────────────────────────────┤
│  LAYER 3: Pairwise & Relational Statistics                              │
│  Kinship (KING-robust, GRM) · IBS/IBD distances · pairwise F_ST        │
│  Haplotype sharing · ancestry segment detection · ROH                   │
│  LSH/MinHash for sparse pair finding at scale [v3 optimization]         │
├─────────────────────────────────────────────────────────────────────────┤
│  LAYER 2: Summary Statistics Engine                                     │
│  SFS · π · θ_W · Tajima's D · Fay & Wu's H · Fu & Li's D              │
│  H_o · H_e · F_IS · F_ST · D_xy · iHS · XP-EHH · r² · D'             │
│  genome_scan() with windowed materialization                            │
│  FAST PATH (allele count arrays) + FULL PATH (CARRIES traversal)        │
├─────────────────────────────────────────────────────────────────────────┤
│  LAYER 1: Data Foundation                                               │
│  Graph schema (Stratified Aggregation Model) · VCF import/export        │
│  Functional annotation lift · VectorOps SIMD core (dense + sparse)      │
│  CARRIES relationships · Index strategy · Storage management            │
└─────────────────────────────────────────────────────────────────────────┘
```

**Key Principle:** Each layer only calls downward. Layer 7 (AI agents) orchestrates calls to any layer.

---

## 5. Layer-by-Layer Design

### 5.1 LAYER 1: Data Foundation

**Status:** Mostly designed. Import pipeline specified. Needs implementation.

#### 5.1.1 Graph Schema (v3 — Stratified Aggregation Model)

**Node Types:**

```
(:Variant {
    // Identity
    id: STRING,                  // "chr1:12345:A:T"
    chr: STRING,
    pos: LONG,                   // LONG not INT for large genomes
    ref: STRING,
    alt: STRING,
    variant_type: STRING,        // "SNP", "INDEL", "SV"

    // Ancestral state
    ancestral: STRING,
    is_polarized: BOOLEAN,

    // Population allele count arrays (compact: K populations)
    pop_ids: [STRING],           // ["AFR", "AMR", "EAS", "EUR", "MID", "SAS"]
    ac: [INT],                   // alt allele count per pop
    an: [INT],                   // total allele count per pop (2 × n_samples minus missing)
    af: [FLOAT],                 // allele frequency per pop
    het_count: [INT],            // heterozygote count per pop
    hom_alt_count: [INT],        // homozygous alt count per pop

    // Global summaries
    ac_total: INT,
    an_total: INT,
    af_total: FLOAT,

    // Quality
    call_rate: FLOAT,

    // Pre-computed per-pop statistics (materialized on demand)
    het_exp: [FLOAT],            // expected heterozygosity per pop = 2·p·(1-p)

    // Functional annotation summary (from VEP)
    consequence: STRING,         // most severe consequence
    impact: STRING,              // HIGH, MODERATE, LOW, MODIFIER
    gene_symbol: STRING,
    cadd_score: FLOAT,
})

(:Sample {
    id: STRING,
    population: STRING,
    sex: STRING,
    ancestry_proportions: [FLOAT],
    pca_coords: [FLOAT],
    phenotypes: MAP,             // {trait1: value1, ...}
    embedding: [FLOAT],          // computed in Layer 4
})

(:Population {
    id: STRING,
    name: STRING,
    n_samples: INT,
    a_n: FLOAT,                  // harmonic number for Tajima's D
    a_n2: FLOAT,                 // sum(1/i²)
    b_n: FLOAT,                  // Tajima's D constant
})

(:Chromosome {
    id: STRING,
    length: LONG,
    n_variants: LONG,
    centromere_start: LONG,
    centromere_end: LONG,
})

(:GenomicWindow {
    id: STRING,                  // "chr1:1000000-2000000:EUR:run001"
    chr: STRING,
    start: LONG,
    end: LONG,
    population: STRING,
    run_id: STRING,
    n_variants: INT,
    n_segregating: INT,
    pi: FLOAT,
    theta_w: FLOAT,
    tajima_d: FLOAT,
    fay_wu_h: FLOAT,
    het_obs: FLOAT,
    het_exp: FLOAT,
    fis: FLOAT,
    fst: FLOAT,
    dxy: FLOAT,
    da: FLOAT,
})

(:DenseWindow {
    // Materialized on demand for common-variant analyses (LD, PCA)
    chr: STRING,
    start: LONG,
    end: LONG,
    population: STRING,
    variant_ids: [STRING],
    sample_ids: [STRING],
    genotype_matrix: [INT],      // flattened row-major, common variants only (AF > 0.05)
})

(:Gene {
    id: STRING,
    symbol: STRING,
    chr: STRING,
    start: LONG,
    end: LONG,
    strand: STRING,
    biotype: STRING,
})

(:RegulatoryElement {
    id: STRING,
    type: STRING,                // "promoter", "enhancer", "TFBS"
    chr: STRING,
    start: LONG,
    end: LONG,
})

(:Pathway {
    id: STRING,
    name: STRING,
    source: STRING,              // "KEGG", "Reactome", "GO"
})

(:GOTerm {
    id: STRING,
    name: STRING,
    namespace: STRING,           // "biological_process", "molecular_function", "cellular_component"
})

(:InferenceResult {
    id: STRING,
    method: STRING,
    run_id: STRING,
    parameters: MAP,
    timestamp: DATETIME,
})
```

**Relationship Types:**

```
// === Genomic topology ===
(:Variant)-[:NEXT {distance_bp: LONG, distance_cM: FLOAT}]->(:Variant)
    // Positional chain per chromosome. ~36 GB for 1.5B variants.

(:Variant)-[:ON_CHROMOSOME]->(:Chromosome)

// === Individual genotypes (THE core relationship — sparse) ===
(:Sample)-[:CARRIES {gt: INT, phase: INT}]->(:Variant)
    // gt: 1=HET, 2=HOM_ALT. Homozygous REF implicit (no edge).
    // ~4.5 trillion edges for 1M samples. ~135 TB. Infinigraph cluster.

// === Linkage disequilibrium (sparse, above threshold) ===
(:Variant)-[:LD {r2: FLOAT, dprime: FLOAT, population: STRING}]->(:Variant)

// === Functional annotation ===
(:Variant)-[:HAS_CONSEQUENCE {type: STRING, impact: STRING}]->(:Gene)
(:Gene)-[:IN_PATHWAY]->(:Pathway)
(:Gene)-[:HAS_GO_TERM]->(:GOTerm)
(:GOTerm)-[:IS_A]->(:GOTerm)    // GO ontology hierarchy

// === Population membership ===
(:Sample)-[:IN_POPULATION]->(:Population)

// === Pairwise sample statistics (sparse, above threshold) ===
(:Sample)-[:KINSHIP {coefficient: FLOAT, method: STRING}]->(:Sample)
(:Sample)-[:IBD_SEGMENT {chr: STRING, start: LONG, end: LONG, length_cM: FLOAT}]->(:Sample)

// === Inference provenance ===
(:InferenceResult)-[:COMPUTED_FROM]->(:Population)
(:InferenceResult)-[:CONDITIONED_ON]->(:GenomicWindow)
(:GenomicWindow)-[:ON_CHROMOSOME]->(:Chromosome)
```

#### 5.1.2 VCF Import Pipeline (revised for sparse representation)

Single-pass Python parser (cyvcf2) that:

1. Reads filtered/annotated VCF
2. Emits CSV files for bulk neo4j-admin import (nodes + relationships)
3. For each variant: emit Variant node CSV with population allele counts (fast-path data)
4. For each non-reference genotype: emit CARRIES relationship CSV (Sample → Variant)
5. Skip homozygous reference genotypes (they are implicit)
6. Lifts VEP/SnpEff CSQ into Gene/Consequence nodes
7. Computes NEXT edges (sorted positional chain)
8. Creates Population and Chromosome nodes from metadata

**Chunking strategy for large VCFs:**

- Process one chromosome at a time
- Within chromosome, process in blocks of 100k variants
- Write CARRIES CSVs in parallel per chromosome → merge for bulk import
- Use neo4j-admin import with multiple CSV files in parallel

#### 5.1.3 Export Pipeline (bidirectional)

Stored procedures for exporting graph state to external tool formats:

- `graphpop.export.recomb_map(chr, pop, format, path)` — for msprime/SLiM
- `graphpop.export.sfs(chr, pop, path, format)` — for ∂a∂i/moments/fastsimcoal2
- `graphpop.export.vcf(chr, start, end, samples, path)` — for any classical tool
- `graphpop.export.slim_input(chr, start, end, params, path)` — SLiM initialization

#### 5.1.4 VectorOps Core (SIMD-accelerated, dual-mode)

**Dense-mode primitives (array-based):** `alleleCounts`, `pearsonR2`, `dPrime`, `haplotypeHomozygosity`, `conditionalSubset`, `ibsDistance`, `hammingDistance`, `dotProduct`, `cosineSimilarity`, `outerProduct`, `euclideanDistance`

**Sparse-mode primitives (carrier-list-based — v3 addition):** `pearsonR2_sparse`, `haplotypeHomozygosity_sparse`, `ibsDistance_sparse`

All VectorOps methods have BOTH dense and sparse variants.

#### 5.1.5 Index Strategy

- Composite index on `(Variant.chr, Variant.pos)` — fast range queries for windowed operations
- Unique constraint on `Variant.id`
- Index on `Sample.id`, `Population.id`, `Gene.symbol`
- Full-text index on `GOTerm.name` for functional enrichment queries


### 5.2 LAYER 2: Summary Statistics Engine

**Status:** Designed in detail. Ready for implementation.

#### The Dual-Path Computation Engine

**Fast Path (population-level, no CARRIES traversal):**

| Statistic | Input from Variant node | Formula |
|-----------|------------------------|---------|
| π (nucleotide diversity) | af[k] for population k | Σ 2·p·(1-p) / L |
| θ_W (Watterson's theta) | count where ac[k] > 0 and ac[k] < an[k] | S / a_n |
| Tajima's D | π and θ_W | standard formula |
| Fay & Wu's H | af[k], ancestral state | Σ 2·p_i·p_i / (n(n-1)) |
| H_e (expected het.) | af[k] | 2·p·(1-p) |
| H_o (observed het.) | het_count[k] / (an[k]/2) | direct |
| F_IS | H_o and H_e | 1 - H_o/H_e |
| F_ST (Hudson) | ac[k], an[k] for two pops | standard from allele counts |
| SFS | ac[k], an[k] | histogram of derived allele counts |
| Joint SFS | ac[k1], ac[k2] | 2D histogram |
| D_xy | af[k1], af[k2] | Σ (p1·(1-p2) + p2·(1-p1)) / L |

**Full Path (individual-level, CARRIES traversal):**

| Statistic | Why individual data needed | CARRIES traversal pattern |
|-----------|---------------------------|--------------------------|
| r² / D' (LD) | Pairwise haplotype correlation | Load CARRIES for two variants, join by sample |
| iHS, XP-EHH | Extended haplotype homozygosity | Walk NEXT chain, load CARRIES at each position |
| Kinship (KING) | Pairwise IBS counts | Load CARRIES for all variants in window |
| PCA / GRM | Covariance matrix | Load CARRIES, accumulate outer products |
| IBD detection | Shared haplotype segments | Walk NEXT chain, compare CARRIES between pairs |
| ROH detection | Homozygous stretches | Walk NEXT chain for one sample's CARRIES |

**Procedures:**

| Procedure | Compute Path | Output |
|-----------|-------------|--------|
| `graphpop.diversity(chr, start, end, pop)` | FAST | π, θ_W, Tajima's D, Fay & Wu's H, Fu & Li's D, H_o, H_e, F_IS |
| `graphpop.sfs(chr, start, end, pop)` | FAST | Folded/unfolded SFS vector |
| `graphpop.joint_sfs(chr, start, end, pop1, pop2)` | FAST | 2D SFS matrix |
| `graphpop.divergence(chr, start, end, pop1, pop2)` | FAST | F_ST (Hudson/WC), D_xy, D_a |
| `graphpop.ld(chr, start, end, pop, max_dist, r2_threshold)` | FULL | Sparse LD edges written to graph |
| `graphpop.ihs(chr, pop)` | FULL | iHS scores on Variant nodes |
| `graphpop.xpehh(chr, pop1, pop2)` | FULL | XP-EHH scores |
| `graphpop.genome_scan(chr, pop, window, step)` | FAST (+ optional FULL) | GenomicWindow nodes with all stats |

**Conditioned queries — the graph's unique strength:**

```cypher
// Tajima's D only for missense variants in drought-response genes
CALL graphpop.diversity('chr1', 1000000, 2000000, 'Indica',
     {consequence: 'missense_variant', pathway: 'drought_response'})
```

One Cypher call replaces: extract positions from annotation → filter VCF → re-index → run vcftools → parse output.

**Critical insight at scale:** The most frequently used statistics (π, θ_W, Tajima's D, F_ST, SFS, H_e, H_o, F_IS, Fay & Wu's H) are ALL fast-path. With 1M samples, reading a 50-element allele count array is ~20,000× smaller I/O than a 1M-element dosage vector.


### 5.3 LAYER 3: Pairwise & Relational Statistics

**Status:** Designed conceptually. Needs detailed procedure specification.

**5.3.1 Kinship and Relatedness**

- Procedure: `graphpop.kinship(chr, pop, method)` — methods: KING-robust, GRM, Yang's
- **At scale:** LSH/MinHash on CARRIES edge sets to identify candidate related pairs → exact kinship only for candidates
- Output: `(:Sample)-[:KINSHIP {coefficient, method, run_id}]->(:Sample)` edges

**5.3.2 IBD Segment Detection and Integration**

- Import from external callers (RaPID, hap-IBD, IBDNe) → `(:Sample)-[:IBD_SEGMENT]->(:Sample)` edges
- Graph-native detection (future): walk NEXT chain, detect shared CARRIES patterns between sample pairs

**5.3.3 Pairwise F_ST Matrix**

- Procedure: `graphpop.pairwise_fst(chr, start, end, populations[])` → `(:Population)-[:FST]->(:Population)` edges

**5.3.4 Runs of Homozygosity (ROH)**

- Procedure: `graphpop.roh(sample, chr, min_length, max_het_rate)` → `(:ROHSegment)-[:IN_SAMPLE]->(:Sample)` nodes

**5.3.5 IBS Distance Matrix**

- Procedure: `graphpop.ibs_matrix(chr, start, end, samples[])` → matrix or Sample–Sample edges


### 5.4 LAYER 4: Graph Algorithms & Embeddings

**Status:** Conceptually outlined. Leverages Neo4j GDS library.

**5.4.1 Population Structure via Community Detection**

- Procedure: `graphpop.structure.community(pop, method, resolution)` — Leiden/Louvain via GDS
- Advantage over ADMIXTURE: incremental, no fixed K, handles continuous structure

**5.4.2 PCA via Graph-Native GRM**

- Procedure: `graphpop.structure.pca(pop, n_components, region)` — eigendecomposition via LAPACK/Apache Commons Math

**5.4.3 Graph Embeddings**

| Method | GDS Built-in? | Best For |
|--------|---------------|----------|
| node2vec | Yes | Haplotype block structure, local LD patterns |
| GraphSAGE | Yes (via GDS ML) | Supervised structure prediction |
| Spectral embedding | Partial | Classical MDS-like, interpretable axes |
| Personalized PageRank | Yes | Ancestry proportion estimation |
| FastRP | Yes | Scalable unsupervised embedding |

- Procedure: `graphpop.embed(graph_projection, method, dimensions, params)`
- Embedding targets: Variant-Variant (NEXT + LD), Sample-Sample (KINSHIP + IBS + IBD), Heterogeneous (R-GCN/HGT)

**5.4.4 Label Propagation for Semi-Supervised Classification**

- Procedure: `graphpop.structure.label_propagation(seed_labels, iterations)` — powerful for landrace/wild-relative collections

**5.4.5 Manifold Analysis (Novel)**

- Procedure: `graphpop.manifold.analyze(embedding_property, method)` — local dimensionality, curvature, density ridge detection
- Outputs per-node manifold statistics → new summary statistics for selection scan

**5.4.6 Bipartite Community Detection (v3 new capability)**

- Project the CARRIES bipartite graph → jointly cluster samples and variants
- Degree distribution of Variant nodes = allele frequency spectrum (SFS)
- Degree distribution of Sample nodes = individual heterozygosity
- Shared neighborhood of two Sample nodes = IBS/kinship
- Community structure of bipartite projection = population structure
- This bipartite graph IS the object that population genetics studies


### 5.5 LAYER 5: Inference Engines

**Status:** Conceptual. Consumes outputs from Layers 2–4.

**5.5.1 Demographic Inference**

- **Bridge mode:** `graphpop.demography.bridge(tool, model, params, sfs_source)` — export SFS → ∂a∂i/moments/fastsimcoal2 → re-import InferenceResult
- **Native mode:** lightweight 1-pop/2-pop model fitting in Java (for ABC inner loops)

**5.5.2 Selection Inference — Multi-Signal Integration**

- Procedure: `graphpop.selection.multi_signal_scan(chr, pop, signals[], method)`
- Methods: (a) composite likelihood, (b) random forest on labeled windows, (c) embedding anomaly detection (novel)
- Procedure: `graphpop.selection.functional_enrichment(selection_result, annotation_type)` — single query replaces clusterProfiler workflow

**5.5.3 Recombination Landscape Inference**

- Procedure: `graphpop.recombination.rate_map(chr, pop, method)` — LD-decay fitting or LDhat/pyrho import
- Recombination topology graph via community detection on LD graph

**5.5.4 Admixture & Ancestry Inference**

- Procedure: `graphpop.admixture.local_ancestry(chr, sample, ref_pops)` — graph diffusion from reference seeds
- `(:AncestryTract)-[:IN_SAMPLE]->(:Sample)` nodes — queryable ancestry-switch enrichment


### 5.6 LAYER 6: Simulation Integration

**Status:** Conceptual. Critical for hypothesis testing.

**5.6.1 Simulation Conditioning — Graph → Simulator**

```
Graph property              → Simulation input
──────────────────────────────────────────────────
SFS (from Layer 2)          → ∂a∂i/moments fit → Ne trajectory for msprime
Recomb rate on NEXT edges   → Recombination map for msprime/SLiM
Allele frequencies per pop  → Initial allele frequencies for SLiM
Selection scores (Layer 5)  → Selection coefficients for SLiM
Population structure (L4)   → Migration matrix for msprime
```

- Procedure: `graphpop.simulate.condition(chr, start, end, pop, simulator, output_path)`

**5.6.2 Simulation Execution**

- Procedure: `graphpop.simulate.run(config_path, n_replicates)` — external process via Java ProcessBuilder
- Re-imports as separate database or nodes with `simulated: true` label

**5.6.3 Observed vs Simulated Comparison**

- Procedure: `graphpop.simulate.compare(observed_run, simulated_run, statistics[])`

**5.6.4 ABC Framework**

- Procedure: `graphpop.abc.run(prior_ranges, n_simulations, summary_statistics[], tolerance)`
- Posterior distributions stored as InferenceResult nodes


### 5.7 LAYER 7: AI-Augmented Inference

**Status:** Long-term vision. Depends on all other layers.

**5.7.1 GraphRAG over the Variation Graph**

- Vector index on embedding properties (Variant + Sample embeddings)
- Cypher-based structured retrieval for relational queries
- Combined vector + graph retrieval for hybrid questions

Example: "What variants show signs of selection in drought-response pathways in indica rice?" → GraphRAG retrieves high-iHS variants → traverses IN_GENE → IN_PATHWAY → filters by GO term → augments with embedding analysis → generates ranked candidate list with evidence chains.

**5.7.2 Agentic Hypothesis-Testing Loop**

LangGraph agent with access to all Layer 2–6 procedures as tools:

1. **Observe:** query graph for anomalous patterns
2. **Hypothesize:** "selective sweep for drought adaptation in indica"
3. **Design:** propose SLiM simulation with specific selection coefficient
4. **Execute:** run simulation via Layer 6
5. **Evaluate:** compare simulated vs observed statistics
6. **Refine:** adjust parameters or reject hypothesis
7. **Report:** generate summary with evidence chain

**5.7.3 Natural Language Interface**

Translates natural language to Cypher + procedure calls:

- "Show me diversity across chromosome 1 for all populations" → `genome_scan()` + visualization
- "Which regions have extreme Tajima's D and what genes are there?" → windowed query + functional traversal
- "Compare this population's SFS to neutral expectation" → `sfs()` + simulation + comparison

---

## 6. Cross-Cutting Concerns

### 6.1 Analysis Versioning

Every materialized output carries a `run_id` property:

```
run_id = "{method}_{population}_{parameters_hash}_{timestamp}"
```

Example: `diversity_YRI_w100k_s50k_20260213T1430`

### 6.2 Storage Architecture — Fully Graph-Native (v3)

**No external columnar store.** All data lives in Neo4j:

- Population-level statistics → Variant node properties (allele count arrays)
- Individual genotypes → CARRIES relationships (sparse)
- Common-variant dense windows → DenseWindow nodes (on-demand materialization)
- Pairwise statistics → Sample–Sample edges (sparse, above threshold)

### 6.3 Benchmarking Framework

Built in from the start. Every procedure records wall-clock time and peak memory. Benchmark suite compares GraphPop to classical tools (vcftools, scikit-allel, PLINK, selscan, pixy, ∂a∂i) on identical data.

**Ten benchmark scenarios (A–J):** from simple windowed π to conditioned multi-signal selection scans. Scenario E (conditioned scan: Tajima's D only for missense variants) is the flagship demonstration of relational advantage.

### 6.4 Error Handling and Data Quality

- Missing data: -1 in dosage contexts; absent CARRIES edge is reference homozygote
- Quality metrics stored during import (call rate per variant, per sample)
- Procedures emit warnings for low-quality windows
- Post-import validation: orphaned nodes, NEXT chain completeness, allele frequency consistency

### 6.5 Extensibility Convention

All new methods follow: Java class → VectorOps primitives → Neo4j API → results as properties/edges → `@Procedure` registration with `graphpop.{module}.{method}` naming. Plugin is modular; third-party extensions follow same convention.

---

## 7. Why This Design Is Superior

### 7.1 Storage Scales Sub-Linearly with Sample Size

New sample → ~4.5M CARRIES edges + incremental allele count updates. Does NOT touch 1.5B Variant nodes.

### 7.2 Fast-Path Statistics: O(V × K) Independent of Sample Size

Whether 1,000 or 1,000,000 samples, allele count arrays are the same size (K elements for K populations). ~20,000× smaller I/O per node vs dosage vectors.

### 7.3 Rare Variants Are Efficiently Represented

>90% of All of Us variants are rare (AF < 0.01). Each rare variant has only a handful of CARRIES edges vs. a million-element vector that is >99% zeros.

### 7.4 The Graph Topology Becomes the Computation

```cypher
// Find samples carrying missense BRCA1 AND regulatory TP53 variants
MATCH (s:Sample)-[:CARRIES]->(v1:Variant)-[:HAS_CONSEQUENCE {type: 'missense'}]->(g1:Gene {symbol: 'BRCA1'})
MATCH (s)-[:CARRIES]->(v2:Variant)-[:HAS_CONSEQUENCE {type: 'regulatory_region'}]->(g2:Gene {symbol: 'TP53'})
RETURN s.id, v1.id, v2.id
```

### 7.5 Incremental Updates Are Trivial

- New sample → add CARRIES edges + update allele counts
- New annotations → add Gene/Pathway/GOTerm nodes + edges, no genotype changes

### 7.6 The CARRIES Bipartite Graph Is a Novel Computational Object

This bipartite graph IS what population genetics studies. Making it explicit as graph topology rather than implicit in a matrix is the conceptual core of the GraphPop paradigm.

---

## 8. Deployment Architecture — Three Tiers

### Tier 1: Development & Prototyping

- **Machine:** i9/64GB/RTX 4090, 4 TB NVMe
- **Data:** Single chromosome, 5,000–50,000 samples
- **Neo4j:** Community Edition 5.x on WSL2, page cache 20–30 GB
- **Graph:** ~50M variants, ~200M–1B CARRIES, 5–50 GB total (~140 GB for chr with 10k samples)

### Tier 2: Production Research

- **Cluster:** 3–5 machines, Enterprise + Infinigraph
- **Data:** Whole genome, 50,000–200,000 samples (e.g., single ancestry from All of Us)
- **Storage:** 10–50 TB, ~1B variants, ~200B–500B CARRIES
- **This is the "publishable paper" scale**

### Tier 3: Full All of Us Scale

- **Cluster:** 10–20+ machines, Enterprise Infinigraph
- **Data:** Whole genome, 500K–1M+ samples, all ancestries
- **Storage:** 100–300 TB, ~1.5B variants, ~4.5T CARRIES
- **Long-term target**

**Key Principle:** Schema, stored procedures, and query patterns are IDENTICAL across all tiers. Code on laptop works unchanged on cluster.

### Arithmetic Summary (Tier 3)

| Component | Scale | Storage | Feasible? |
|-----------|-------|---------|-----------|
| Variant nodes | 1.5B × ~2.5 KB | ~3.7 TB | Yes (Infinigraph) |
| NEXT edges | 1.5B × 24 bytes | ~36 GB | Easily |
| CARRIES edges | ~4.5T × 30 bytes | ~135 TB | Yes (cluster) |
| Sample nodes | 1M × ~2 KB | ~2 GB | Easily |
| LD edges | ~50B × 20 bytes | ~1 TB | Yes |
| Kinship edges | ~10B × 16 bytes | ~160 GB | Yes |
| Gene/Pathway/GO | ~1M nodes + 5M edges | ~1 GB | Trivially |
| GenomicWindow | ~3M × 200 bytes | ~600 MB | Trivially |
| **TOTAL** | | **~140 TB** | **Yes** |

---

## 9. Novel Scientific Contributions

### 9.1 Graph-Derived Summary Statistics

Community structure coefficients, embedding dimensionality, manifold curvature, centrality distributions — new summary statistics that classical population genetics does not have. The key question: do they contain information about evolutionary processes that SFS-based and LD-based statistics miss? Must be tested rigorously with simulation-based power analyses.

### 9.2 Conditioned Statistics as Default

Computing any statistic conditioned on functional annotation, recombination rate, or population structure in a single query. Making conditioned computation routine may change how the field interprets diversity patterns.

### 9.3 Embedding Geometry and Selection Regimes

- **Compression** for positive sweeps (reduced variance, lower local dimensionality)
- **Ridges** for balancing selection
- **Corridors** for purifying selection
- Rigorous test: simulations under known regimes + quantitative geometric statistics + power analysis

### 9.4 Relational Variant Graphs as Models of Evolutionary Constraint

Functional projection (variants → genes → pathways) may encode constraint topology. Variants with high betweenness centrality under stronger purifying selection. Testable by correlating graph centrality with conservation scores and dN/dS.

### 9.5 The CARRIES Bipartite Graph as Population-Genetic Object (v3 novel)

SFS = degree distribution of Variant nodes. Heterozygosity = degree distribution of Sample nodes. Kinship = shared neighborhoods. Population structure = bipartite community structure. This makes explicit what was always implicit.

---

## 10. Risk Assessment and Mitigation

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| 4.5T CARRIES edges infeasible in Infinigraph | Medium | High | Start at Tier 2 (200k samples, ~27 TB). Validate as Infinigraph matures. |
| Full-path operations (LD, kinship) too slow via CARRIES | Medium | Medium | DenseWindow materialization for common variants; pre-computed LD edges |
| Infinigraph property sharding doesn't keep traversal fast | Medium | Medium | CARRIES traversals local to topology shard; benchmark early |
| Embedding-based selection has no power beyond classical | Medium | High | Test early with simulation-based null models before investing further |
| Scope creep | High | High | Each phase is self-contained + publishable. Stop if foundations fail. |
| Java stored procedure performance disappointing | Low | Medium | SIMD benchmarked early. Fall back to C++ for bottlenecks. |
| Community adoption barrier (unfamiliar with Neo4j) | Medium | Medium | Docker + example datasets. Python wrapper API (graphpop-py). |
| Neo4j Enterprise + Infinigraph cost | Medium | Medium | Academic licensing. Portable design (Memgraph/TigerGraph hedge). |
| Phased haplotype data complexity | Medium | Low | CARRIES `phase` property; dedicated `(:Haplotype)-[:CONTAINS]->(:Variant)` for large-scale phased data |

---

## 11. Deliverables Summary

| Deliverable | Target Phase | Type |
|------------|-------------|------|
| graphpop-procedures.jar (Neo4j plugin) | Phase 1–4 | Software |
| graphpop-import (Python VCF→Neo4j pipeline) | Phase 1 | Software |
| graphpop-export (graph→simulator pipeline) | Phase 4 | Software |
| Benchmark suite + results | Phase 1 | Data + Paper 1 |
| Graph embedding for population structure | Phase 2 | Paper 2 |
| Manifold geometry as selection statistic | Phase 3 | Paper 3 (core scientific contribution) |
| Simulation-conditioned inference framework | Phase 4 | Paper 4 |
| AI-augmented evolutionary inference | Phase 5 | Paper 5 |
| Docker container + tutorial + example datasets | Phase 1 onward | Community resource |
| Python wrapper API (graphpop-py) | Phase 2 onward | Software |

---

## 12. The Competitive Landscape

### At Moderate Scale (≤10K samples)

- **vcftools / PLINK / scikit-allel / PopGenome**: Flat-file statistics. No relational queries. No functional conditioning. No persistent state.
- **tskit**: Brilliant for coalescent simulations. Not for empirical data without tree inference. No annotation layer.
- **vg / minigraph / PGGB**: Pan-genome sequence graphs. Different object (reference nodes vs variant nodes). Complementary.

### At All of Us Scale (100K–1M+ samples)

- **Hail (Broad Institute)**: Distributed genotype matrix on Spark. Powerful for GWAS. Not graph-native.
- **GVS (Genomic Variant Store)**: All of Us's variant storage. Optimized for joint calling. Not a graph database.
- **glow (Databricks)**: Distributed genomics on Spark. Matrix operations. Not graph-native.

**GraphPop occupies a unique niche:** the only system that natively represents the relational structure of population variation and computes evolutionary statistics directly from that structure. Hail can compute Tajima's D faster on a Spark cluster, but it cannot answer "find all individuals sharing an IBD segment overlapping a high-F_ST missense variant in a drought-response pathway" without building an ad-hoc pipeline for each such query.

No one is building this. That is both the opportunity and the risk.
