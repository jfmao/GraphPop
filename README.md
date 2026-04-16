# GraphPop

[![PyPI](https://img.shields.io/pypi/v/graphpop-cli)](https://pypi.org/project/graphpop-cli/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/jfmao/GraphPop/blob/main/LICENSE)
[![Zenodo](https://zenodo.org/badge/DOI/10.5281/zenodo.19471963.svg)](https://doi.org/10.5281/zenodo.19471963)

**Graph-native population genomics with annotation-conditioned queries**

GraphPop is a population genomics platform that replaces the flat-file, matrix-based paradigm with a graph-native architecture. Variants, genotypes, genes, pathways, and computed statistics coexist as a single queryable structure within a Neo4j graph database. This enables three capabilities absent from classical tools: annotation-conditioned computation (any statistic restricted to a functional variant class in one parameter), a persistent analytical record (computed results stored as queryable node properties), and multi-statistic composition (cross-query independently computed statistics via graph traversal).

## Key Features

- **12 population genetics procedures** — diversity (π, θ_W, Tajima's D, Fay & Wu's H), divergence (Fst, D_xy, PBS), SFS, LD, iHS, XP-EHH, nSL, ROH, Garud's H, genome scan
- **Annotation conditioning** — any procedure conditioned on VEP consequence, Reactome pathway, or gene via a single `--consequence` / `--pathway` / `--gene` parameter
- **Persistent analytical record** — computed statistics written to graph nodes; multi-statistic convergence queries without re-computation
- **62-command CLI** — complete command-line interface with default TSV output, no Python/Neo4j/Cypher knowledge required
- **11 visualization types** — publication-ready figures following Nature Methods guidelines
- **MCP server** — 21 tools for AI agent access via the Model Context Protocol
- **63–327× speedup** over scikit-allel; constant ~160 MB memory

## Quick Start

```bash
# Install the CLI
cd graphpop-cli && pip install -e .

# Set up Neo4j (downloads, configures, and starts automatically)
graphpop setup --password mypassword
graphpop start

# Import your data
graphpop import --vcf data.vcf.gz --panel panel.txt --database myproject

# Run an analysis
graphpop diversity chr1 1 50000000 EUR -o diversity.tsv

# Annotation-conditioned analysis (piN)
graphpop diversity chr1 1 50000000 EUR --consequence missense_variant -o piN.tsv

# Selection scan with persistent results
graphpop ihs chr22 EUR --persist -o ihs.tsv

# Multi-statistic convergence detection
graphpop converge --stats ihs,xpehh,h12 --thresholds 2.0,2.0,0.3 --pop EUR

# Rank genes by composite selection evidence
graphpop rank-genes --pop EUR --top 50 -o top_genes.tsv

# Visualize
graphpop plot manhattan ihs.tsv --stat ihs --threshold 2.5 -o fig_ihs.png
graphpop plot diversity-bar results/diversity/ -o fig_diversity.png

# Generate automated report
graphpop report -o analysis_report.html
```

## Architecture

GraphPop organises data as a labelled property graph:

```
Variant ──HAS_CONSEQUENCE──> Gene ──IN_PATHWAY──> Pathway
   │                          │
   ├── ac[], an[], af[]        ├── IN_PATHWAY ──> GOTerm
   ├── gt_packed (2 bits/sample)
   ├── ihs_EUR, xpehh_EUR_AFR  (persisted selection statistics)
   │
   ├──NEXT──> Variant  (genomic order)
   └──ON_CHROMOSOME──> Chromosome
```

Two computational paths:

| Path | Procedures | Data source | Complexity |
|------|-----------|------------|------------|
| **FAST PATH** | diversity, divergence, SFS, genome scan, pop summary | Pre-aggregated allele count arrays | O(V × K), independent of sample count |
| **FULL PATH** | LD, iHS, XP-EHH, nSL, ROH, Garud's H | Bit-packed haplotype matrices (1 bit/haplotype) | O(V × N), 87% memory reduction |

## Installation

### Quick install (recommended)

```bash
pip install graphpop-cli                  # Install the CLI
graphpop setup --password mypassword      # Downloads Neo4j + procedures plugin
graphpop start                            # Start the database
graphpop doctor                           # Verify installation health
```

**Prerequisites:** Python 3.10+ and Java 21+ (required by Neo4j runtime).
Java can be installed via `conda install -c conda-forge openjdk=21` (no admin
privileges needed). The `graphpop setup` command automatically downloads
Neo4j Community Edition and the pre-compiled GraphPop procedures plugin from
GitHub Releases. No Maven or admin privileges are required.

For conda install, development install, Docker, HPC cluster deployment,
using an existing Neo4j installation, offline/air-gapped install, and
troubleshooting, see the **[full Installation Guide](docs/INSTALL.md)**.

## CLI Command Overview

62 commands across 11 categories:

```
graphpop
├── Setup & Server        setup, start, stop, status
├── Database Management   db {list,create,switch,drop,info}, import, dump, load
├── Configuration         config {init,show,set,path}, validate, inventory
├── FAST PATH Procedures  diversity, divergence, sfs, joint-sfs, genome-scan, pop-summary
├── FULL PATH Procedures  ihs, xpehh, nsl, roh, garud-h, ld
├── Conditioned Analysis  filter (post-hoc annotation filtering)
├── Annotation Lookup     lookup {gene,pathway,variant,region}, neighbors
├── Multi-Stat Integration converge, rank-genes, compare
├── Orchestration         run-all, aggregate, batch, report, export-windows
├── Data Extraction       extract {variants,samples,genotypes}, export-bed
├── Visualization         plot {diversity-bar,fst-heatmap,manhattan,pinpis,
│                               sfs-plot,roh-landscape,gene-zoom,chromosome,
│                               pop-tree,pca-scatter,heatmap}
└── Utilities             query
```

See [graphpop-cli/COMMANDS.md](graphpop-cli/COMMANDS.md) for the full command tree with all options.

## Annotation Conditioning

GraphPop's signature capability: condition any FAST PATH statistic on functional annotations via graph traversal.

```bash
# Standard diversity
graphpop diversity chr1 1 43270923 GJ-tmp

# Missense-only diversity (piN) — one parameter change
graphpop diversity chr1 1 43270923 GJ-tmp --consequence missense_variant

# Pathway-restricted genome scan
graphpop genome-scan chr1 GJ-tmp 100000 50000 --pathway "Starch biosynthesis"

# For FULL PATH statistics, use compute-then-filter:
graphpop ihs chr1 GJ-tmp --persist
graphpop filter ihs chr1 GJ-tmp --consequence missense_variant --min-score 2.0
```

## Example: End-to-End Rice 3K Analysis

```bash
# Import
graphpop import --vcf NB_final_snp.vcf.gz --panel rice_panel.txt \
    --database rice3k --vep rice_snpeff.vcf --pathways plant_reactome.tsv

# Full-genome analysis (all 12 pops × 12 chromosomes)
graphpop run-all -d results/

# Cost of domestication (piN/piS across all subpopulations)
for pop in GJ-tmp GJ-trp XI-1A XI-1B cA-Aus; do
    graphpop diversity Chr1 1 43270923 $pop --consequence missense_variant
    graphpop diversity Chr1 1 43270923 $pop --consequence synonymous_variant
done

# Multi-statistic convergence at domestication loci
graphpop converge --stats ihs,xpehh,h12 --thresholds 2.0,2.0,0.3 --pop GJ-tmp

# Gene ranking
graphpop rank-genes --pop GJ-tmp --top 50 -o top_genes.tsv

# Explore a candidate gene
graphpop lookup gene GW5
graphpop plot gene-zoom GW5 --pop GJ-tmp -o fig_gw5.pdf

# Summary
graphpop aggregate -d results/ -o tables/
graphpop report -o rice3k_report.html
graphpop dump --database rice3k -o rice3k_v1.dump
```

## Documentation

| Resource | Description |
|----------|-------------|
| [Command reference](graphpop-cli/docs/index.md) | Per-command manuals (53 files, R-style format) |
| [Rice 3K vignette](graphpop-cli/vignettes/rice-3k-analysis.md) | End-to-end tutorial: VCF → biological insight |
| [Human 1000G vignette](graphpop-cli/vignettes/human-1000g-analysis.md) | Full-genome population genomics tutorial |
| [Command tree](graphpop-cli/COMMANDS.md) | All 62 commands with options and examples |
| [MCP server](graphpop-cli/docs/commands/mcp-server.md) | AI agent access via Model Context Protocol |
| [Design document](docs/GraphPop_Compiled_Holistic_Design.md) | Architecture and implementation details |

## Pre-built Databases

Pre-built Neo4j graph databases with all computed statistics are available for download:

| Dataset | Samples | Variants | Size | Download |
|---------|---------|----------|------|----------|
| Human 1000 Genomes (22 autosomes) | 3,202 | 70.7M | ~31 GB | [Zenodo](https://doi.org/10.5281/zenodo.19472010) |
| Rice 3K Genomes (12 chromosomes) | 3,024 | 29.6M | ~15 GB | [Zenodo](https://doi.org/10.5281/zenodo.19471968) |

Load a pre-built database:

```bash
graphpop load --dump-file rice3k_v1.dump --database rice3k
graphpop start
graphpop db info
```

## Repository Structure

```
GraphPop/
├── graphpop-procedures/    Java compute engine (12 Neo4j stored procedures)
│   ├── src/main/java/      Procedure implementations + VectorOps SIMD kernels
│   └── src/test/java/      147 unit tests
├── graphpop-cli/           Python CLI (pip install -e .)
│   ├── src/graphpop_cli/   62 commands across 11 categories
│   ├── docs/commands/      53 per-command manuals
│   └── vignettes/          2 end-to-end tutorials (rice 3K + human 1000G)
├── graphpop-mcp/           MCP server for AI agent access
│   ├── src/graphpop_mcp/   21 tools (procedures + analytical + lookup)
│   └── tooluniverse/       ToolUniverse registration
├── graphpop-import/        VCF-to-graph import pipeline
│   └── src/graphpop_import/  VCF parser, CSV emitter, annotation loaders
├── scripts/                Analysis and utility scripts
│   ├── rice_investigate_*  21 rice deep investigation scripts
│   ├── benchmark-*.py      Performance benchmarking
│   └── cluster/            HPC job scripts (SLURM/PBS)
├── docs/                   Design documentation
└── benchmarks/             Benchmark descriptions
```

## Performance

Benchmarked on 1000 Genomes chr22 (3,202 samples, 1,066,555 variants):

| Statistic | GraphPop | scikit-allel | Speedup |
|-----------|----------|-------------|---------|
| π / θ_W / Tajima's D | 12.1 s | 1,757 s | **146×** |
| Fst / D_xy | 8.8 s | 2,165 s | **245×** |
| SFS | 5.9 s | 1,937 s | **327×** |
| iHS | 10.9 s | 1,948 s | **179×** |
| XP-EHH | 16.7 s | 2,280 s | **136×** |

Memory: GraphPop Neo4j path maintains constant ~160 MB regardless of analysis type.

## Data Availability

Pre-built Neo4j databases and analysis data are available on Zenodo:

| Resource | DOI | Size |
|----------|-----|------|
| Analysis results, source data, scripts | [10.5281/zenodo.19471963](https://doi.org/10.5281/zenodo.19471963) | 29 MB |
| Rice 3K pre-built database | [10.5281/zenodo.19471968](https://doi.org/10.5281/zenodo.19471968) | 14 GB |
| Human 1000G pre-built database | [10.5281/zenodo.19472010](https://doi.org/10.5281/zenodo.19472010) | 31 GB |

## Sister Project: GraphMana

[**GraphMana**](https://github.com/jfmao/GraphMana) is the companion data management platform for population genomics projects. While GraphPop provides the compute engine (12 stored procedures, annotation conditioning, multi-statistic composition), GraphMana handles the project lifecycle: incremental sample addition, cohort management, data provenance, multi-format export (17 formats including TreeMix, EIGENSOFT, BayPass), and quality control.

> Estaji, E., Zhao, S.-W., Chen, Z.-Y., Nie, S. & Mao, J.-F. GraphMana: graph-native data management for population genomics projects. *bioRxiv* (2026). [doi:10.64898/2026.04.11.717925](https://doi.org/10.64898/2026.04.11.717925)

## Citation

If you use GraphPop in your research, please cite:

> Mao, J.-F. GraphPop: graph-native computation decouples population genomics complexity from sample count. *Preprint* (2026).

Software: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19471963.svg)](https://doi.org/10.5281/zenodo.19471963)

## License

MIT License. See [LICENSE](LICENSE) for details.

## Contributing

Contributions are welcome. Please open an issue for bug reports or feature requests, or submit a pull request for code contributions.
