# GraphPop CLI — User Manual

Version 0.1.0 | 60 commands across 11 categories

## Table of Contents

### Setup & Server
- [graphpop setup](commands/setup.md) — Download, configure, and initialize Neo4j
- [graphpop start](commands/start.md) — Start the Neo4j database server
- [graphpop stop](commands/stop.md) — Stop the Neo4j database server
- [graphpop status](commands/status.md) — Check server status and configuration

### Database Management
- [graphpop db list](commands/db-list.md) — List all databases
- [graphpop db create](commands/db-create.md) — Create a new database
- [graphpop db switch](commands/db-switch.md) — Set the active database
- [graphpop db drop](commands/db-drop.md) — Drop a database
- [graphpop db info](commands/db-info.md) — Show database details
- [graphpop import](commands/import.md) — Import VCF data into a graph database
- [graphpop dump](commands/dump.md) — Dump a database to file
- [graphpop load](commands/load.md) — Restore a database from dump

### Configuration & Validation
- [graphpop config](commands/config.md) — Manage configuration (init, show, set, path)
- [graphpop config init](commands/config-init.md) — Create config file
- [graphpop validate](commands/validate.md) — Check database integrity
- [graphpop inventory](commands/inventory.md) — Show what has been computed

### Population Genetics Procedures — FAST PATH
- [graphpop diversity](commands/diversity.md) — Nucleotide diversity, Tajima's D, Fay & Wu's H
- [graphpop divergence](commands/divergence.md) — Fst, Dxy, Da, PBS
- [graphpop sfs](commands/sfs.md) — Site frequency spectrum
- [graphpop joint-sfs](commands/joint-sfs.md) — Joint site frequency spectrum
- [graphpop genome-scan](commands/genome-scan.md) — Sliding-window genome scan
- [graphpop pop-summary](commands/pop-summary.md) — Population summary statistics

### Population Genetics Procedures — FULL PATH
- [graphpop ihs](commands/ihs.md) — Integrated haplotype score
- [graphpop xpehh](commands/xpehh.md) — Cross-population extended haplotype homozygosity
- [graphpop nsl](commands/nsl.md) — Number of segregating sites by length
- [graphpop roh](commands/roh.md) — Runs of homozygosity
- [graphpop garud-h](commands/garud-h.md) — Garud's haplotype homozygosity statistics
- [graphpop ld](commands/ld.md) — Linkage disequilibrium

### Conditioned Analysis & Filtering
- [graphpop filter](commands/filter.md) — Query persisted statistics with annotation filters

### Annotation Lookup & Exploration
- [graphpop lookup](commands/lookup.md) — Look up genes, pathways, variants, regions
- [graphpop neighbors](commands/neighbors.md) — Explore graph neighborhood around a gene

### Multi-Statistic Integration
- [graphpop converge](commands/converge.md) — Find convergent multi-statistic signals
- [graphpop rank-genes](commands/rank-genes.md) — Rank genes by composite selection evidence
- [graphpop compare](commands/compare.md) — Compare statistics between populations

### Orchestration & Export
- [graphpop run-all](commands/run-all.md) — Full-genome analysis
- [graphpop aggregate](commands/aggregate.md) — Generate summary tables
- [graphpop batch](commands/batch.md) — Run commands across multiple pops/chrs
- [graphpop report](commands/report.md) — Generate automated HTML report
- [graphpop export-windows](commands/export-windows.md) — Export GenomicWindow nodes

### Data Extraction
- [graphpop extract](commands/extract.md) — Extract variants, samples, or genotypes
- [graphpop export-bed](commands/export-bed.md) — Export regions as BED format

### Visualization
- [graphpop plot diversity-bar](commands/plot-diversity-bar.md) — Population diversity ranking
- [graphpop plot fst-heatmap](commands/plot-fst-heatmap.md) — Pairwise Fst matrix
- [graphpop plot manhattan](commands/plot-manhattan.md) — Genome-wide statistic scan
- [graphpop plot pinpis](commands/plot-pinpis.md) — piN/piS ratios
- [graphpop plot sfs-plot](commands/plot-sfs-plot.md) — Site frequency spectrum
- [graphpop plot roh-landscape](commands/plot-roh-landscape.md) — FROH distributions
- [graphpop plot gene-zoom](commands/plot-gene-zoom.md) — Multi-track regional view
- [graphpop plot chromosome](commands/plot-chromosome.md) — Chromosome-wide ideogram
- [graphpop plot pop-tree](commands/plot-pop-tree.md) — Population phylogeny
- [graphpop plot pca-scatter](commands/plot-pca-scatter.md) — PCA scatter plot
- [graphpop plot heatmap](commands/plot-heatmap.md) — General-purpose heatmap

### Utilities
- [graphpop query](commands/query.md) — Run arbitrary Cypher queries

### Programmatic Access (MCP Server)
- [graphpop-mcp](commands/mcp-server.md) — MCP server for AI agent access (21 tools)

### HPC Cluster Computing
- [Cluster Computing Guide](cluster-guide.md) — Deployment guide for SLURM/PBS clusters
- [Cluster Script Reference](commands/cluster/cluster-index.md) — All 11 job templates
  - SLURM: [setup](commands/cluster/slurm-setup-neo4j.md), [prepare-csv](commands/cluster/slurm-prepare-csv.md), [load-csv](commands/cluster/slurm-load-csv.md), [ingest-single](commands/cluster/slurm-ingest-single.md), [analysis](commands/cluster/slurm-analysis.md), [fullgenome-array](commands/cluster/slurm-fullgenome-array.md), [pairwise-array](commands/cluster/slurm-pairwise-array.md), [interactive](commands/cluster/slurm-interactive.md)
  - PBS: [prepare-csv](commands/cluster/pbs-prepare-csv.md), [analysis](commands/cluster/pbs-analysis.md), [fullgenome-array](commands/cluster/pbs-fullgenome-array.md)

## Vignettes
- [Rice 3K Genome Analysis](../vignettes/rice-3k-analysis.md) — Complete workflow from VCF to biological insight
- [Human 1000 Genomes Analysis](../vignettes/human-1000g-analysis.md) — Full-genome population genomics tutorial
