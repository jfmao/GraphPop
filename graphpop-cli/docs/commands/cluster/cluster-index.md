# GraphPop HPC Cluster Scripts

## Description

Collection of job submission scripts for running GraphPop workflows on
high-performance computing clusters. Scripts are provided for both SLURM and
PBS/Torque schedulers, covering the full pipeline from VCF ingestion through
genome-wide population genetics analysis.

## Workflow Overview

A typical cluster workflow proceeds in three phases:

1. **Setup** -- install Neo4j and deploy the GraphPop plugin (once per cluster).
2. **Ingest** -- convert VCF to CSV, then bulk-import into Neo4j.
3. **Analysis** -- run per-population and pairwise statistics as array jobs.

## Script Inventory

| Script | Scheduler | Purpose | Doc |
|---|---|---|---|
| `slurm_setup_neo4j.sh` | SLURM | One-time Neo4j installation and configuration | [slurm-setup-neo4j](slurm-setup-neo4j.md) |
| `slurm_prepare_csv.sh` | SLURM | VCF to CSV conversion (no database needed) | [slurm-prepare-csv](slurm-prepare-csv.md) |
| `slurm_load_csv.sh` | SLURM | Bulk import CSVs into Neo4j | [slurm-load-csv](slurm-load-csv.md) |
| `slurm_ingest_single.sh` | SLURM | Combined CSV generation + import in one job | [slurm-ingest-single](slurm-ingest-single.md) |
| `slurm_analysis.sh` | SLURM | Run arbitrary analysis scripts as batch jobs | [slurm-analysis](slurm-analysis.md) |
| `slurm_fullgenome_array.sh` | SLURM | Per-chromosome array job for population statistics | [slurm-fullgenome-array](slurm-fullgenome-array.md) |
| `slurm_pairwise_array.sh` | SLURM | Per-chromosome array job for pairwise statistics | [slurm-pairwise-array](slurm-pairwise-array.md) |
| `slurm_interactive.sh` | SLURM | Interactive session with Neo4j (source, do not submit) | [slurm-interactive](slurm-interactive.md) |
| `pbs_prepare_csv.sh` | PBS | VCF to CSV conversion | [pbs-prepare-csv](pbs-prepare-csv.md) |
| `pbs_analysis.sh` | PBS | Run analysis scripts as batch jobs | [pbs-analysis](pbs-analysis.md) |
| `pbs_fullgenome_array.sh` | PBS | Per-chromosome array job for population statistics | [pbs-fullgenome-array](pbs-fullgenome-array.md) |

## Quick Start (SLURM)

```bash
# 1. Setup Neo4j (once)
sbatch slurm_setup_neo4j.sh

# 2. Prepare CSVs from VCF
sbatch slurm_prepare_csv.sh data.vcf.gz populations.tsv

# 3. Load into Neo4j
sbatch slurm_load_csv.sh /scratch/$USER/graphpop_csv ~/neo4j

# 4. Run genome-wide analysis
export GRAPHPOP_URI=bolt://db-node:7687
sbatch --array=1-22 slurm_fullgenome_array.sh
```

## Quick Start (PBS)

```bash
# 1. Prepare CSVs
qsub -v VCF_INPUT=data.vcf.gz,POP_MAP=populations.tsv pbs_prepare_csv.sh

# 2. Run genome-wide analysis
qsub -v GRAPHPOP_URI=bolt://db-node:7687,GRAPHPOP_PASSWORD=mypass pbs_fullgenome_array.sh
```

## Environment Setup

All scripts expect a conda environment named `graphpop` with Python 3.11+,
neo4j-driver, cyvcf2, and the graphpop packages installed. Java 21+ is required
for Neo4j and is loaded via `module load java/21` where available.

## See Also

- [cluster-guide.md](../../cluster-guide.md) -- full cluster deployment guide
- [GraphPop CLI reference](../README.md) -- command-line interface documentation
