# slurm_prepare_csv.sh

## VCF to CSV Conversion for Neo4j Bulk Import (SLURM)

## Description

Converts a VCF file and population map into CSV files suitable for
`neo4j-admin import`. This is the first step of the two-step ingestion
pipeline. The job runs on any compute node and does not require a running
Neo4j instance. The operation is embarrassingly parallel by chromosome.

Scheduler: **SLURM**

## Usage

```bash
sbatch slurm_prepare_csv.sh <vcf> <population_map> [output_dir]
```

## Arguments

| Position | Name | Required | Default | Description |
|---|---|---|---|---|
| 1 | `vcf` | Yes | -- | Path to input VCF or VCF.gz file |
| 2 | `population_map` | Yes | -- | TSV mapping sample IDs to population labels |
| 3 | `output_dir` | No | `/scratch/$USER/graphpop_csv` | Directory for generated CSV files |

## SBATCH Directives

| Directive | Value | Description |
|---|---|---|
| `--job-name` | `graphpop-csv` | Job name in SLURM queue |
| `--cpus-per-task` | `16` | CPU cores for parallel processing |
| `--mem` | `64G` | Memory allocation |
| `--time` | `4:00:00` | Wall-time limit |
| `--output` | `graphpop_csv_%j.log` | Log file (`%j` = job ID) |

## Details

The script invokes `python -m graphpop_import.csv_emitter` which reads the
VCF file and population map, then writes the following CSV files to the
output directory:

- `variant_header.csv`, `variants.csv` -- Variant nodes with allele count arrays
- `sample_header.csv`, `samples.csv` -- Sample nodes
- `population_header.csv`, `populations.csv` -- Population nodes
- `chromosome_header.csv`, `chromosomes.csv` -- Chromosome nodes
- `window_header.csv`, `windows.csv` -- GenomicWindow nodes
- `next_header.csv`, `next_edges.csv` -- NEXT relationships (variant chain)
- `belongs_to_header.csv`, `belongs_to_edges.csv` -- BELONGS_TO relationships

Each CSV pair consists of a header file defining the Neo4j import schema and
a data file with the actual records. The population map should be a
tab-separated file with columns `sample_id` and `population`.

The thread count is automatically read from `$SLURM_CPUS_PER_TASK`. Upon
completion the script prints a directory listing of the generated files with
sizes.

## Examples

```bash
# Basic: 1000 Genomes chr22
sbatch slurm_prepare_csv.sh data/1kg_chr22.vcf.gz data/1kg_populations.tsv

# Specify output directory
sbatch slurm_prepare_csv.sh rice3k.vcf.gz rice_pops.tsv /scratch/$USER/rice_csv

# Request more time for large datasets
sbatch --time=8:00:00 --mem=128G slurm_prepare_csv.sh whole_genome.vcf.gz pops.tsv

# Chain with load step (submit load after CSV completes)
CSV_JOB=$(sbatch --parsable slurm_prepare_csv.sh data.vcf.gz pops.tsv)
sbatch --dependency=afterok:$CSV_JOB slurm_load_csv.sh /scratch/$USER/graphpop_csv ~/neo4j
```

## See Also

- [slurm-load-csv](slurm-load-csv.md) -- step 2: import CSVs into Neo4j
- [slurm-ingest-single](slurm-ingest-single.md) -- combined CSV + import
- [pbs-prepare-csv](pbs-prepare-csv.md) -- PBS equivalent
- [cluster-index](cluster-index.md) -- overview of all cluster scripts
