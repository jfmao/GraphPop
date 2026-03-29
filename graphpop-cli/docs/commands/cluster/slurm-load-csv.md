# slurm_load_csv.sh

## Neo4j Bulk Import from Pre-Generated CSVs (SLURM)

## Description

Loads pre-generated CSV files into Neo4j using `neo4j-admin database import
full`. This is the second step of the two-step ingestion pipeline. The job
should run on a node with local SSD or fast scratch storage for the Neo4j
data directory. Requires Neo4j to be installed but not running (the import
operates directly on the data directory).

Scheduler: **SLURM**

## Usage

```bash
sbatch slurm_load_csv.sh <csv_dir> <neo4j_home> [data_dir]
```

## Arguments

| Position | Name | Required | Default | Description |
|---|---|---|---|---|
| 1 | `csv_dir` | Yes | -- | Directory containing generated CSV files |
| 2 | `neo4j_home` | Yes | -- | Neo4j installation directory |
| 3 | `data_dir` | No | `/scratch/$USER/graphpop_db` | Target directory for Neo4j database files |

## SBATCH Directives

| Directive | Value | Description |
|---|---|---|
| `--job-name` | `graphpop-load` | Job name in SLURM queue |
| `--cpus-per-task` | `4` | CPU cores allocated |
| `--mem` | `128G` | Memory allocation (high for import) |
| `--time` | `2:00:00` | Wall-time limit |
| `--output` | `graphpop_load_%j.log` | Log file (`%j` = job ID) |

## Details

The script performs the following steps:

1. Validates that the target data directory is on a local filesystem (not a
   network mount) using `graphpop_import.cluster.check_neo4j_data_dir`.
   The job aborts with an error if a network filesystem is detected, since
   Neo4j performance degrades severely on NFS/Lustre.
2. Creates the data directory if it does not exist.
3. Runs `neo4j-admin database import full` with `--overwrite-destination`,
   importing all node types (Variant, Sample, Population, Chromosome,
   GenomicWindow) and relationship types (NEXT, BELONGS_TO).

The import uses comma delimiters and semicolon array delimiters. The
`--overwrite-destination` flag replaces any existing database, so this is
not a resumable operation.

### Imported Node Types

| Label | Source Files |
|---|---|
| `Variant` | `variant_header.csv`, `variants.csv` |
| `Sample` | `sample_header.csv`, `samples.csv` |
| `Population` | `population_header.csv`, `populations.csv` |
| `Chromosome` | `chromosome_header.csv`, `chromosomes.csv` |
| `GenomicWindow` | `window_header.csv`, `windows.csv` |

### Imported Relationship Types

| Type | Source Files |
|---|---|
| `NEXT` | `next_header.csv`, `next_edges.csv` |
| `BELONGS_TO` | `belongs_to_header.csv`, `belongs_to_edges.csv` |

## Examples

```bash
# Basic import from scratch CSV directory
sbatch slurm_load_csv.sh /scratch/$USER/graphpop_csv ~/neo4j

# Custom data directory on local SSD
sbatch slurm_load_csv.sh /scratch/$USER/csv ~/neo4j /local/$USER/neo4j_db

# Chain after CSV generation
CSV_JOB=$(sbatch --parsable slurm_prepare_csv.sh data.vcf.gz pops.tsv)
sbatch --dependency=afterok:$CSV_JOB slurm_load_csv.sh /scratch/$USER/graphpop_csv ~/neo4j

# Request more memory for large datasets
sbatch --mem=256G slurm_load_csv.sh /scratch/$USER/csv ~/neo4j
```

## See Also

- [slurm-prepare-csv](slurm-prepare-csv.md) -- step 1: generate CSVs
- [slurm-ingest-single](slurm-ingest-single.md) -- combined CSV + import
- [slurm-setup-neo4j](slurm-setup-neo4j.md) -- Neo4j installation
- [cluster-index](cluster-index.md) -- overview of all cluster scripts
