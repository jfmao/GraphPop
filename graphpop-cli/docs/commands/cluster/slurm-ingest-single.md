# slurm_ingest_single.sh

## Combined CSV Generation and Neo4j Import in One Job (SLURM)

## Description

Runs the full ingestion pipeline -- VCF to CSV conversion, `neo4j-admin`
bulk import, plugin deployment, and Neo4j startup -- in a single SLURM job.
This is a convenience wrapper that combines the separate prepare-csv and
load-csv steps. Use this when you want a hands-off pipeline; use the
two-step approach when debugging or when CSV generation and loading run on
different nodes.

Scheduler: **SLURM**

## Usage

```bash
sbatch slurm_ingest_single.sh <vcf> <population_map> [neo4j_home] [data_dir] [csv_dir]
```

## Arguments

| Position | Name | Required | Default | Description |
|---|---|---|---|---|
| 1 | `vcf` | Yes | -- | Path to input VCF or VCF.gz file |
| 2 | `population_map` | Yes | -- | TSV mapping sample IDs to population labels |
| 3 | `neo4j_home` | No | `$HOME/neo4j` | Neo4j installation directory |
| 4 | `data_dir` | No | `/scratch/$USER/graphpop_db` | Target directory for Neo4j database files |
| 5 | `csv_dir` | No | `/scratch/$USER/graphpop_csv` | Intermediate directory for generated CSVs |

## SBATCH Directives

| Directive | Value | Description |
|---|---|---|
| `--job-name` | `graphpop-ingest` | Job name in SLURM queue |
| `--cpus-per-task` | `16` | CPU cores (used in CSV generation) |
| `--mem` | `128G` | Memory allocation |
| `--time` | `8:00:00` | Wall-time limit (longer for combined pipeline) |
| `--output` | `graphpop_ingest_%j.log` | Log file (`%j` = job ID) |

## Details

The script executes three sequential steps:

**Step 1: Generate CSVs.** Calls `graphpop_import.csv_emitter` to convert
the VCF and population map into Neo4j-compatible CSV files. Output goes to
the intermediate CSV directory.

**Step 2: Bulk import.** Runs `neo4j-admin database import full` with
`--overwrite-destination` to load all node and relationship CSVs. Before
importing, validates that the data directory is on a local filesystem (not
NFS/Lustre) and aborts if a network mount is detected.

**Step 3: Deploy and start.** Searches for the GraphPop procedures JAR in
the standard build location, copies it to the Neo4j plugins directory, and
starts Neo4j via `graphpop_import.cluster.start_neo4j`. Waits up to 120
seconds for the database to become available.

After completion, Neo4j is running and ready to accept analysis queries.
The imported node types are Variant, Sample, Population, Chromosome, and
GenomicWindow. Relationship types are NEXT and BELONGS_TO.

## Examples

```bash
# Basic: ingest 1000 Genomes chr22
sbatch slurm_ingest_single.sh data/1kg_chr22.vcf.gz data/1kg_populations.tsv

# Custom Neo4j and data locations
sbatch slurm_ingest_single.sh rice3k.vcf.gz rice_pops.tsv \
    /opt/neo4j /local/$USER/neo4j_db /scratch/$USER/rice_csv

# Override time for whole-genome dataset
sbatch --time=24:00:00 --mem=256G slurm_ingest_single.sh \
    whole_genome.vcf.gz pops.tsv

# Follow with analysis
INGEST_JOB=$(sbatch --parsable slurm_ingest_single.sh data.vcf.gz pops.tsv)
sbatch --dependency=afterok:$INGEST_JOB slurm_analysis.sh scripts/full_analysis.py --phase all
```

## See Also

- [slurm-prepare-csv](slurm-prepare-csv.md) -- standalone CSV step
- [slurm-load-csv](slurm-load-csv.md) -- standalone import step
- [slurm-setup-neo4j](slurm-setup-neo4j.md) -- initial Neo4j setup
- [cluster-index](cluster-index.md) -- overview of all cluster scripts
