# pbs_prepare_csv.sh

## VCF to CSV Conversion for Neo4j Bulk Import (PBS/Torque)

## Description

Converts a VCF file and population map into CSV files suitable for
`neo4j-admin import`. This is the PBS/Torque equivalent of
`slurm_prepare_csv.sh`. The job runs on any compute node and does not
require a running Neo4j instance.

Scheduler: **PBS/Torque**

## Usage

```bash
qsub -v VCF_INPUT=/path/to/data.vcf.gz,POP_MAP=/path/to/pops.tsv pbs_prepare_csv.sh
qsub -v VCF_INPUT=data.vcf.gz,POP_MAP=pops.tsv,OUTPUT_DIR=/scratch/$USER/csv pbs_prepare_csv.sh
```

## Environment Variables (passed via `qsub -v`)

| Variable | Required | Default | Description |
|---|---|---|---|
| `VCF_INPUT` | **Yes** | -- | Path to input VCF or VCF.gz file |
| `POP_MAP` | **Yes** | -- | Path to TSV mapping sample IDs to populations |
| `OUTPUT_DIR` | No | `/scratch/$USER/graphpop_csv` | Directory for generated CSV files |
| `NCPUS` | No | `16` | Thread count (set by PBS) |

## PBS Directives

| Directive | Value | Description |
|---|---|---|
| `-N` | `graphpop-csv` | Job name |
| `-l ncpus` | `16` | CPU cores allocated |
| `-l mem` | `64GB` | Memory allocation |
| `-l walltime` | `4:00:00` | Wall-time limit |
| `-o` | `graphpop_csv.log` | Output log file |
| `-j oe` | -- | Merge stdout and stderr |

## Details

The script changes to the PBS submission directory (`$PBS_O_WORKDIR`) and
invokes `python -m graphpop_import.csv_emitter` with the VCF, population
map, and output directory. It produces the same CSV files as the SLURM
equivalent:

- Variant, Sample, Population, Chromosome, GenomicWindow node CSVs
- NEXT and BELONGS_TO relationship CSVs
- Corresponding header files for each CSV

Each CSV pair consists of a header file defining the Neo4j import schema
and a data file containing the records. The population map should be a
tab-separated file with `sample_id` and `population` columns.

### Differences from SLURM Version

- Arguments are passed via `qsub -v` environment variables instead of
  positional arguments.
- The script changes to `$PBS_O_WORKDIR` so relative paths work correctly.
- Thread count is read from `$NCPUS` (PBS) instead of
  `$SLURM_CPUS_PER_TASK`.

## Examples

```bash
# Basic: 1000 Genomes chr22
qsub -v VCF_INPUT=data/1kg_chr22.vcf.gz,POP_MAP=data/1kg_populations.tsv \
    pbs_prepare_csv.sh

# Custom output directory
qsub -v VCF_INPUT=rice3k.vcf.gz,POP_MAP=rice_pops.tsv,OUTPUT_DIR=/scratch/$USER/rice_csv \
    pbs_prepare_csv.sh

# Override resources for large dataset
qsub -v VCF_INPUT=whole_genome.vcf.gz,POP_MAP=pops.tsv \
    -l walltime=8:00:00,mem=128GB pbs_prepare_csv.sh

# Chain with analysis (PBS dependency)
CSV_JOB=$(qsub -v VCF_INPUT=data.vcf.gz,POP_MAP=pops.tsv pbs_prepare_csv.sh)
qsub -W depend=afterok:$CSV_JOB -v SCRIPT=scripts/full_analysis.py pbs_analysis.sh
```

## See Also

- [slurm-prepare-csv](slurm-prepare-csv.md) -- SLURM equivalent
- [pbs-analysis](pbs-analysis.md) -- run analyses under PBS
- [pbs-fullgenome-array](pbs-fullgenome-array.md) -- per-chromosome array job
- [cluster-index](cluster-index.md) -- overview of all cluster scripts
