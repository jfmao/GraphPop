# slurm_analysis.sh

## Run Analysis Scripts as SLURM Batch Jobs

## Description

Submits an arbitrary GraphPop analysis Python script as a SLURM batch job.
Neo4j must already be running with the GraphPop database loaded and
procedures deployed. This is a thin wrapper that activates the environment
and passes the script and arguments to the Python interpreter.

Scheduler: **SLURM**

## Usage

```bash
sbatch slurm_analysis.sh <analysis_script> [args...]
```

## Arguments

| Position | Name | Required | Default | Description |
|---|---|---|---|---|
| 1 | `analysis_script` | Yes | -- | Path to the Python analysis script |
| 2+ | `args` | No | -- | Arguments passed through to the analysis script |

## SBATCH Directives

| Directive | Value | Description |
|---|---|---|
| `--job-name` | `graphpop-analysis` | Job name in SLURM queue |
| `--cpus-per-task` | `8` | CPU cores allocated |
| `--mem` | `64G` | Memory allocation |
| `--time` | `12:00:00` | Wall-time limit |
| `--output` | `graphpop_analysis_%j.log` | Log file (`%j` = job ID) |

## Details

The script performs the following:

1. Activates the Python/conda environment (uncomment the appropriate
   `module load` and `source activate` lines for your cluster).
2. Runs `python <analysis_script> <args...>`.
3. Logs start and completion timestamps.

This script is intentionally minimal and generic. It delegates all logic to
the Python analysis script being invoked. The script does not start or stop
Neo4j; the database must be available before submission.

### Common Analysis Scripts

| Script | Purpose |
|---|---|
| `scripts/human_full_analysis.py` | Full human 1000G analysis pipeline |
| `scripts/rice_full_analysis.py` | Full rice 3K analysis pipeline |
| `scripts/human_interpret.py` | Interpretation of human results |
| `scripts/rice_investigate_*.py` | Individual rice investigation scripts |

### Error Handling

The script runs with `set -euo pipefail`, so it will exit immediately on
any error. Check the SLURM log file for diagnostics. Common failures
include Neo4j not running, incorrect `GRAPHPOP_URI`, or missing Python
dependencies.

## Examples

```bash
# Full human analysis
sbatch slurm_analysis.sh scripts/human_full_analysis.py --phase all

# Rice analysis, specific phase
sbatch slurm_analysis.sh scripts/rice_full_analysis.py --phase diversity

# Interpretation script with subcommand
sbatch slurm_analysis.sh scripts/human_interpret.py pinsps

# Override resources for long-running analysis
sbatch --time=24:00:00 --mem=128G slurm_analysis.sh scripts/rice_full_analysis.py --phase all

# Chain after ingestion
INGEST_JOB=$(sbatch --parsable slurm_ingest_single.sh data.vcf.gz pops.tsv)
sbatch --dependency=afterok:$INGEST_JOB slurm_analysis.sh scripts/human_full_analysis.py --phase all
```

## See Also

- [slurm-fullgenome-array](slurm-fullgenome-array.md) -- per-chromosome array analysis
- [slurm-pairwise-array](slurm-pairwise-array.md) -- pairwise statistics
- [pbs-analysis](pbs-analysis.md) -- PBS equivalent
- [cluster-index](cluster-index.md) -- overview of all cluster scripts
