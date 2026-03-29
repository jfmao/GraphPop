# pbs_analysis.sh

## Run Analysis Scripts as PBS Batch Jobs

## Description

Submits an arbitrary GraphPop analysis Python script as a PBS/Torque batch
job. Neo4j must already be running with the GraphPop database loaded and
procedures deployed. This is the PBS equivalent of `slurm_analysis.sh`.

Scheduler: **PBS/Torque**

## Usage

```bash
qsub -v SCRIPT=scripts/analysis.py pbs_analysis.sh
qsub -v SCRIPT=scripts/analysis.py,ARGS="--phase all" pbs_analysis.sh
```

## Environment Variables (passed via `qsub -v`)

| Variable | Required | Default | Description |
|---|---|---|---|
| `SCRIPT` | **Yes** | -- | Path to the Python analysis script |
| `ARGS` | No | (empty) | Arguments passed through to the script |
| `GRAPHPOP_URI` | No | -- | Neo4j bolt URI (if needed by the script) |
| `GRAPHPOP_PASSWORD` | No | -- | Neo4j password (if needed by the script) |

## PBS Directives

| Directive | Value | Description |
|---|---|---|
| `-N` | `graphpop-analysis` | Job name |
| `-l ncpus` | `8` | CPU cores allocated |
| `-l mem` | `64GB` | Memory allocation |
| `-l walltime` | `12:00:00` | Wall-time limit |
| `-o` | `graphpop_analysis.log` | Output log file |
| `-j oe` | -- | Merge stdout and stderr |

## Details

The script performs the following:

1. Changes to the PBS submission directory (`$PBS_O_WORKDIR`) so that
   relative paths in the `SCRIPT` variable resolve correctly.
2. Activates the Python/conda environment (uncomment the appropriate
   `module load` and `source activate` lines for your cluster).
3. Runs `python $SCRIPT $ARGS`.
4. Logs start and completion timestamps.

This script is intentionally minimal. It delegates all analysis logic to the
Python script being invoked. The script does not start or stop Neo4j; the
database must be accessible before job submission.

### Common Analysis Scripts

| Script | Purpose |
|---|---|
| `scripts/human_full_analysis.py` | Full human 1000G analysis pipeline |
| `scripts/rice_full_analysis.py` | Full rice 3K analysis pipeline |
| `scripts/human_interpret.py` | Interpretation of human results |

### Error Handling

The script runs with `set -euo pipefail`, so it exits immediately on any
error. Check the PBS log file for diagnostics. The `-j oe` directive merges
stdout and stderr into a single log for easier debugging.

### Differences from SLURM Version

- Script path and arguments are passed via `qsub -v` environment variables
  instead of positional arguments.
- The script changes to `$PBS_O_WORKDIR` for correct path resolution.
- Resource requests use PBS syntax (`-l ncpus`, `-l mem`, `-l walltime`).

## Examples

```bash
# Full human analysis
qsub -v SCRIPT=scripts/human_full_analysis.py,ARGS="--phase all" pbs_analysis.sh

# Rice analysis, specific phase
qsub -v SCRIPT=scripts/rice_full_analysis.py,ARGS="--phase diversity" pbs_analysis.sh

# Interpretation script
qsub -v SCRIPT=scripts/human_interpret.py,ARGS="pinsps" pbs_analysis.sh

# Override resources for long-running analysis
qsub -v SCRIPT=scripts/rice_full_analysis.py,ARGS="--phase all" \
    -l walltime=24:00:00,mem=128GB pbs_analysis.sh

# With Neo4j connection info
qsub -v SCRIPT=scripts/human_full_analysis.py,ARGS="--phase all",GRAPHPOP_URI=bolt://db:7687 \
    pbs_analysis.sh
```

## See Also

- [slurm-analysis](slurm-analysis.md) -- SLURM equivalent
- [pbs-prepare-csv](pbs-prepare-csv.md) -- VCF to CSV conversion
- [pbs-fullgenome-array](pbs-fullgenome-array.md) -- per-chromosome array job
- [cluster-index](cluster-index.md) -- overview of all cluster scripts
