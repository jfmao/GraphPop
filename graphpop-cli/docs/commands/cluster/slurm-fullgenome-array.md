# slurm_fullgenome_array.sh

## Per-Chromosome Array Job for All Per-Population Procedures (SLURM)

## Description

Runs the full suite of per-population GraphPop procedures across all
chromosomes using SLURM array jobs. Each array task processes one
chromosome, iterating over all specified populations. Computes diversity
(fast path), iHS, nSL, Garud's H, and ROH (full path) for every
population-chromosome combination.

Scheduler: **SLURM**

## Usage

```bash
sbatch --array=1-12 slurm_fullgenome_array.sh [population]
sbatch --array=1-22 slurm_fullgenome_array.sh [population]
```

## Arguments

| Position | Name | Required | Default | Description |
|---|---|---|---|---|
| 1 | `population` | No | Auto-detected from database | Comma-separated list of population IDs, or single population |

## SBATCH Directives

| Directive | Value | Description |
|---|---|---|
| `--job-name` | `graphpop-fullgenome` | Job name in SLURM queue |
| `--array` | `1-12` | Array range (adjust per organism) |
| `--cpus-per-task` | `4` | CPU cores per chromosome task |
| `--mem` | `16G` | Memory per task |
| `--time` | `4:00:00` | Wall-time per task |
| `--output` | `graphpop_%A_%a.log` | Log file (`%A` = array job, `%a` = task) |

## Environment Variables

| Variable | Default | Description |
|---|---|---|
| `GRAPHPOP_URI` | -- | Neo4j bolt URI (e.g., `bolt://db-node:7687`) |
| `GRAPHPOP_PASSWORD` | -- | Neo4j password |
| `GRAPHPOP_DATABASE` | -- | Neo4j database name (e.g., `rice3k`) |
| `GRAPHPOP_POPULATIONS` | Auto-detected | Comma-separated population IDs |
| `GRAPHPOP_OUTDIR` | `results` | Base output directory |
| `CHR_PREFIX` | `Chr` | Chromosome name prefix (`chr` for human, `Chr` for rice) |

## Details

Each array task constructs a chromosome name from `CHR_PREFIX` and
`SLURM_ARRAY_TASK_ID` (e.g., `Chr1`, `Chr2`, ..., `Chr12`). It then loops
over all populations and runs five procedures per population:

| Procedure | Path | Persists | Output |
|---|---|---|---|
| `graphpop diversity` | Fast | No | `results/diversity/{POP}_{CHR}.tsv` |
| `graphpop ihs` | Full | Yes | `results/ihs/{POP}_{CHR}.tsv` |
| `graphpop nsl` | Full | Yes | `results/nsl/{POP}_{CHR}.tsv` |
| `graphpop garud-h` | Full | No | `results/garud_h/{POP}_{CHR}.tsv` |
| `graphpop roh` | Full | No | `results/roh/{POP}_{CHR}.tsv` |

### Output Directory Structure

```
results/
  diversity/GJ-tmp_Chr1.tsv, GJ-tmp_Chr2.tsv, ...
  ihs/GJ-tmp_Chr1.tsv, ...
  nsl/GJ-tmp_Chr1.tsv, ...
  garud_h/GJ-tmp_Chr1.tsv, ...
  roh/GJ-tmp_Chr1.tsv, ...
```

### Population Auto-Detection

If no populations are specified via argument or `GRAPHPOP_POPULATIONS`,
the script queries the database for all populations with more than one
sample using `graphpop query`. The job exits with an error if detection
fails.

### Error Handling

Individual procedure failures are caught and logged as warnings (e.g.,
`WARN: ihs GJ-tmp Chr5 failed`). The job continues to the next
procedure/population rather than aborting, ensuring maximum throughput.

## Examples

```bash
# Rice 3K: all 12 chromosomes, auto-detect populations
export GRAPHPOP_URI=bolt://compute-042:7687
export GRAPHPOP_PASSWORD=graphpop
sbatch --array=1-12 slurm_fullgenome_array.sh

# Human 1000G: 22 autosomes with chr prefix
export CHR_PREFIX=chr
sbatch --array=1-22 slurm_fullgenome_array.sh

# Single population only
sbatch --array=1-12 slurm_fullgenome_array.sh GJ-tmp

# Custom output and specific chromosomes
export GRAPHPOP_OUTDIR=/scratch/$USER/rice_results
sbatch --array=1,3,5 slurm_fullgenome_array.sh

# Multiple specific populations
export GRAPHPOP_POPULATIONS="GJ-tmp,GJ-trp,XI-1A"
sbatch --array=1-12 slurm_fullgenome_array.sh
```

## See Also

- [slurm-pairwise-array](slurm-pairwise-array.md) -- pairwise statistics array job
- [slurm-analysis](slurm-analysis.md) -- single-script analysis
- [pbs-fullgenome-array](pbs-fullgenome-array.md) -- PBS equivalent
- [cluster-index](cluster-index.md) -- overview of all cluster scripts
