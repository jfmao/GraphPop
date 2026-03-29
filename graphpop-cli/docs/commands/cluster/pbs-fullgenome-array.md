# pbs_fullgenome_array.sh

## Per-Chromosome Array Job for All Per-Population Procedures (PBS)

## Description

Runs the full suite of per-population GraphPop procedures across all
chromosomes using PBS array jobs. Each array task processes one chromosome,
iterating over all specified populations. Computes diversity (fast path),
iHS, nSL, Garud's H, and ROH (full path) for every population-chromosome
combination. This is the PBS/Torque equivalent of
`slurm_fullgenome_array.sh`.

Scheduler: **PBS/Torque**

## Usage

```bash
qsub -v GRAPHPOP_URI=bolt://db-node:7687,GRAPHPOP_PASSWORD=mypass pbs_fullgenome_array.sh
qsub -v GRAPHPOP_URI=bolt://db-node:7687,GRAPHPOP_POPULATIONS=GJ-tmp,GJ-trp pbs_fullgenome_array.sh
```

## Environment Variables (passed via `qsub -v`)

| Variable | Required | Default | Description |
|---|---|---|---|
| `GRAPHPOP_URI` | No | -- | Neo4j bolt URI |
| `GRAPHPOP_PASSWORD` | No | -- | Neo4j password |
| `GRAPHPOP_POPULATIONS` | No | All rice 3K populations | Comma-separated population IDs |
| `GRAPHPOP_OUTDIR` | No | `results` | Base output directory |
| `CHR_PREFIX` | No | `Chr` | Chromosome name prefix (`chr` for human, `Chr` for rice) |

## PBS Directives

| Directive | Value | Description |
|---|---|---|
| `-N` | `graphpop-fullgenome` | Job name |
| `-J` | `1-12` | Array range (one task per chromosome) |
| `-l ncpus` | `4` | CPU cores per task |
| `-l mem` | `16GB` | Memory per task |
| `-l walltime` | `4:00:00` | Wall-time per task |
| `-o` | `graphpop_fullgenome.log` | Output log file |
| `-j oe` | -- | Merge stdout and stderr |

## Details

Each array task reads `PBS_ARRAY_INDEX` to determine which chromosome to
process (e.g., index 1 becomes `Chr1`). It then iterates over all
populations and runs five procedures per population:

| Procedure | Path | Persists | Output |
|---|---|---|---|
| `graphpop diversity` | Fast | No | `results/diversity/{POP}_{CHR}.tsv` |
| `graphpop ihs` | Full | Yes | `results/ihs/{POP}_{CHR}.tsv` |
| `graphpop nsl` | Full | Yes | `results/nsl/{POP}_{CHR}.tsv` |
| `graphpop garud-h` | Full | No | `results/garud_h/{POP}_{CHR}.tsv` |
| `graphpop roh` | Full | No | `results/roh/{POP}_{CHR}.tsv` |

### Default Populations

When `GRAPHPOP_POPULATIONS` is not set, the script uses the full rice 3K
population list:

```
GJ-tmp, GJ-trp, XI-1A, XI-1B, XI-2, XI-3, XI-adm, cA-Aus, cB-Bas,
GJ-adm, GJ-sbtrp, admix
```

Override this for human data or to analyze a subset.

### Output Directory Structure

```
results/
  diversity/GJ-tmp_Chr1.tsv, GJ-tmp_Chr2.tsv, ...
  ihs/GJ-tmp_Chr1.tsv, ...
  nsl/GJ-tmp_Chr1.tsv, ...
  garud_h/GJ-tmp_Chr1.tsv, ...
  roh/GJ-tmp_Chr1.tsv, ...
```

### Differences from SLURM Version

- Array index uses `PBS_ARRAY_INDEX` instead of `SLURM_ARRAY_TASK_ID`.
- Array range is specified with `-J 1-12` instead of `--array=1-12`.
- Populations default to the rice 3K list rather than auto-detecting from
  the database (since `graphpop query` may not be available on all PBS
  clusters).
- Individual procedure failures are silently caught (`|| true`) rather than
  printing explicit warnings.

### Error Handling

Individual procedure failures are caught with `|| true`, allowing the job
to continue processing remaining populations and procedures. Check the PBS
log file for any error output from failed procedures.

## Examples

```bash
# Rice 3K: all populations, all 12 chromosomes
qsub -v GRAPHPOP_URI=bolt://db-node:7687,GRAPHPOP_PASSWORD=mypass \
    pbs_fullgenome_array.sh

# Specific populations only
qsub -v GRAPHPOP_URI=bolt://db:7687,GRAPHPOP_POPULATIONS="GJ-tmp,XI-1A" \
    pbs_fullgenome_array.sh

# Human 1000G: override prefix and array range
qsub -v GRAPHPOP_URI=bolt://db:7687,CHR_PREFIX=chr \
    -J 1-22 pbs_fullgenome_array.sh

# Custom output directory
qsub -v GRAPHPOP_URI=bolt://db:7687,GRAPHPOP_OUTDIR=/scratch/$USER/rice_results \
    pbs_fullgenome_array.sh

# Extended time for large chromosomes
qsub -v GRAPHPOP_URI=bolt://db:7687 -l walltime=8:00:00 pbs_fullgenome_array.sh
```

## See Also

- [slurm-fullgenome-array](slurm-fullgenome-array.md) -- SLURM equivalent
- [pbs-analysis](pbs-analysis.md) -- single-script analysis
- [pbs-prepare-csv](pbs-prepare-csv.md) -- VCF to CSV conversion
- [cluster-index](cluster-index.md) -- overview of all cluster scripts
