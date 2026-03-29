# slurm_pairwise_array.sh

## Per-Chromosome Array Job for Pairwise Statistics (SLURM)

## Description

Runs pairwise population statistics -- XP-EHH and divergence -- across all
chromosomes using SLURM array jobs. Each array task processes one chromosome,
iterating over all specified population pairs. XP-EHH results are persisted
back to the graph for downstream queries.

Scheduler: **SLURM**

## Usage

```bash
export GRAPHPOP_PAIRS="POP1:POP2,POP3:POP4"
sbatch --array=1-12 slurm_pairwise_array.sh
```

## SBATCH Directives

| Directive | Value | Description |
|---|---|---|
| `--job-name` | `graphpop-pairwise` | Job name in SLURM queue |
| `--array` | `1-12` | Array range (one task per chromosome) |
| `--cpus-per-task` | `4` | CPU cores per task |
| `--mem` | `32G` | Memory per task (higher than fullgenome) |
| `--time` | `8:00:00` | Wall-time per task |
| `--output` | `graphpop_pair_%A_%a.log` | Log file (`%A` = array job, `%a` = task) |

## Environment Variables

| Variable | Required | Default | Description |
|---|---|---|---|
| `GRAPHPOP_PAIRS` | **Yes** | -- | Population pairs, colon-separated within pair, comma-separated between pairs (e.g., `GJ-tmp:GJ-trp,XI-1A:XI-2`) |
| `GRAPHPOP_URI` | No | -- | Neo4j bolt URI |
| `GRAPHPOP_PASSWORD` | No | -- | Neo4j password |
| `GRAPHPOP_OUTDIR` | No | `results` | Base output directory |
| `CHR_PREFIX` | No | `Chr` | Chromosome name prefix |

## Details

Each array task constructs a chromosome name from `CHR_PREFIX` and
`SLURM_ARRAY_TASK_ID`. It then parses the `GRAPHPOP_PAIRS` variable and
runs two procedures for each population pair:

| Procedure | Path | Persists | Output |
|---|---|---|---|
| `graphpop xpehh` | Full | Yes | `results/xpehh/{POP1}_vs_{POP2}_{CHR}.tsv` |
| `graphpop divergence` | Fast | No | `results/divergence/{POP1}_vs_{POP2}_{CHR}.tsv` |

### Pair Specification Format

Pairs are specified as `POP1:POP2` with colons separating the two
populations within each pair. Multiple pairs are comma-separated:

```
GRAPHPOP_PAIRS="GJ-tmp:GJ-trp,GJ-tmp:XI-1A,XI-1A:cA-Aus"
```

### Output Directory Structure

```
results/
  xpehh/GJ-tmp_vs_GJ-trp_Chr1.tsv, ...
  divergence/GJ-tmp_vs_GJ-trp_Chr1.tsv, ...
```

### Error Handling

Individual procedure failures are caught and logged as warnings. The job
continues to the next pair/procedure rather than aborting. Memory is set
higher than the fullgenome array (32G vs 16G) because XP-EHH requires
loading haplotype data for two populations simultaneously.

### Resource Considerations

XP-EHH is a full-path procedure that traverses CARRIES relationships for
both populations, making it more memory- and time-intensive than single-
population statistics. For datasets with many pairs, consider submitting
separate jobs per pair to improve scheduling throughput.

## Examples

```bash
# Rice: three key population pairs across 12 chromosomes
export GRAPHPOP_URI=bolt://compute-042:7687
export GRAPHPOP_PAIRS="GJ-tmp:GJ-trp,GJ-tmp:XI-1A,XI-1A:cA-Aus"
sbatch --array=1-12 slurm_pairwise_array.sh

# Human: all 22 autosomes
export CHR_PREFIX=chr
export GRAPHPOP_PAIRS="CEU:YRI,CEU:CHB,YRI:CHB"
sbatch --array=1-22 slurm_pairwise_array.sh

# Single pair, specific chromosomes
export GRAPHPOP_PAIRS="GJ-tmp:XI-1A"
sbatch --array=1,5,10 slurm_pairwise_array.sh

# Custom output directory and extended time
export GRAPHPOP_OUTDIR=/scratch/$USER/pairwise_results
sbatch --array=1-12 --time=12:00:00 slurm_pairwise_array.sh
```

## See Also

- [slurm-fullgenome-array](slurm-fullgenome-array.md) -- per-population array job
- [slurm-analysis](slurm-analysis.md) -- single-script analysis
- [cluster-index](cluster-index.md) -- overview of all cluster scripts
- [cluster-guide.md](../../cluster-guide.md) -- full deployment guide
