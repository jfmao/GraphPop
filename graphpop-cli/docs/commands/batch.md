# Run a command across multiple populations and chromosomes in parallel

## Description

`graphpop batch` executes any GraphPop procedure across multiple population and
chromosome combinations in parallel. Instead of writing shell loops to run
`graphpop ihs` for each population and chromosome, `batch` handles the
cross-product internally with configurable worker count and progress reporting.

The command constructs the cross-product of the provided populations and
chromosomes, then dispatches each combination as an independent job. Outputs are
written to a structured directory with `{command}_{pop}_{chr}` naming. Failed
jobs are logged and can be retried.

## Usage

```
graphpop batch COMMAND --pops POP1,POP2[,...] --chrs CHR1,CHR2[,...] [OPTIONS] [-- COMMAND_OPTIONS]
```

## Arguments

| Name | Type | Required | Default | Description |
|------|------|----------|---------|-------------|
| `COMMAND` | string | yes | -- | GraphPop subcommand to execute (e.g., `ihs`, `diversity`, `genome-scan`, `sfs`, `nsl`). |

## Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--pops` | string | required | Comma-separated list of population identifiers. |
| `--chrs` | string | required | Comma-separated list of chromosome names. |
| `--workers` | integer | 4 | Number of parallel worker processes. |
| `-d`, `--output-dir` | path | `./batch_results/` | Base directory for output files. Subdirectories are created per command. |
| `--retry` | integer | 0 | Number of retry attempts for failed jobs. |
| `--dry-run` | flag | false | Print the jobs that would be executed without running them. |
| `--format` | choice | `tsv` | Output format passed to each job. |

All options after `--` are passed through to the underlying command. For example,
`-- --min-af 0.05 --persist` passes those flags to each iHS invocation.

## Value

Each job produces an output file in `{output-dir}/{command}/{command}_{pop}_{chr}.{ext}`.
When all jobs complete, a summary is printed:

```
Batch complete: 30/30 jobs succeeded (0 failed)
Results in: ./batch_results/ihs/
  ihs_EUR_chr1.tsv   (2.4 MB, 3m12s)
  ihs_EUR_chr22.tsv  (0.8 MB, 1m05s)
  ...
```

If any jobs fail, a `failed_jobs.log` file is written to the output directory
with the command, error message, and stderr for each failure.

## Details

### Job scheduling

Jobs are dispatched to a worker pool. Each worker invokes `graphpop {command}`
as a subprocess with the appropriate population and chromosome arguments.
Workers are OS-level processes, allowing full Neo4j connection parallelism.

### Argument placement

The COMMAND argument determines how population and chromosome are passed:
- For commands with positional `CHR POPULATION` arguments (e.g., `ihs`,
  `genome-scan`), batch inserts them in the correct positions.
- For commands with `--pop` and `--chr` options, batch passes them as flags.
- Pass-through options (after `--`) are appended to every invocation.

### Resource management

Each worker opens its own Neo4j connection. With `--workers 4` and default
Neo4j connection pool settings, 4 concurrent queries execute against the
database. Increase workers for I/O-bound commands (diversity, sfs) and keep
workers low for memory-intensive commands (ihs, xpehh).

### Dry run

Use `--dry-run` to preview the full set of commands before execution. This is
recommended when constructing complex batch runs to verify argument placement.

## Examples

```bash
# Run iHS for 5 populations across 2 chromosomes
graphpop batch ihs --pops EUR,AFR,EAS,SAS,AMR --chrs chr1,chr22 \
    --workers 4 -d results/ -- --min-af 0.05 --persist

# Genome scan for all rice populations on chr1
graphpop batch genome-scan --pops GJ-tmp,GJ-trp,XI-1A,XI-1B,XI-2,XI-3 \
    --chrs chr1 --workers 6 -d rice_results/ -- 50000 25000

# Dry run to preview commands before executing
graphpop batch diversity --pops CEU,YRI,CHB --chrs chr1,chr22 \
    --dry-run -- 1 50818468

# SFS for multiple populations with retry on failure
graphpop batch sfs --pops CEU,YRI,CHB,GIH,JPT --chrs chr22 \
    --workers 3 --retry 2 -d sfs_results/

# Batch XP-EHH (pass --pop2 as a through option)
graphpop batch xpehh --pops EUR,EAS,SAS --chrs chr1,chr22 \
    --workers 3 -d xpehh_results/ -- --pop2 AFR --min-af 0.05 --persist
```

## See Also

- `graphpop run-all` -- orchestrated full analysis pipeline (predefined phases)
- `graphpop ihs` -- integrated haplotype score
- `graphpop genome-scan` -- sliding-window genome scan
- `graphpop xpehh` -- cross-population extended haplotype homozygosity
