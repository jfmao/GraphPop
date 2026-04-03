"""graphpop batch — run any procedure across multiple populations and chromosomes."""
from __future__ import annotations

from pathlib import Path

import click

from ..cli import pass_ctx
from ..config import build_options_map, build_cypher
from ..formatters import format_output


# Map command names to their procedure and yield columns
COMMAND_REGISTRY = {
    "diversity": {
        "procedure": "graphpop.diversity",
        "args": lambda chr, pop, **kw: [f"'{chr}'", "1", "999999999", f"'{pop}'"],
        "yield": ["pi", "theta_w", "tajima_d", "fay_wu_h", "fay_wu_h_norm",
                   "het_exp", "het_obs", "fis", "n_variants", "n_segregating",
                   "n_polarized"],
    },
    "divergence": {
        "procedure": "graphpop.divergence",
        "args": lambda chr, pop, pop2=None, **kw: [
            f"'{chr}'", "1", "999999999", f"'{pop}'", f"'{pop2}'"
        ] if pop2 else None,
        "yield": ["fst_hudson", "fst_wc", "dxy", "da", "n_variants"],
    },
    "ihs": {
        "procedure": "graphpop.ihs",
        "args": lambda chr, pop, **kw: [f"'{chr}'", f"'{pop}'"],
        "yield": ["n_variants", "n_computed", "n_significant"],
    },
    "xpehh": {
        "procedure": "graphpop.xpehh",
        "args": lambda chr, pop, pop2=None, **kw: [
            f"'{chr}'", f"'{pop}'", f"'{pop2}'"
        ] if pop2 else None,
        "yield": ["n_variants", "n_computed", "n_significant"],
    },
    "nsl": {
        "procedure": "graphpop.nsl",
        "args": lambda chr, pop, **kw: [f"'{chr}'", f"'{pop}'"],
        "yield": ["n_variants", "n_computed", "n_significant"],
    },
    "sfs": {
        "procedure": "graphpop.sfs",
        "args": lambda chr, pop, **kw: [f"'{chr}'", "1", "999999999", f"'{pop}'"],
        "yield": ["sfs", "n_variants", "n_segregating"],
    },
    "roh": {
        "procedure": "graphpop.roh",
        "args": lambda chr, pop, **kw: [f"'{chr}'", f"'{pop}'"],
        "yield": ["n_samples", "mean_froh", "median_froh", "n_roh_segments"],
    },
    "garud-h": {
        "procedure": "graphpop.garud_h",
        "args": lambda chr, pop, **kw: [f"'{chr}'", f"'{pop}'"],
        "yield": ["n_windows", "mean_h12", "max_h12"],
    },
}


@click.command("batch")
@click.argument("command")
@click.option("--pops", required=True, help="Comma-separated population list")
@click.option("--chrs", required=True, help="Comma-separated chromosome list")
@click.option("--pop2", help="Second population (for divergence/xpehh, applied to all)")
@click.option("--workers", type=int, default=1,
              help="Parallel workers (default: 1, currently sequential)")
@click.option("-d", "--output-dir", required=True,
              type=click.Path(), help="Output directory (one TSV per pop-chr combo)")
@click.option("--persist", is_flag=True, help="Pass --persist to underlying command")
@click.option("--consequence", help="Filter by VEP consequence type")
@click.option("--pathway", help="Filter by pathway name")
@click.option("--format", "fmt", default="tsv", type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def batch(ctx, command, pops, chrs, pop2, workers, output_dir, persist,
          consequence, pathway, fmt):
    """Run a GraphPop procedure across multiple populations and/or chromosomes.

    COMMAND is the procedure name: diversity, divergence, ihs, xpehh, nsl,
    sfs, roh, garud-h.

    Creates one output file per (population, chromosome) combination in
    the output directory, named {command}_{pop}_{chr}.{ext}.

    \b
    Examples:
      graphpop batch diversity --pops EUR,AFR,EAS --chrs chr1,chr2,chr22 -d output/
      graphpop batch ihs --pops EUR,AFR --chrs chr22 --persist -d output/
      graphpop batch divergence --pops EUR --pop2 AFR --chrs chr22 -d output/
    """
    if command not in COMMAND_REGISTRY:
        available = ", ".join(sorted(COMMAND_REGISTRY.keys()))
        click.echo(f"Error: unknown command '{command}'. Available: {available}", err=True)
        raise SystemExit(1)

    spec = COMMAND_REGISTRY[command]
    pop_list = [p.strip() for p in pops.split(",")]
    chr_list = [c.strip() for c in chrs.split(",")]

    # Create output directory
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ext = fmt
    opts = build_options_map(consequence=consequence, pathway=pathway)
    if persist:
        opts["persist"] = True

    total = len(pop_list) * len(chr_list)
    completed = 0
    failed = 0

    for pop in pop_list:
        for chr_name in chr_list:
            completed += 1
            label = f"[{completed}/{total}] {command} {pop} {chr_name}"
            click.echo(f"{label} ...", err=True)

            try:
                args = spec["args"](chr=chr_name, pop=pop, pop2=pop2)
                if args is None:
                    click.echo(f"  Skipping: missing required argument (e.g., --pop2).", err=True)
                    failed += 1
                    continue

                cypher = build_cypher(
                    spec["procedure"],
                    args,
                    options=opts if opts else None,
                    yield_cols=spec["yield"],
                )
                records = ctx.run(cypher)

                if not records:
                    click.echo(f"  No results.", err=True)
                    failed += 1
                    continue

                out_file = out_dir / f"{command}_{pop}_{chr_name}.{ext}"
                format_output(records, str(out_file), fmt, f"batch {command}",
                              {"population": pop, "chr": chr_name})
                click.echo(f"  -> {out_file} ({len(records)} rows)", err=True)

            except SystemExit:
                click.echo(f"  FAILED (query error).", err=True)
                failed += 1
            except Exception as e:
                click.echo(f"  FAILED: {e}", err=True)
                failed += 1

    click.echo(f"\nBatch complete: {completed - failed}/{total} succeeded, "
               f"{failed} failed.", err=True)
