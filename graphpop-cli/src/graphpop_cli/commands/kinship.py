"""graphpop kinship — pairwise kinship suite (KING-robust + branch GRM)."""
from __future__ import annotations

import click

from ..cli import pass_ctx
from ..config import build_cypher
from ..formatters import format_output


@click.group()
def kinship():
    """Pairwise kinship procedures.

    Subcommands:
      king           KING-robust on packed genotypes (validation baseline)
      bgrm           Branch GRM via ARG traversal (Phase 4 ARG layer required)
      bgrm-posterior Posterior-aware branch GRM (over multiple ARG runs)
    """


@kinship.command("king")
@click.argument("chr")
@click.argument("population")
@click.option("--start", type=int, default=1, show_default=True,
              help="Start position in bp")
@click.option("--end", type=int, default=None,
              help="End position in bp (default: end of chromosome)")
@click.option("--min-snp", type=int, default=1000, show_default=True,
              help="Skip pairs with fewer informative variants")
@click.option("--min-phi", type=float, default=-0.5, show_default=True,
              help="Skip pairs with kinship below this threshold")
@click.option("--include-self", is_flag=True,
              help="Emit (s, s) self-pair rows")
@click.option("--samples", multiple=True,
              help="Restrict to a subset of sample IDs (overrides population)")
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def king(ctx, chr, population, start, end, min_snp, min_phi, include_self,
         samples, output_path, fmt):
    """Compute KING-robust pairwise kinship for samples in POPULATION on CHR.

    KING-robust (Manichaikul et al. 2010) is the matrix-based validation
    baseline of the GraphPop kinship suite. For graph-native, ARG-based
    kinship use ``graphpop kinship bgrm`` instead.
    """
    opts: dict[str, object] = {
        "min_snp": min_snp,
        "min_phi": min_phi,
        "include_self": include_self,
        "start": start,
    }
    if end is not None:
        opts["end"] = end
    if samples:
        opts["samples"] = list(samples)

    cypher = build_cypher(
        "graphpop.kinship.king",
        [f"'{chr}'", f"'{population}'"],
        options=opts,
        yield_cols=["sample_a", "sample_b", "phi", "ibs0", "het_het",
                    "n_snp", "n_aa_min", "method"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "kinship-king",
                  {"chr": chr, "pop": population})


@kinship.command("bgrm")
@click.argument("run_id")
@click.option("--start", type=int, default=None,
              help="Region start in bp (default: 0)")
@click.option("--end", type=int, default=None,
              help="Region end in bp (default: full ARG sequence)")
@click.option("--no-self", is_flag=True,
              help="Omit (s, s) self-pair rows from output")
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def bgrm(ctx, run_id, start, end, no_self, output_path, fmt):
    """Branch GRM via ARG traversal (Fan, Mancuso & Chiang 2022).

    Validation baseline: matches the egrm.varGRM Python package element-wise
    to relative error < 1e-6. Requires an :ARGRun with the given RUN_ID
    already ingested via ``graphpop arg ingest``.
    """
    opts: dict[str, object] = {}
    if start is not None:
        opts["start"] = start
    if end is not None:
        opts["end"] = end
    if no_self:
        opts["include_self"] = False

    cypher = build_cypher(
        "graphpop.kinship.branch_grm",
        [f"'{run_id}'"],
        options=opts if opts else None,
        yield_cols=["sample_a", "sample_b", "phi", "ibs0", "het_het",
                    "n_snp", "n_aa_min", "method"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "kinship-bgrm",
                  {"run_id": run_id})


@kinship.command("bgrm-posterior")
@click.argument("run_ids", nargs=-1, required=True)
def bgrm_posterior(run_ids):
    """Posterior-aware branch GRM over multiple ARG runs (e.g. SINGER posterior).

    Not yet implemented — pending the M4.A ARG ingest and the
    branch_grm primary implementation.
    """
    raise click.ClickException(
        "graphpop.kinship.branch_grm_posterior is not yet implemented. "
        "Awaiting M4.A and branch_grm; see tasks/phase4_pairwise_plan.md."
    )
