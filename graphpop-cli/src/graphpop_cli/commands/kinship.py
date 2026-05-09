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


@kinship.command("ibs")
@click.argument("chr")
@click.argument("population")
@click.option("--start", type=int, default=1, show_default=True)
@click.option("--end", type=int, default=None,
              help="End position (default: end of chromosome)")
@click.option("--min-snp", type=int, default=1000, show_default=True)
@click.option("--min-ibs", type=float, default=0.0, show_default=True)
@click.option("--include-self", is_flag=True)
@click.option("--samples", multiple=True,
              help="Restrict to a subset of sample IDs")
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def ibs(ctx, chr, population, start, end, min_snp, min_ibs, include_self,
        samples, output_path, fmt):
    """Identity-by-state pairwise statistic on packed genotypes (M4.2).

    Matrix-only sibling of `kinship king`; no ARG required. IBS is in
    [0, 1]; identical samples = 1, opposite homozygotes average toward 0.
    """
    opts: dict[str, object] = {
        "min_snp": min_snp,
        "min_ibs": min_ibs,
        "include_self": include_self,
        "start": start,
    }
    if end is not None:
        opts["end"] = end
    if samples:
        opts["samples"] = list(samples)

    cypher = build_cypher(
        "graphpop.kinship.ibs",
        [f"'{chr}'", f"'{population}'"],
        options=opts,
        yield_cols=["sample_a", "sample_b", "phi", "ibs0", "het_het",
                    "n_snp", "n_aa_min", "method"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "kinship-ibs",
                  {"chr": chr, "pop": population})


@kinship.command("bgrm")
@click.argument("run_id")
@click.option("--start", type=int, default=None,
              help="Region start in bp (default: 0)")
@click.option("--end", type=int, default=None,
              help="Region end in bp (default: full ARG sequence)")
@click.option("--no-self", is_flag=True,
              help="Omit (s, s) self-pair rows from output")
@click.option("--pathway",
              help="Restrict to branches whose mutations land in this Pathway "
                   "(closes G2: annotation co-residence)")
@click.option("--consequence",
              help="Restrict to branches whose mutations have this consequence "
                   "(e.g. 'missense_variant'; closes G2)")
@click.option("--time-window-start", type=float,
              help="Lower bound (generations) of the time window")
@click.option("--time-window-end", type=float,
              help="Upper bound (generations) of the time window. Both "
                   "--time-window-start and --time-window-end must be set.")
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def bgrm(ctx, run_id, start, end, no_self, pathway, consequence,
         time_window_start, time_window_end, output_path, fmt):
    """Branch GRM via ARG traversal (Fan, Mancuso & Chiang 2022).

    Unconditional mode validates against egrm.varGRM (rel err < 1e-6).
    Conditional predicates (pathway / consequence / time-window) are the
    GraphPop-novel statistics that close gaps G1 (time-stratification)
    and G2 (annotation co-residence) in the ARG literature.

    Examples:

      graphpop kinship bgrm tsinfer_chr22_v1
      graphpop kinship bgrm tsinfer_chr22_v1 --pathway GO:0006281
      graphpop kinship bgrm tsinfer_chr22_v1 --consequence missense_variant
      graphpop kinship bgrm tsinfer_chr22_v1 --time-window-start 0 --time-window-end 1000
      graphpop kinship bgrm tsinfer_chr22_v1 --pathway P_test --time-window-start 0 --time-window-end 0.5

    Requires an :ARGRun with the given RUN_ID already ingested via
    ``graphpop arg ingest``.
    """
    opts: dict[str, object] = {}
    if start is not None:
        opts["start"] = start
    if end is not None:
        opts["end"] = end
    if no_self:
        opts["include_self"] = False
    if pathway:
        opts["restrict_to_pathway"] = pathway
    if consequence:
        opts["mutation_filter"] = consequence
    if (time_window_start is None) ^ (time_window_end is None):
        raise click.ClickException(
            "--time-window-start and --time-window-end must be set together"
        )
    if time_window_start is not None and time_window_end is not None:
        if time_window_end <= time_window_start:
            raise click.ClickException(
                "--time-window-end must be greater than --time-window-start"
            )
        opts["time_window"] = [time_window_start, time_window_end]

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


@kinship.command("bgrm-apply")
@click.argument("run_id")
@click.argument("vector_tsv", type=click.Path(exists=True, dir_okay=False))
@click.option("--start", type=int, default=None)
@click.option("--end", type=int, default=None)
@click.option("--pathway",
              help="Restrict to branches whose mutations land in this Pathway")
@click.option("--consequence",
              help="Restrict to branches whose mutations have this consequence")
@click.option("--time-window-start", type=float)
@click.option("--time-window-end", type=float)
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def bgrm_apply(ctx, run_id, vector_tsv, start, end,
                pathway, consequence,
                time_window_start, time_window_end,
                output_path, fmt):
    """Apply the branch GRM to a vector via Algorithm V (G * v).

    Vector TSV: one float per line, length = number of samples in the
    run (haplotype-ordered). Conditional predicates compose unchanged.

    Tractable at biobank scale; no full matrix is materialised.
    """
    v: list[float] = []
    with open(vector_tsv) as fh:
        for lineno, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            try:
                v.append(float(line))
            except ValueError as exc:
                raise click.ClickException(
                    f"line {lineno}: not a float: {raw!r}") from exc

    opts: dict[str, object] = {}
    if start is not None:
        opts["start"] = start
    if end is not None:
        opts["end"] = end
    if pathway:
        opts["restrict_to_pathway"] = pathway
    if consequence:
        opts["mutation_filter"] = consequence
    if (time_window_start is None) ^ (time_window_end is None):
        raise click.ClickException(
            "--time-window-start and --time-window-end must be set together")
    if time_window_start is not None and time_window_end is not None:
        opts["time_window"] = [time_window_start, time_window_end]

    cypher = ("CALL graphpop.kinship.branch_grm_apply($run_id, $v, $opts) "
              "YIELD sample_id, col, value, n_branches, method")
    records = ctx.run(cypher, {"run_id": run_id, "v": v, "opts": opts})
    format_output(records, output_path, fmt, "kinship-bgrm-apply",
                  {"run_id": run_id, "vector_length": len(v)})


@kinship.command("bgrm-by-ancestry")
@click.argument("run_id")
@click.option("--start", type=int, default=None,
              help="Region start in bp")
@click.option("--end", type=int, default=None,
              help="Region end in bp")
@click.option("--no-self", is_flag=True,
              help="Omit (s, s) self-pair rows")
@click.option("--painter", default=None,
              help="Painter id (e.g. 'majority_vote', 'rfmix') when "
                   "the run has multiple paintings")
@click.option("--pathway",
              help="Restrict to branches whose mutations land in this Pathway")
@click.option("--consequence",
              help="Restrict to branches whose mutations have this consequence")
@click.option("--time-window-start", type=float)
@click.option("--time-window-end", type=float)
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def bgrm_by_ancestry(ctx, run_id, start, end, no_self, painter,
                      pathway, consequence,
                      time_window_start, time_window_end,
                      output_path, fmt):
    """Branch GRM decomposed by local ancestry painting (closes G3).

    Returns one row per (sample_a, sample_b, ancestry) triple. Sum
    across ancestries reproduces the unconditional branch_grm matrix
    (modulo unpainted nodes). Conditional predicates (pathway,
    consequence, time_window) compose unchanged.

    Requires :HAS_ANCESTRY edges ingested via ``graphpop ancestry
    ingest`` or ``graphpop ancestry ingest-from-samples``.
    """
    opts: dict[str, object] = {}
    if start is not None:
        opts["start"] = start
    if end is not None:
        opts["end"] = end
    if no_self:
        opts["include_self"] = False
    if painter:
        opts["painter"] = painter
    if pathway:
        opts["restrict_to_pathway"] = pathway
    if consequence:
        opts["mutation_filter"] = consequence
    if (time_window_start is None) ^ (time_window_end is None):
        raise click.ClickException(
            "--time-window-start and --time-window-end must be set together")
    if time_window_start is not None and time_window_end is not None:
        if time_window_end <= time_window_start:
            raise click.ClickException(
                "--time-window-end must be greater than --time-window-start")
        opts["time_window"] = [time_window_start, time_window_end]

    cypher = build_cypher(
        "graphpop.kinship.branch_grm_by_ancestry",
        [f"'{run_id}'"],
        options=opts if opts else None,
        yield_cols=["sample_a", "sample_b", "ancestry",
                    "b_ij_component", "n_branches", "method"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "kinship-bgrm-by-ancestry",
                  {"run_id": run_id})


@kinship.command("bgrm-posterior")
@click.argument("run_ids", nargs=-1, required=True)
@click.option("--start", type=int, default=None,
              help="Region start in bp (same for every run)")
@click.option("--end", type=int, default=None,
              help="Region end in bp (same for every run)")
@click.option("--no-self", is_flag=True,
              help="Omit (s, s) self-pair rows from output")
@click.option("--pathway",
              help="Restrict to branches whose mutations land in this Pathway")
@click.option("--consequence",
              help="Restrict to branches whose mutations have this consequence")
@click.option("--time-window-start", type=float,
              help="Lower bound (generations) of the time window")
@click.option("--time-window-end", type=float,
              help="Upper bound (generations) of the time window")
@click.option("-o", "--output", "output_path",
              help="Output file (default: stdout)")
@click.option("--format", "fmt", default="tsv",
              type=click.Choice(["tsv", "csv", "json"]))
@pass_ctx
def bgrm_posterior(ctx, run_ids, start, end, no_self, pathway, consequence,
                    time_window_start, time_window_end, output_path, fmt):
    """Posterior-aware branch GRM over multiple ARG runs (e.g. SINGER posterior).

    Aggregates per-run branch_grm with element-wise Welford. Returns
    posterior mean and standard error of the mean. Conditional
    predicates (pathway / consequence / time-window) compose unchanged
    — the posterior is computed *after* the predicate.

    Closes the ARG-inference uncertainty (G1) gap. Combine with
    --pathway to get "posterior-mean kinship through pathway X with
    credible interval" — the GraphPop v2 paper's headline statistic.

    Run IDs may be supplied either positionally or comma-separated:

      graphpop kinship bgrm-posterior r0 r1 r2 r3 r4
      graphpop kinship bgrm-posterior r0,r1,r2,r3,r4
    """
    expanded: list[str] = []
    for rid in run_ids:
        expanded.extend(p for p in rid.split(",") if p)
    if not expanded:
        raise click.ClickException("at least one run-id required")

    opts: dict[str, object] = {}
    if start is not None:
        opts["start"] = start
    if end is not None:
        opts["end"] = end
    if no_self:
        opts["include_self"] = False
    if pathway:
        opts["restrict_to_pathway"] = pathway
    if consequence:
        opts["mutation_filter"] = consequence
    if (time_window_start is None) ^ (time_window_end is None):
        raise click.ClickException(
            "--time-window-start and --time-window-end must be set together")
    if time_window_start is not None and time_window_end is not None:
        if time_window_end <= time_window_start:
            raise click.ClickException(
                "--time-window-end must be greater than --time-window-start")
        opts["time_window"] = [time_window_start, time_window_end]

    rids_literal = "[" + ",".join(f"'{r}'" for r in expanded) + "]"
    cypher = build_cypher(
        "graphpop.kinship.branch_grm_posterior",
        [rids_literal],
        options=opts if opts else None,
        yield_cols=["sample_a", "sample_b", "b_ij_mean", "b_ij_sd",
                    "n_runs", "method"],
    )
    records = ctx.run(cypher)
    format_output(records, output_path, fmt, "kinship-bgrm-posterior",
                  {"n_runs": len(expanded)})
