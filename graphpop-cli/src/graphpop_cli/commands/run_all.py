"""graphpop run-all — orchestrate full-genome analysis across populations and chromosomes."""
from __future__ import annotations

import json
import time
from pathlib import Path

import click

from ..cli import pass_ctx
from ..config import build_cypher


# Default procedures for each phase
PHASE1_PROCEDURES = [
    "diversity", "sfs", "pop_summary", "roh", "ihs", "nsl", "garud_h",
]
PHASE2_PROCEDURES = ["xpehh", "divergence"]

YIELD_COLS = {
    "diversity": ["pi", "theta_w", "tajima_d", "fay_wu_h", "fay_wu_h_norm",
                   "het_exp", "het_obs", "fis", "n_variants", "n_segregating"],
    "sfs": ["sfs", "n_variants", "max_ac"],
    "pop_summary": ["pi", "theta_w", "tajima_d", "n_variants", "n_segregating"],
    "roh": ["sampleId", "n_roh", "total_length", "froh", "mean_length", "max_length"],
    "ihs": ["variantId", "pos", "af", "ihs_unstd", "ihs"],
    "nsl": ["variantId", "pos", "af", "nsl_unstd", "nsl"],
    "garud_h": ["chr", "start", "end", "population", "h1", "h12", "h2_h1",
                 "hap_diversity", "n_haplotypes", "n_variants"],
    "xpehh": ["variantId", "pos", "af_pop1", "af_pop2", "xpehh_unstd", "xpehh"],
    "divergence": ["fst_hudson", "fst_wc", "dxy", "da", "pbs", "n_variants"],
}


def get_chromosome_lengths(ctx) -> dict[str, int]:
    """Query chromosome lengths from the graph."""
    cypher = "MATCH (c:Chromosome) RETURN c.chromosomeId AS chr, c.length AS length ORDER BY c.chromosomeId"
    records = ctx.run(cypher)
    return {r["chr"]: r["length"] for r in records}


def get_populations(ctx) -> list[str]:
    """Query population IDs from the graph."""
    cypher = "MATCH (p:Population) WHERE p.n_samples > 1 RETURN p.populationId AS pop ORDER BY p.n_samples DESC"
    records = ctx.run(cypher)
    return [r["pop"] for r in records]


def get_xpehh_pairs(populations: list[str], max_pairs: int = 20) -> list[tuple[str, str]]:
    """Generate representative XP-EHH population pairs."""
    pairs = []
    for i, p1 in enumerate(populations):
        for p2 in populations[i + 1:]:
            pairs.append((p1, p2))
            if len(pairs) >= max_pairs:
                return pairs
    return pairs


@click.command("run-all")
@click.option("--phase", type=click.Choice(["1", "2", "all"]), default="all",
              help="Phase 1 (per-pop), Phase 2 (pairwise), or all")
@click.option("--output-dir", "-d", type=click.Path(), default="graphpop_results",
              help="Output directory for result files")
@click.option("--resume/--no-resume", default=True,
              help="Skip already-completed tasks (default: resume)")
@click.option("--json-output", type=click.Path(),
              help="Accumulated JSON results file (default: <output-dir>/results.json)")
@click.option("--persist/--no-persist", default=True,
              help="Write results to graph nodes (default: yes)")
@click.option("--populations", "pop_list",
              help="Comma-separated population list (default: auto-detect)")
@click.option("--chromosomes", "chr_list",
              help="Comma-separated chromosome list (default: auto-detect)")
@click.option("--xpehh-pairs", type=int, default=20,
              help="Max number of XP-EHH population pairs")
@click.option("--workers", type=int, default=1,
              help="Parallel workers (experimental)")
@pass_ctx
def run_all(ctx, phase, output_dir, resume, json_output, persist,
            pop_list, chr_list, xpehh_pairs, workers):
    """Run full-genome analysis across all populations and chromosomes.

    Phase 1: Per-population statistics (diversity, SFS, iHS, nSL, ROH, Garud's H)
    for each population × chromosome combination.

    Phase 2: Pairwise statistics (XP-EHH, divergence) for representative
    population pairs × chromosomes.

    Results are saved as TSV files in the output directory and optionally
    persisted to graph nodes.
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = Path(json_output) if json_output else out_dir / "results.json"

    # Load or initialize results
    results = {}
    if resume and json_path.exists():
        with open(json_path) as f:
            results = json.load(f)
        click.echo(f"Resuming from {json_path} ({len(results)} entries)")

    # Auto-detect populations and chromosomes
    click.echo("Querying graph for populations and chromosomes...")
    chr_lens = get_chromosome_lengths(ctx)
    chromosomes = sorted(chr_lens.keys()) if not chr_list else chr_list.split(",")
    populations = get_populations(ctx) if not pop_list else pop_list.split(",")

    click.echo(f"Populations: {len(populations)} ({', '.join(populations[:5])}...)")
    click.echo(f"Chromosomes: {len(chromosomes)} ({', '.join(chromosomes[:3])}...)")

    # Phase 1: Per-population
    if phase in ("1", "all"):
        click.echo(f"\n=== Phase 1: Per-population ({len(populations)} pops × {len(chromosomes)} chrs) ===")
        total = len(populations) * len(chromosomes) * len(PHASE1_PROCEDURES)
        done = 0
        t0 = time.time()

        for pop in populations:
            for chrom in chromosomes:
                for proc in PHASE1_PROCEDURES:
                    key = f"{pop}_{chrom}_{proc}"
                    if resume and key in results:
                        done += 1
                        continue

                    try:
                        if proc in ("diversity", "sfs", "pop_summary"):
                            length = chr_lens.get(chrom, 300_000_000)
                            cypher = build_cypher(
                                f"graphpop.{proc}",
                                [f"'{chrom}'", "1", str(length), f"'{pop}'"],
                                yield_cols=YIELD_COLS.get(proc),
                            )
                        elif proc in ("ihs", "nsl"):
                            cypher = build_cypher(
                                f"graphpop.{proc}",
                                [f"'{chrom}'", f"'{pop}'"],
                                options={"min_af": 0.05},
                                yield_cols=["variantId", "pos", "af",
                                            f"{proc}_unstd", proc],
                            )
                        elif proc == "roh":
                            cypher = build_cypher(
                                f"graphpop.{proc}",
                                [f"'{chrom}'", f"'{pop}'"],
                                yield_cols=YIELD_COLS["roh"],
                            )
                        elif proc == "garud_h":
                            cypher = build_cypher(
                                f"graphpop.{proc}",
                                [f"'{chrom}'", f"'{pop}'", "100000", "50000"],
                                yield_cols=YIELD_COLS["garud_h"],
                            )
                        else:
                            continue

                        records = ctx.run(cypher)
                        results[key] = {
                            "population": pop, "chr": chrom, "procedure": proc,
                            "n_records": len(records),
                            "summary": records[0] if len(records) == 1 else f"{len(records)} rows",
                        }

                        # Write per-procedure TSV
                        tsv_dir = out_dir / proc
                        tsv_dir.mkdir(exist_ok=True)
                        tsv_path = tsv_dir / f"{pop}_{chrom}.tsv"
                        if records:
                            _write_tsv(tsv_path, records)

                    except Exception as e:
                        results[key] = {"error": str(e)}
                        click.echo(f"  ERROR {key}: {e}", err=True)

                    done += 1
                    elapsed = time.time() - t0
                    rate = done / elapsed if elapsed > 0 else 0
                    eta = (total - done) / rate if rate > 0 else 0
                    if done % 10 == 0:
                        click.echo(
                            f"  [{done}/{total}] {key} "
                            f"({elapsed:.0f}s elapsed, ~{eta:.0f}s remaining)"
                        )

                # Save checkpoint after each chromosome
                _save_json(json_path, results)

    # Phase 2: Pairwise
    if phase in ("2", "all"):
        pairs = get_xpehh_pairs(populations, xpehh_pairs)
        click.echo(f"\n=== Phase 2: Pairwise ({len(pairs)} pairs × {len(chromosomes)} chrs) ===")

        for pop1, pop2 in pairs:
            for chrom in chromosomes:
                for proc in PHASE2_PROCEDURES:
                    key = f"{pop1}_vs_{pop2}_{chrom}_{proc}"
                    if resume and key in results:
                        continue

                    try:
                        if proc == "xpehh":
                            cypher = build_cypher(
                                "graphpop.xpehh",
                                [f"'{chrom}'", f"'{pop1}'", f"'{pop2}'"],
                                options={"min_af": 0.05},
                                yield_cols=YIELD_COLS["xpehh"],
                            )
                        elif proc == "divergence":
                            length = chr_lens.get(chrom, 300_000_000)
                            cypher = build_cypher(
                                "graphpop.divergence",
                                [f"'{chrom}'", "1", str(length),
                                 f"'{pop1}'", f"'{pop2}'"],
                                yield_cols=YIELD_COLS["divergence"],
                            )
                        else:
                            continue

                        records = ctx.run(cypher)
                        results[key] = {
                            "pop1": pop1, "pop2": pop2, "chr": chrom,
                            "procedure": proc, "n_records": len(records),
                            "summary": records[0] if len(records) == 1 else f"{len(records)} rows",
                        }

                        tsv_dir = out_dir / proc
                        tsv_dir.mkdir(exist_ok=True)
                        tsv_path = tsv_dir / f"{pop1}_vs_{pop2}_{chrom}.tsv"
                        if records:
                            _write_tsv(tsv_path, records)

                    except Exception as e:
                        results[key] = {"error": str(e)}
                        click.echo(f"  ERROR {key}: {e}", err=True)

                _save_json(json_path, results)

    # Final save
    _save_json(json_path, results)
    n_ok = sum(1 for v in results.values() if "error" not in v)
    n_err = sum(1 for v in results.values() if "error" in v)
    click.echo(f"\nDone. {n_ok} succeeded, {n_err} failed.")
    click.echo(f"Results: {json_path}")
    click.echo(f"TSV files: {out_dir}/")


def _write_tsv(path: Path, records: list[dict]):
    """Write records to a TSV file."""
    if not records:
        return
    keys = list(records[0].keys())
    with open(path, "w") as f:
        f.write("\t".join(keys) + "\n")
        for rec in records:
            vals = []
            for k in keys:
                v = rec[k]
                if isinstance(v, float):
                    vals.append(f"{v:.6g}")
                elif isinstance(v, list):
                    vals.append(",".join(str(x) for x in v))
                elif v is None:
                    vals.append("NA")
                else:
                    vals.append(str(v))
            f.write("\t".join(vals) + "\n")


def _save_json(path: Path, data: dict):
    """Save results as JSON with atomic write."""
    tmp = path.with_suffix(".tmp")
    with open(tmp, "w") as f:
        json.dump(data, f, indent=2, default=str)
    tmp.rename(path)
