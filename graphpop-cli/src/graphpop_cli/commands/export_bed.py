"""graphpop export-bed — export high-scoring regions as BED format."""
from __future__ import annotations

import sys

import click

from ..cli import pass_ctx


# Statistics stored on GenomicWindow nodes vs Variant nodes
WINDOW_STATS = {"fst", "pi", "tajima_d", "h12"}
VARIANT_STATS = {"ihs", "xpehh", "nsl"}


@click.command("export-bed")
@click.option("--stat", required=True,
              type=click.Choice(["fst", "pi", "tajima_d", "h12", "ihs", "xpehh", "nsl"]),
              help="Statistic to threshold")
@click.option("--threshold", required=True, type=float,
              help="Minimum value to include (absolute value for ihs/xpehh/nsl)")
@click.option("--pop", "population", required=True, help="Population name")
@click.option("--pop2", help="Second population (required for xpehh)")
@click.option("--chr", "chromosome", help="Chromosome filter (optional)")
@click.option("--merge-distance", type=int, default=100000,
              help="Merge variants within this distance into intervals (variant-based stats, default: 100000)")
@click.option("-o", "--output", "output_path", required=True, help="Output BED file")
@pass_ctx
def export_bed(ctx, stat, threshold, population, pop2, chromosome,
               merge_distance, output_path):
    """Export regions exceeding a statistic threshold as BED format.

    For window-based stats (fst, pi, tajima_d, h12), queries GenomicWindow
    nodes directly. For variant-based stats (ihs, xpehh, nsl), merges
    consecutive high-scoring variants into intervals using --merge-distance.

    Output: standard 5-column BED (chr, start, end, name, score).

    \b
    Examples:
      graphpop export-bed --stat fst --threshold 0.5 --pop EUR -o high_fst.bed
      graphpop export-bed --stat ihs --threshold 2.5 --pop EUR --chr chr22 -o ihs_peaks.bed
      graphpop export-bed --stat xpehh --threshold 3.0 --pop EUR --pop2 AFR -o xpehh.bed
      graphpop export-bed --stat tajima_d --threshold -2.0 --pop GJ-tmp -o tajimad.bed
    """
    if stat == "xpehh" and not pop2:
        click.echo("Error: --pop2 is required for xpehh.", err=True)
        raise SystemExit(1)

    bed_name = f"{stat}_{population}" if not pop2 else f"{stat}_{population}_{pop2}"

    if stat in WINDOW_STATS:
        records = _query_window_stat(ctx, stat, threshold, population, pop2, chromosome)
        bed_lines = _windows_to_bed(records, stat, bed_name)
    else:
        records = _query_variant_stat(ctx, stat, threshold, population, pop2, chromosome)
        bed_lines = _merge_variants_to_bed(records, merge_distance, bed_name)

    if not bed_lines:
        click.echo(f"No regions found exceeding threshold {threshold} for {stat}.", err=True)
        return

    with open(output_path, "w") as f:
        for line in bed_lines:
            f.write(line + "\n")

    click.echo(f"Wrote {len(bed_lines)} BED intervals to {output_path}.", err=True)


def _query_window_stat(ctx, stat, threshold, population, pop2, chromosome):
    """Query GenomicWindow nodes for window-based statistics."""
    where_parts = [f"w.population = '{population}'"]
    if chromosome:
        where_parts.append(f"w.chr = '{chromosome}'")

    prop = stat
    if stat == "fst" and pop2:
        prop = f"fst_{population}_{pop2}"

    # For Tajima's D, extreme negative values indicate selection
    if stat == "tajima_d":
        where_parts.append(f"w.{prop} <= {threshold}")
    else:
        where_parts.append(f"w.{prop} >= {threshold}")

    cypher = (
        f"MATCH (w:GenomicWindow) "
        f"WHERE {' AND '.join(where_parts)} "
        f"RETURN w.chr AS chr, w.start AS start, w.end AS end, "
        f"w.{prop} AS score "
        f"ORDER BY w.chr, w.start"
    )
    return ctx.run(cypher)


def _query_variant_stat(ctx, stat, threshold, population, pop2, chromosome):
    """Query Variant nodes for variant-based statistics."""
    if stat == "xpehh":
        prop = f"xpehh_{population}_{pop2}"
    else:
        prop = f"{stat}_{population}"

    where_parts = [f"v.{prop} IS NOT NULL", f"abs(v.{prop}) >= {threshold}"]
    if chromosome:
        where_parts.append(f"v.chr = '{chromosome}'")

    cypher = (
        f"MATCH (v:Variant) "
        f"WHERE {' AND '.join(where_parts)} "
        f"RETURN v.chr AS chr, v.pos AS pos, v.{prop} AS score "
        f"ORDER BY v.chr, v.pos"
    )
    return ctx.run(cypher)


def _windows_to_bed(records, stat, bed_name):
    """Convert window records directly to BED lines."""
    lines = []
    for r in records:
        chrom = r.get("chr", "")
        start = r.get("start", 0)
        end = r.get("end", 0)
        score = r.get("score", 0)
        score_str = f"{score:.6g}" if isinstance(score, float) else str(score)
        lines.append(f"{chrom}\t{start}\t{end}\t{bed_name}\t{score_str}")
    return lines


def _merge_variants_to_bed(records, merge_distance, bed_name):
    """Merge consecutive high-scoring variants into BED intervals."""
    if not records:
        return []

    intervals = []
    current_chr = None
    current_start = None
    current_end = None
    current_scores = []

    for r in records:
        chrom = r.get("chr", "")
        pos = r.get("pos", 0)
        score = r.get("score", 0)

        if (current_chr is None or chrom != current_chr
                or pos - current_end > merge_distance):
            # Emit previous interval
            if current_chr is not None:
                mean_score = sum(current_scores) / len(current_scores)
                intervals.append((current_chr, current_start, current_end, mean_score))
            # Start new interval
            current_chr = chrom
            current_start = pos
            current_end = pos
            current_scores = [abs(score) if isinstance(score, (int, float)) else 0]
        else:
            current_end = pos
            current_scores.append(abs(score) if isinstance(score, (int, float)) else 0)

    # Emit last interval
    if current_chr is not None:
        mean_score = sum(current_scores) / len(current_scores)
        intervals.append((current_chr, current_start, current_end, mean_score))

    lines = []
    for chrom, start, end, score in intervals:
        lines.append(f"{chrom}\t{start}\t{end}\t{bed_name}\t{score:.6g}")
    return lines
