"""graphpop report -- generate an automated HTML analysis summary report."""
from __future__ import annotations

import click

from ..cli import pass_ctx


@click.command()
@click.option("-o", "--output", "output_path", required=True,
              help="Output HTML file path")
@click.option("--database", help="Override database name for the report title")
@pass_ctx
def report(ctx, output_path, database):
    """Generate a self-contained HTML analysis report.

    Queries the graph database for dataset overview, per-population diversity,
    pairwise Fst, top selection signals, and annotation summary.

    \b
    Examples:
      graphpop report -o report.html
      graphpop report --database rice3k -o rice3k_report.html
    """
    db_name = database or ctx.database or "GraphPop"
    click.echo(f"Generating report for database: {db_name} ...", err=True)

    # ---- 1. Dataset overview ------------------------------------------------
    overview_rows = ctx.run(
        "MATCH (v:Variant) "
        "WITH count(v) AS n_variants "
        "OPTIONAL MATCH (s:Sample) "
        "WITH n_variants, count(DISTINCT s) AS n_samples "
        "OPTIONAL MATCH (g:Gene) "
        "WITH n_variants, n_samples, count(DISTINCT g) AS n_genes "
        "RETURN n_variants, n_samples, n_genes"
    )
    overview = overview_rows[0] if overview_rows else {
        "n_variants": 0, "n_samples": 0, "n_genes": 0,
    }

    pop_list_rows = ctx.run(
        "MATCH (v:Variant) WHERE v.pop_ids IS NOT NULL "
        "RETURN v.pop_ids AS pids LIMIT 1"
    )
    populations = pop_list_rows[0]["pids"] if pop_list_rows else []

    chr_rows = ctx.run(
        "MATCH (v:Variant) "
        "RETURN DISTINCT v.chr AS chr ORDER BY chr"
    )
    chromosomes = [r["chr"] for r in chr_rows]

    # ---- 2. Per-population diversity ----------------------------------------
    diversity_rows = ctx.run(
        "MATCH (w:GenomicWindow) "
        "WHERE w.pi IS NOT NULL "
        "RETURN w.population AS population, "
        "       avg(w.pi) AS mean_pi, "
        "       avg(w.theta_w) AS mean_theta_w, "
        "       avg(w.tajima_d) AS mean_tajima_d "
        "ORDER BY mean_pi DESC"
    )

    # ---- 3. Pairwise Fst (top 10 pairs) ------------------------------------
    fst_rows = ctx.run(
        "MATCH (w:GenomicWindow) "
        "WHERE w.fst IS NOT NULL AND w.pop_pair IS NOT NULL "
        "RETURN w.pop_pair AS pop_pair, avg(w.fst) AS mean_fst "
        "ORDER BY mean_fst DESC LIMIT 10"
    )

    # ---- 4. Top selection signals -------------------------------------------
    ihs_rows = ctx.run(
        "MATCH (v:Variant) "
        "WHERE any(k IN keys(v) WHERE k STARTS WITH 'ihs_') "
        "WITH v, [k IN keys(v) WHERE k STARTS WITH 'ihs_'] AS ks "
        "UNWIND ks AS k "
        "WITH v.variantId AS variant, v.chr AS chr, v.pos AS pos, "
        "     k AS stat, v[k] AS value "
        "WHERE abs(value) > 2.0 "
        "RETURN variant, chr, pos, stat, value "
        "ORDER BY abs(value) DESC LIMIT 20"
    )

    xpehh_rows = ctx.run(
        "MATCH (v:Variant) "
        "WHERE any(k IN keys(v) WHERE k STARTS WITH 'xpehh_') "
        "WITH v, [k IN keys(v) WHERE k STARTS WITH 'xpehh_'] AS ks "
        "UNWIND ks AS k "
        "WITH v.variantId AS variant, v.chr AS chr, v.pos AS pos, "
        "     k AS stat, v[k] AS value "
        "WHERE abs(value) > 2.0 "
        "RETURN variant, chr, pos, stat, value "
        "ORDER BY abs(value) DESC LIMIT 20"
    )

    selection_rows = ihs_rows + xpehh_rows
    selection_rows.sort(key=lambda r: abs(r.get("value", 0)), reverse=True)
    selection_rows = selection_rows[:20]

    # ---- 5. Annotation summary ----------------------------------------------
    annot_rows = ctx.run(
        "OPTIONAL MATCH (g:Gene) "
        "WITH count(DISTINCT g) AS n_genes "
        "OPTIONAL MATCH (pw:Pathway) "
        "WITH n_genes, count(DISTINCT pw) AS n_pathways "
        "OPTIONAL MATCH (go:GOTerm) "
        "RETURN n_genes, n_pathways, count(DISTINCT go) AS n_go_terms"
    )
    annot = annot_rows[0] if annot_rows else {
        "n_genes": 0, "n_pathways": 0, "n_go_terms": 0,
    }

    # ---- Build HTML ---------------------------------------------------------
    def _table(headers, rows_data):
        """Build an HTML table from headers and list-of-lists."""
        lines = ['<table>', '<tr>' + ''.join(f'<th>{h}</th>' for h in headers) + '</tr>']
        for row in rows_data:
            lines.append('<tr>' + ''.join(f'<td>{_fmt(v)}</td>' for v in row) + '</tr>')
        lines.append('</table>')
        return '\n'.join(lines)

    def _fmt(v):
        if v is None:
            return "NA"
        if isinstance(v, float):
            return f"{v:.6g}"
        if isinstance(v, list):
            return ", ".join(str(x) for x in v)
        return str(v)

    # Overview table
    overview_html = _table(
        ["Metric", "Value"],
        [
            ["Variants", overview.get("n_variants", 0)],
            ["Samples", overview.get("n_samples", 0)],
            ["Genes", overview.get("n_genes", 0)],
            ["Populations", len(populations)],
            ["Population IDs", ", ".join(str(p) for p in populations)],
            ["Chromosomes", ", ".join(str(c) for c in chromosomes)],
        ],
    )

    # Diversity table
    if diversity_rows:
        div_html = _table(
            ["Population", "Mean pi", "Mean theta_W", "Mean Tajima's D"],
            [[r["population"], r["mean_pi"], r["mean_theta_w"], r["mean_tajima_d"]]
             for r in diversity_rows],
        )
    else:
        div_html = "<p>No GenomicWindow diversity data found.</p>"

    # Fst table
    if fst_rows:
        fst_html = _table(
            ["Population Pair", "Mean Fst"],
            [[r["pop_pair"], r["mean_fst"]] for r in fst_rows],
        )
    else:
        fst_html = "<p>No pairwise Fst data found in GenomicWindow nodes.</p>"

    # Selection signals table
    if selection_rows:
        sel_html = _table(
            ["Variant", "Chr", "Pos", "Statistic", "Value"],
            [[r["variant"], r["chr"], r["pos"], r["stat"], r["value"]]
             for r in selection_rows],
        )
    else:
        sel_html = "<p>No iHS or XP-EHH signals above threshold found.</p>"

    # Annotation summary table
    annot_html = _table(
        ["Annotation Type", "Count"],
        [
            ["Genes", annot.get("n_genes", 0)],
            ["Pathways", annot.get("n_pathways", 0)],
            ["GO Terms", annot.get("n_go_terms", 0)],
        ],
    )

    html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>GraphPop Report: {db_name}</title>
<style>
body {{
    font-family: Arial, Helvetica, sans-serif;
    max-width: 960px;
    margin: 40px auto;
    padding: 0 20px;
    color: #333;
    line-height: 1.5;
}}
h1 {{
    color: #0072B2;
    border-bottom: 2px solid #0072B2;
    padding-bottom: 8px;
}}
h2 {{
    color: #555;
    margin-top: 32px;
}}
table {{
    border-collapse: collapse;
    width: 100%;
    margin: 12px 0 24px 0;
    font-size: 13px;
}}
th, td {{
    border: 1px solid #ddd;
    padding: 6px 10px;
    text-align: left;
}}
th {{
    background-color: #0072B2;
    color: white;
    font-weight: 600;
}}
tr:nth-child(even) {{
    background-color: #f9f9f9;
}}
tr:hover {{
    background-color: #e9f3fb;
}}
footer {{
    margin-top: 48px;
    padding-top: 12px;
    border-top: 1px solid #ddd;
    font-size: 11px;
    color: #999;
}}
</style>
</head>
<body>
<h1>GraphPop Analysis Report: {db_name}</h1>

<h2>Dataset Overview</h2>
{overview_html}

<h2>Population Diversity</h2>
{div_html}

<h2>Pairwise Fst (Top 10)</h2>
{fst_html}

<h2>Top Selection Signals</h2>
{sel_html}

<h2>Annotation Summary</h2>
{annot_html}

<footer>Generated by GraphPop CLI v0.1.0</footer>
</body>
</html>
"""

    with open(output_path, "w") as f:
        f.write(html)
    click.echo(f"Report saved to: {output_path}")
