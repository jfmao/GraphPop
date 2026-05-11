"""Paper 2 Fig 5 figure-gen — graph-native query plane (G5).

Two single-column vector PDFs:

- **fig5a.pdf** — the composed Cypher query rendered as a code
  block, with each predicate annotated to the right.
- **fig5d.pdf** — ecosystem map: a 4-quadrant box layout
  showing Java procedures, competitor wrappers, paper-2
  drivers, and the CLI surface.

No msprime / wrapper dep — reads only the JSON+text artefacts
produced by `fig5_panels.py`.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np


_RCPARAMS = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.titlesize": 7,
    "axes.labelsize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}

SINGLE_COL = (3.5, 4.5)         # slightly taller for code panel
ONE_HALF_COL = (5.0, 3.75)      # ecosystem map fits 1.5-col

WONG = {
    "blue":       "#0072B2",
    "orange":     "#E69F00",
    "green":      "#009E73",
    "vermillion": "#D55E00",
    "purple":     "#CC79A7",
    "sky":        "#56B4E9",
    "grey":       "#999999",
}


# ---------------------------------------------------------------------------
# Fig 5a — Cypher code block
# ---------------------------------------------------------------------------

def make_fig5a(
    template_path: Path,
    annotations_path: Path,
    out_pdf: Path,
) -> dict:
    """Render the Cypher query as a code panel."""
    template = template_path.read_text().rstrip("\n")
    annotations = json.loads(annotations_path.read_text())

    with plt.rc_context(_RCPARAMS):
        fig, ax = plt.subplots(figsize=SINGLE_COL)
        ax.set_axis_off()
        # Header.
        ax.text(
            0.02, 0.98,
            "Composed predicate query — Cypher",
            transform=ax.transAxes,
            fontsize=8, fontweight="bold",
            va="top", ha="left",
        )
        # Code block (monospace).
        ax.text(
            0.02, 0.92, template,
            transform=ax.transAxes,
            family="monospace",
            fontsize=5.0, va="top", ha="left",
            linespacing=1.25,
            bbox=dict(
                facecolor="#f7f7f7", edgecolor="#cccccc",
                linewidth=0.4, boxstyle="round,pad=0.4",
            ),
        )
        # Annotations beneath.
        annot_lines = [
            f"• {a['token']} — {a['explanation']}"
            for a in annotations
        ]
        ax.text(
            0.02, 0.20, "Composition annotations:\n" + "\n".join(annot_lines),
            transform=ax.transAxes,
            fontsize=5.5, va="top", ha="left",
        )
        ax.text(
            -0.05, 1.02, "a",
            transform=ax.transAxes,
            fontsize=9, fontweight="bold",
            va="top", ha="left",
        )
        out_pdf.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_pdf, format="pdf")
        plt.close(fig)

    return {
        "n_lines": len(template.splitlines()),
        "n_annotations": len(annotations),
        "out_pdf": str(out_pdf),
    }


# ---------------------------------------------------------------------------
# Fig 5d — Ecosystem map
# ---------------------------------------------------------------------------

def _draw_quadrant(
    ax, x, y, w, h, title, items, color,
) -> None:
    """Draw a single quadrant rectangle with a title + items list."""
    rect = mpatches.FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.01,rounding_size=0.012",
        linewidth=0.6, facecolor=color, alpha=0.10,
        edgecolor=color, transform=ax.transAxes,
    )
    ax.add_patch(rect)
    ax.text(
        x + 0.01, y + h - 0.025, title,
        transform=ax.transAxes,
        fontsize=7, fontweight="bold",
        color=color, va="top", ha="left",
    )
    body = "\n".join(items)
    ax.text(
        x + 0.012, y + h - 0.06, body,
        transform=ax.transAxes,
        family="monospace",
        fontsize=4.5, va="top", ha="left",
        linespacing=1.30,
    )


def make_fig5d(
    eco_json: Path, out_pdf: Path,
    *, max_procedures: int = 16,
) -> dict:
    """Render the ecosystem JSON as a 4-quadrant map."""
    eco = json.loads(eco_json.read_text())

    # Java procedures: group by module, present module: name list.
    proc_by_module: Dict[str, List[str]] = {}
    for p in eco["java_procedures"]:
        proc_by_module.setdefault(p["module"], []).append(p["name"])
    java_lines: List[str] = []
    total_procs = 0
    for module in sorted(proc_by_module):
        names = sorted(proc_by_module[module])
        java_lines.append(f"[{module}]")
        for n in names:
            java_lines.append(f"  {n}")
            total_procs += 1
            if total_procs >= max_procedures:
                java_lines.append(f"  …(+{len(eco['java_procedures']) - max_procedures} more)")
                break
        if total_procs >= max_procedures:
            break

    wrapper_lines = [
        f"{w['label']} → {w['dependency']}"
        for w in eco["competitor_wrappers"]
    ]
    driver_lines = [
        f"{d['figure']}: {d['module'].split('.')[-1]}"
        for d in eco["paper2_drivers"]
    ]
    cli_lines = eco["cli_commands"]

    with plt.rc_context(_RCPARAMS):
        fig, ax = plt.subplots(figsize=ONE_HALF_COL)
        ax.set_axis_off()
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)

        # Header strip.
        ax.text(
            0.5, 0.98,
            "GraphPop ecosystem — Paper 2 query plane",
            transform=ax.transAxes,
            fontsize=8, fontweight="bold",
            va="top", ha="center",
        )
        ax.text(
            -0.03, 1.02, "d",
            transform=ax.transAxes,
            fontsize=9, fontweight="bold",
            va="top", ha="left",
        )

        # 4 quadrants.
        _draw_quadrant(
            ax, 0.02, 0.50, 0.46, 0.42,
            f"Java procedures (Cypher, n={len(eco['java_procedures'])})",
            java_lines, WONG["blue"])
        _draw_quadrant(
            ax, 0.52, 0.50, 0.46, 0.42,
            f"Competitor wrappers (n={len(eco['competitor_wrappers'])})",
            wrapper_lines, WONG["green"])
        _draw_quadrant(
            ax, 0.02, 0.05, 0.46, 0.42,
            f"Paper-2 drivers (n={len(eco['paper2_drivers'])})",
            driver_lines, WONG["orange"])
        _draw_quadrant(
            ax, 0.52, 0.05, 0.46, 0.42,
            f"CLI / module entries (n={len(eco['cli_commands'])})",
            cli_lines, WONG["purple"])

        out_pdf.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_pdf, format="pdf")
        plt.close(fig)

    return {
        "n_java_procedures": len(eco["java_procedures"]),
        "n_competitor_wrappers": len(eco["competitor_wrappers"]),
        "n_paper2_drivers": len(eco["paper2_drivers"]),
        "n_cli_commands": len(eco["cli_commands"]),
        "out_pdf": str(out_pdf),
    }


# ---------------------------------------------------------------------------
# Fig 5b — cryptic-pair recall vs injected IBD fraction
# ---------------------------------------------------------------------------

def make_fig5b(
    panel_csv: Path, out_pdf: Path,
    *, success_threshold: float = 0.8,
) -> dict:
    """Recall vs. injected IBD fraction (per-replicate scatter +
    mean curve). Horizontal at success_threshold; legend lists the
    phi threshold used."""
    import csv as _csv
    import numpy as _np

    by_fraction: dict = {}
    for row in _csv.DictReader(open(panel_csv)):
        f = float(row["fraction"])
        by_fraction.setdefault(f, []).append(float(row["recall"]))

    if not by_fraction:
        raise ValueError(f"no rows in {panel_csv}")

    fractions = sorted(by_fraction.keys())
    means = [_np.mean(by_fraction[f]) for f in fractions]
    stds = [_np.std(by_fraction[f], ddof=1) if len(by_fraction[f]) > 1
            else 0.0 for f in fractions]

    with plt.rc_context(_RCPARAMS):
        fig, ax = plt.subplots(figsize=SINGLE_COL)
        # Per-replicate scatter (jitter on x for visibility).
        rng = _np.random.default_rng(0)
        for f, recalls in by_fraction.items():
            xs = _np.full(len(recalls), f) * (
                1 + rng.normal(scale=0.02, size=len(recalls)))
            ax.scatter(
                xs, recalls, s=8, alpha=0.55,
                color=WONG["sky"], edgecolor="black",
                linewidth=0.25, zorder=2,
            )
        # Mean line.
        ax.errorbar(
            fractions, means, yerr=stds,
            marker="o", linestyle="-", color=WONG["blue"],
            markeredgecolor="black", markeredgewidth=0.3,
            ecolor=WONG["blue"], elinewidth=0.6, capsize=2,
            label="mean ± std", zorder=3,
        )
        ax.axhline(
            success_threshold, color=WONG["green"],
            linewidth=0.5, linestyle="--",
            label=f"benchmark target ({success_threshold:.2f})",
            zorder=1,
        )
        ax.set_xscale("log")
        ax.set_xlabel("injected IBD fraction f")
        ax.set_ylabel("recall at φ > 0.0625")
        ax.set_ylim(-0.05, 1.05)
        ax.set_title("Cryptic-pair recall (simulation proxy)")
        ax.legend(frameon=False, loc="lower right", fontsize=5.5)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.text(
            -0.18, 1.05, "b",
            transform=ax.transAxes,
            fontsize=9, fontweight="bold",
            va="top", ha="left",
        )
        out_pdf.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_pdf, format="pdf")
        plt.close(fig)
    return {
        "fractions": list(fractions),
        "mean_recall": list(means),
        "out_pdf": str(out_pdf),
    }


# ---------------------------------------------------------------------------
# Fig 5c — Cypher vs Pipeline LOC comparison
# ---------------------------------------------------------------------------

def make_fig5c(
    panel_csv: Path,
    out_pdf: Path,
) -> dict:
    """Grouped bar chart: pipeline-vs-Cypher on (files, code-lines, runtime)."""
    import csv as _csv

    pipeline = {"files": 0, "code_lines": 0, "runtime_s": 0.0}
    graphpop = {"files": 0, "code_lines": 0, "runtime_s": 0.0}
    with open(panel_csv) as fh:
        for row in _csv.DictReader(fh):
            target = pipeline if row["source"] == "pipeline" else graphpop
            target["files"] += 1
            target["code_lines"] += int(row["code_lines"])
            target["runtime_s"] += float(row["runtime_s"])

    metrics = ["files", "code_lines", "runtime_s"]
    labels = ["files", "code lines", "runtime (s)"]

    with plt.rc_context(_RCPARAMS):
        fig, ax = plt.subplots(figsize=SINGLE_COL)
        x = np.arange(len(metrics))
        bar_w = 0.36
        pipeline_vals = [pipeline[m] for m in metrics]
        graphpop_vals = [graphpop[m] for m in metrics]

        ax.bar(x - bar_w / 2, pipeline_vals, width=bar_w,
               color=WONG["vermillion"], edgecolor="black",
               linewidth=0.4,
               label="Python pipeline (7 tools)")
        ax.bar(x + bar_w / 2, graphpop_vals, width=bar_w,
               color=WONG["blue"], edgecolor="black",
               linewidth=0.4,
               label="GraphPop Cypher (1 query)")

        # Per-bar value annotation.
        for i, (pv, gv) in enumerate(zip(pipeline_vals,
                                          graphpop_vals)):
            ax.text(i - bar_w / 2, pv * 1.04,
                    f"{int(pv)}" if pv >= 10 else f"{pv:.1f}",
                    ha="center", va="bottom", fontsize=5.5)
            ax.text(i + bar_w / 2, gv * 1.04,
                    f"{int(gv)}" if gv >= 10 else f"{gv:.1f}",
                    ha="center", va="bottom", fontsize=5.5)

        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=6)
        ax.set_yscale("log")
        ax.set_ylabel("count / seconds (log)")
        ax.set_title("Pipeline-LOC + runtime vs Cypher one-liner")
        ax.legend(frameon=False, loc="upper left", fontsize=5.5)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.text(
            -0.18, 1.05, "c",
            transform=ax.transAxes,
            fontsize=9, fontweight="bold",
            va="top", ha="left",
        )
        out_pdf.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_pdf, format="pdf")
        plt.close(fig)

    return {
        "pipeline": pipeline,
        "graphpop": graphpop,
        "out_pdf": str(out_pdf),
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

DEFAULT_PANELS = Path(
    "/mnt/data/GraphPop/paper/paper2_kinship_arg/benchmarks/fig5_out")
DEFAULT_FIG_DIR = Path(
    "/mnt/data/GraphPop/paper/paper2_kinship_arg/figures")


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--fig5a-template", type=Path,
                   default=DEFAULT_PANELS / "fig5a_cypher_template.txt")
    p.add_argument("--fig5a-annotations", type=Path,
                   default=DEFAULT_PANELS / "fig5a_cypher_annotations.json")
    p.add_argument("--fig5b-csv", type=Path,
                   default=DEFAULT_PANELS / "fig5b_panel_data.csv")
    p.add_argument("--fig5c-csv", type=Path,
                   default=DEFAULT_PANELS / "fig5c_panel_data.csv")
    p.add_argument("--fig5d-json", type=Path,
                   default=DEFAULT_PANELS / "fig5d_ecosystem.json")
    p.add_argument("--fig5a-pdf", type=Path,
                   default=DEFAULT_FIG_DIR / "fig5a.pdf")
    p.add_argument("--fig5b-pdf", type=Path,
                   default=DEFAULT_FIG_DIR / "fig5b.pdf")
    p.add_argument("--fig5c-pdf", type=Path,
                   default=DEFAULT_FIG_DIR / "fig5c.pdf")
    p.add_argument("--fig5d-pdf", type=Path,
                   default=DEFAULT_FIG_DIR / "fig5d.pdf")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    print(f"[fig5a] {make_fig5a(args.fig5a_template, args.fig5a_annotations, args.fig5a_pdf)}")
    if args.fig5b_csv.exists():
        print(f"[fig5b] {make_fig5b(args.fig5b_csv, args.fig5b_pdf)}")
    else:
        print(f"[fig5b] (skipped: {args.fig5b_csv} not found — run fig5b_panels first)")
    if args.fig5c_csv.exists():
        print(f"[fig5c] {make_fig5c(args.fig5c_csv, args.fig5c_pdf)}")
    else:
        print(f"[fig5c] (skipped: {args.fig5c_csv} not found — run fig5c_panels first)")
    print(f"[fig5d] {make_fig5d(args.fig5d_json, args.fig5d_pdf)}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
