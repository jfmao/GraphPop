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
    p.add_argument("--fig5d-json", type=Path,
                   default=DEFAULT_PANELS / "fig5d_ecosystem.json")
    p.add_argument("--fig5a-pdf", type=Path,
                   default=DEFAULT_FIG_DIR / "fig5a.pdf")
    p.add_argument("--fig5d-pdf", type=Path,
                   default=DEFAULT_FIG_DIR / "fig5d.pdf")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    print(f"[fig5a] {make_fig5a(args.fig5a_template, args.fig5a_annotations, args.fig5a_pdf)}")
    print(f"[fig5d] {make_fig5d(args.fig5d_json, args.fig5d_pdf)}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
