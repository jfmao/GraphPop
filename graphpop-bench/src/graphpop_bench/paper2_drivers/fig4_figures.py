"""Paper 2 Fig 4 figure-gen — biobank-scale scaling (G4).

Renders two single-column vector PDFs from the panel CSV
produced by `fig4_panels.py`:

- **fig4a.pdf** — wall-clock vs n_haploid, log-log.
- **fig4b.pdf** — RSS peak vs n_haploid, log-log.

Each tool gets one line + replicate-aggregated error band.
A least-squares slope annotation appears in the legend.

No msprime / wrapper / Java dep — reads CSVs only.
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .fig4_panels import log_log_fit


_RCPARAMS = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.titlesize": 7,
    "axes.labelsize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "lines.linewidth": 1.0,
    "lines.markersize": 3,
    "savefig.dpi": 600,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}

SINGLE_COL = (3.5, 2.625)

WONG = {
    "blue":       "#0072B2",
    "orange":     "#E69F00",
    "green":      "#009E73",
    "vermillion": "#D55E00",
    "purple":     "#CC79A7",
    "sky":        "#56B4E9",
}

TOOL_COLORS = {
    "tskit":  WONG["blue"],
    "egrm":   WONG["green"],
    "plink":  WONG["orange"],
    "king":   WONG["vermillion"],
    "graphpop": WONG["purple"],
}

TOOL_LABELS = {
    "tskit":    "tskit branch_grm",
    "egrm":     "egrm.varGRM_C",
    "plink":    "PLINK 2.0 GRM",
    "king":     "KING-robust",
    "graphpop": "GraphPop branch_grm",
}


# ---------------------------------------------------------------------------
# CSV reading + aggregation
# ---------------------------------------------------------------------------

def _read_panel_csv(
    path: Path,
) -> Dict[str, Dict[int, List[Tuple[float, float]]]]:
    """Return {tool: {n_haploid: [(wall, rss), ...]}}."""
    out: Dict[str, Dict[int, List[Tuple[float, float]]]] = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            tool = row["tool"]
            n_hap = int(row["n_haploid"])
            wall = float(row["wall_clock_s"])
            rss = float(row["rss_peak_mb"])
            out.setdefault(tool, {}).setdefault(n_hap, []).append(
                (wall, rss))
    return out


def _aggregate_for_plot(
    raw: Dict[int, List[Tuple[float, float]]],
    *, metric_idx: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (n_haploid, mean, half-std) sorted by n_haploid.

    metric_idx: 0 = wall_clock_s, 1 = rss_peak_mb.
    """
    n_hap_values = sorted(raw.keys())
    means = []
    stds = []
    for n in n_hap_values:
        vals = [tup[metric_idx] for tup in raw[n]]
        means.append(float(np.mean(vals)))
        stds.append(float(np.std(vals, ddof=1) if len(vals) > 1
                          else 0.0))
    return (np.asarray(n_hap_values),
            np.asarray(means),
            np.asarray(stds))


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _render_scaling_panel(
    panel_data: Dict[str, Dict[int, List[Tuple[float, float]]]],
    out_pdf: Path,
    *,
    metric_idx: int,
    ylabel: str,
    title: str,
    panel_label: str,
) -> dict:
    """Plot one panel (wall_clock or RSS) for all tools."""
    with plt.rc_context(_RCPARAMS):
        fig, ax = plt.subplots(figsize=SINGLE_COL)

        slopes: Dict[str, float] = {}
        for tool, raw in panel_data.items():
            xs, mean, std = _aggregate_for_plot(
                raw, metric_idx=metric_idx)
            color = TOOL_COLORS.get(tool, WONG["sky"])
            # Slope fit on (xs, mean) — log-log.
            try:
                slope, _ = log_log_fit(list(xs), list(mean))
            except ValueError:
                slope = float("nan")
            slopes[tool] = slope
            label_base = TOOL_LABELS.get(tool, tool)
            label = (
                f"{label_base} (slope {slope:.2f})"
                if not np.isnan(slope) else label_base
            )
            # Errorbars (mean ± half std for visibility).
            ax.errorbar(
                xs, mean, yerr=std * 0.5,
                marker="o", linestyle="-",
                color=color,
                markeredgecolor="black", markeredgewidth=0.3,
                ecolor=color, elinewidth=0.6, capsize=2,
                label=label,
            )

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("n haploid samples")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(frameon=False, loc="upper left", fontsize=5.5)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.text(
            -0.18, 1.05, panel_label,
            transform=ax.transAxes,
            fontsize=9, fontweight="bold",
            va="top", ha="left",
        )

        out_pdf.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_pdf, format="pdf")
        plt.close(fig)
    return {
        "tools": list(panel_data.keys()),
        "slopes": slopes,
        "out_pdf": str(out_pdf),
    }


def make_fig4a(csv_path: Path, out_pdf: Path) -> dict:
    """Wall-clock scaling panel."""
    panel = _read_panel_csv(csv_path)
    return _render_scaling_panel(
        panel, out_pdf,
        metric_idx=0,
        ylabel="wall-clock (s)",
        title="Branch-GRM wall-clock scaling",
        panel_label="a",
    )


def make_fig4b(csv_path: Path, out_pdf: Path) -> dict:
    """RSS peak scaling panel."""
    panel = _read_panel_csv(csv_path)
    return _render_scaling_panel(
        panel, out_pdf,
        metric_idx=1,
        ylabel="RSS peak (MB)",
        title="Branch-GRM memory scaling",
        panel_label="b",
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

DEFAULT_PANELS = Path(
    "/mnt/data/GraphPop/paper/paper2_kinship_arg/benchmarks/fig4_out")
DEFAULT_FIG_DIR = Path(
    "/mnt/data/GraphPop/paper/paper2_kinship_arg/figures")


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--panel-csv", type=Path,
                   default=DEFAULT_PANELS / "fig4_panel_data.csv")
    p.add_argument("--fig4a-pdf", type=Path,
                   default=DEFAULT_FIG_DIR / "fig4a.pdf")
    p.add_argument("--fig4b-pdf", type=Path,
                   default=DEFAULT_FIG_DIR / "fig4b.pdf")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    print(f"[fig4a] {make_fig4a(args.panel_csv, args.fig4a_pdf)}")
    print(f"[fig4b] {make_fig4b(args.panel_csv, args.fig4b_pdf)}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
