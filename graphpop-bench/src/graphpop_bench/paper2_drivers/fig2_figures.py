"""Paper 2 Fig 2 figure-gen — posterior branch GRM (G1).

Renders three vector PDFs from the panel CSVs produced by
`fig2_panels.py`:

- **fig2a.pdf** — violin plot of per-pair posterior values.
- **fig2b.pdf** — MSE vs N posterior samples; three curves
                  (posterior, MAP, bootstrap-of-MAP).
- **fig2c.pdf** — empirical 95% CI coverage vs N; horizontal
                  reference at nominal 0.95 + the 0.92
                  success-threshold band.

Compliance with `paper/FIGURE_GUIDELINES.md`: vector PDF,
Arial-fallback sans-serif, Wong palette, single-column width.

No msprime / egrm / Java dep — reads CSVs only.
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

import matplotlib
matplotlib.use("Agg")
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


# ---------------------------------------------------------------------------
# Fig 2a — violin plot of per-pair posterior values
# ---------------------------------------------------------------------------

def make_fig2a(csv_path: Path, out_pdf: Path) -> dict:
    """Read fig2a_pair_posteriors.csv → violin plot."""
    pair_values: Dict[str, List[float]] = defaultdict(list)
    with open(csv_path) as fh:
        for row in csv.DictReader(fh):
            pair_values[row["pair_label"]].append(float(row["value"]))
    if not pair_values:
        raise ValueError(f"no rows in {csv_path}")

    labels = sorted(
        pair_values.keys(),
        key=lambda s: np.median(pair_values[s]),
    )
    data = [pair_values[l] for l in labels]

    with plt.rc_context(_RCPARAMS):
        fig, ax = plt.subplots(figsize=SINGLE_COL)
        parts = ax.violinplot(
            data, positions=np.arange(len(labels)),
            showmeans=True, showmedians=False, widths=0.85)
        for body in parts["bodies"]:
            body.set_facecolor(WONG["blue"])
            body.set_edgecolor("black")
            body.set_alpha(0.55)
            body.set_linewidth(0.4)
        for key in ("cmeans", "cmaxes", "cmins", "cbars"):
            if key in parts:
                parts[key].set_color("black")
                parts[key].set_linewidth(0.5)
        ax.set_xticks(np.arange(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=5)
        ax.set_ylabel(r"branch GRM entry $G_{ij}$")
        ax.set_title("Per-pair posterior distribution")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.text(
            -0.18, 1.05, "a",
            transform=ax.transAxes,
            fontsize=9, fontweight="bold",
            va="top", ha="left",
        )
        out_pdf.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_pdf, format="pdf")
        plt.close(fig)
    return {"n_pairs": len(labels),
            "n_draws_per_pair": len(data[0]) if data else 0,
            "out_pdf": str(out_pdf)}


# ---------------------------------------------------------------------------
# Fig 2b — MSE vs N
# ---------------------------------------------------------------------------

def _read_n_value_csv(
    csv_path: Path, value_columns: List[str],
) -> tuple[np.ndarray, Dict[str, np.ndarray]]:
    """Read a CSV with a header column `n` plus arbitrary value
    columns; return (n_arr, {col: arr}).
    """
    ns: List[int] = []
    cols: Dict[str, List[float]] = {c: [] for c in value_columns}
    with open(csv_path) as fh:
        for row in csv.DictReader(fh):
            try:
                ns.append(int(row["n"]))
            except ValueError:
                continue
            for c in value_columns:
                v = row.get(c, "")
                cols[c].append(float("nan") if v in ("", "nan")
                               else float(v))
    return np.asarray(ns), {c: np.asarray(v) for c, v in cols.items()}


def make_fig2b(csv_path: Path, out_pdf: Path) -> dict:
    """MSE-vs-N panel, log-y."""
    ns, cols = _read_n_value_csv(
        csv_path,
        ["mse_posterior", "mse_map", "mse_bootstrap"])
    if len(ns) == 0:
        raise ValueError(f"no rows in {csv_path}")

    with plt.rc_context(_RCPARAMS):
        fig, ax = plt.subplots(figsize=SINGLE_COL)
        ax.plot(ns, cols["mse_posterior"], "-o",
                color=WONG["blue"],
                markeredgecolor="black", markeredgewidth=0.3,
                label="aggregated posterior")
        ax.plot(ns, cols["mse_map"], "--s",
                color=WONG["vermillion"],
                markeredgecolor="black", markeredgewidth=0.3,
                label="MAP (single ARG)")
        ax.plot(ns, cols["mse_bootstrap"], ":^",
                color=WONG["green"],
                markeredgecolor="black", markeredgewidth=0.3,
                label="bootstrap of MAP")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("N posterior samples")
        ax.set_ylabel("MSE vs ground truth")
        ax.set_title("Posterior aggregation reduces MSE")
        ax.legend(frameon=False, loc="upper right")
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
    return {"n_points": int(len(ns)),
            "min_mse_posterior": float(np.nanmin(cols["mse_posterior"])),
            "out_pdf": str(out_pdf)}


# ---------------------------------------------------------------------------
# Fig 2c — coverage vs N
# ---------------------------------------------------------------------------

def make_fig2c(
    csv_path: Path, out_pdf: Path,
    *, nominal: float = 0.95, success_threshold: float = 0.92,
) -> dict:
    """95% CI coverage panel."""
    ns, cols = _read_n_value_csv(
        csv_path, ["coverage_posterior", "coverage_bootstrap"])
    if len(ns) == 0:
        raise ValueError(f"no rows in {csv_path}")

    with plt.rc_context(_RCPARAMS):
        fig, ax = plt.subplots(figsize=SINGLE_COL)
        # Reference lines first (so data overlays).
        ax.axhline(nominal, color="#999999", linewidth=0.5,
                   linestyle="--",
                   label=f"nominal {nominal:.2f}")
        ax.axhspan(success_threshold, 1.0,
                   color=WONG["green"], alpha=0.07, zorder=0)
        ax.plot(ns, cols["coverage_posterior"], "-o",
                color=WONG["blue"],
                markeredgecolor="black", markeredgewidth=0.3,
                label="aggregated posterior")
        ax.plot(ns, cols["coverage_bootstrap"], ":^",
                color=WONG["vermillion"],
                markeredgecolor="black", markeredgewidth=0.3,
                label="bootstrap of MAP")
        ax.set_xscale("log")
        ax.set_ylim(0, 1.02)
        ax.set_xlabel("N posterior samples")
        ax.set_ylabel("empirical 95% CI coverage")
        ax.set_title("Posterior CI recovers nominal coverage")
        ax.legend(frameon=False, loc="lower right")
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
    return {"n_points": int(len(ns)),
            "max_coverage_posterior": float(
                np.nanmax(cols["coverage_posterior"])),
            "out_pdf": str(out_pdf)}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

DEFAULT_PANELS = Path(
    "/mnt/data/GraphPop/paper/paper2_kinship_arg/benchmarks/fig2_out")
DEFAULT_FIG_DIR = Path(
    "/mnt/data/GraphPop/paper/paper2_kinship_arg/figures")


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--fig2a-csv", type=Path,
                   default=DEFAULT_PANELS / "fig2a_pair_posteriors.csv")
    p.add_argument("--fig2b-csv", type=Path,
                   default=DEFAULT_PANELS / "fig2b_mse_vs_n.csv")
    p.add_argument("--fig2c-csv", type=Path,
                   default=DEFAULT_PANELS / "fig2c_coverage_vs_n.csv")
    p.add_argument("--fig2a-pdf", type=Path,
                   default=DEFAULT_FIG_DIR / "fig2a.pdf")
    p.add_argument("--fig2b-pdf", type=Path,
                   default=DEFAULT_FIG_DIR / "fig2b.pdf")
    p.add_argument("--fig2c-pdf", type=Path,
                   default=DEFAULT_FIG_DIR / "fig2c.pdf")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    print(f"[fig2a] {make_fig2a(args.fig2a_csv, args.fig2a_pdf)}")
    print(f"[fig2b] {make_fig2b(args.fig2b_csv, args.fig2b_pdf)}")
    print(f"[fig2c] {make_fig2c(args.fig2c_csv, args.fig2c_pdf)}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
