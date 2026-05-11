"""Paper 2 Fig 1d/1e figure-gen — branch-GRM rel-err validation.

Reads `fig1d_panel_data.csv` and `fig1e_panel_data.csv` (produced
by `paper/paper2_kinship_arg/benchmarks/run_fig1de.py`) and emits
two vector PDFs at `paper/paper2_kinship_arg/figures/`:

- `fig1d.pdf` — unconditional branch-GRM scatter, GraphPop vs
                `egrm.varGRM` reference + rel-err inset.
- `fig1e.pdf` — same for the `restrict_to_pathway` conditional
                predicate (pathway_half fixture).

Compliance with `paper/FIGURE_GUIDELINES.md`:
- Vector PDF, embedded TrueType fonts (pdf.fonttype=42).
- Arial 7 pt default, 6 pt ticks.
- Single-column width (3.5 in × 2.625 in).
- Wong / Okabe-Ito palette (`#0072B2` Fig 1d, `#D55E00` Fig 1e).
- Diagonal reference y=x in grey.
- Validation banner reports max rel-err with `< 1e-6 ✓` when met.

The script has no Java / mvn / egrm dependency — it reads only
CSVs, so reproducibility is captured by the panel CSV provenance.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import List, Tuple

import matplotlib
matplotlib.use("Agg")  # headless backend
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Style (mirrors paper/FIGURE_GUIDELINES.md § 4)
# ---------------------------------------------------------------------------

_RCPARAMS = {
    # Arial is the Nature Methods house font; fall back to a
    # generic sans-serif stack if not installed locally so headless
    # CI renders without warnings. Final figures should be
    # re-rendered on a machine with Arial installed for submission.
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
}


# ---------------------------------------------------------------------------
# CSV reading
# ---------------------------------------------------------------------------

def _read_panel_csv(
    path: Path,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (egrm_ref, graphpop, rel_err, abs_diff) as parallel arrays.

    `abs_diff` is read from the precomputed column rather than
    recomputed from the float-formatted `egrm_ref` - `graphpop`
    columns; CSV formatting truncates each at 10 sig figs which
    introduces ~1e-10 round-off noise that would dominate the true
    sub-1e-9 diff we're trying to display.
    """
    if not path.exists():
        raise FileNotFoundError(f"panel CSV not found: {path}")
    refs: List[float] = []
    gps: List[float] = []
    rels: List[float] = []
    abss: List[float] = []
    with open(path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            refs.append(float(row["egrm_ref"]))
            gps.append(float(row["graphpop_branch_grm"]))
            rels.append(float(row["rel_err"]))
            abss.append(float(row["abs_diff"]))
    return (np.asarray(refs), np.asarray(gps),
            np.asarray(rels), np.asarray(abss))


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _render_scatter_with_inset(
    csv_path: Path, out_pdf: Path,
    *,
    color: str,
    title: str,
    panel_label: str,
    rel_err_gate: float = 1e-6,
) -> dict:
    """Single-figure renderer used by both Fig 1d and Fig 1e.

    Layout: main scatter panel (GraphPop vs egrm reference) +
    inset axes in the lower-right showing a log-scaled histogram
    of per-pair rel-err. Diagonal y=x in grey. Validation banner
    in the upper-left.

    Returns a small dict of summary stats so callers can log them.
    """
    with plt.rc_context(_RCPARAMS):
        ref, gp, rel, abs_diff = _read_panel_csv(csv_path)
        max_rel = float(np.max(rel)) if len(rel) else 0.0
        max_abs = float(np.max(abs_diff)) if len(abs_diff) else 0.0

        fig, ax = plt.subplots(figsize=SINGLE_COL)

        # Main scatter: GraphPop vs egrm reference.
        ax.scatter(
            ref, gp,
            s=8, alpha=0.75,
            facecolor=color, edgecolor="black", linewidth=0.25,
            zorder=3,
        )

        # Diagonal reference y=x in grey.
        lo = float(min(ref.min(), gp.min())) if len(ref) else 0.0
        hi = float(max(ref.max(), gp.max())) if len(ref) else 1.0
        pad = 0.04 * (hi - lo if hi > lo else 1.0)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad],
                color="#999999", linewidth=0.5, zorder=1,
                linestyle="--", label="y = x")
        ax.set_xlim(lo - pad, hi + pad)
        ax.set_ylim(lo - pad, hi + pad)
        ax.set_xlabel(r"egrm.varGRM reference")
        ax.set_ylabel(r"GraphPop branch_grm")
        ax.set_title(title)
        ax.set_aspect("equal", adjustable="box")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        # Panel label (Fig 1d / Fig 1e).
        ax.text(
            -0.18, 1.05, panel_label,
            transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="top", ha="left",
        )

        # Validation banner.
        status = "< 1e-6 ✓" if max_rel < rel_err_gate else "FAIL"
        banner = (
            f"max rel-err = {max_rel:.2e}\n"
            f"max abs diff = {max_abs:.2e}\n"
            f"n pairs = {len(ref)}   {status}"
        )
        ax.text(
            0.02, 0.98, banner,
            transform=ax.transAxes,
            fontsize=5.5, va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.3",
                      facecolor="white", edgecolor="#cccccc",
                      linewidth=0.4),
        )

        # Rel-err histogram inset.
        ins = fig.add_axes([0.61, 0.20, 0.30, 0.25])
        # Add a small floor so log10(0) doesn't appear.
        floor = max(1e-18, float(np.percentile(
            rel[rel > 0], 5)) if (rel > 0).any() else 1e-18)
        positive = np.clip(rel, floor, None)
        ins.hist(
            np.log10(positive), bins=12,
            color=color, edgecolor="black", linewidth=0.3,
        )
        ins.axvline(np.log10(rel_err_gate), color="#999999",
                    linestyle="--", linewidth=0.5)
        ins.set_xlabel(r"log$_{10}$ rel-err", fontsize=5.5)
        ins.set_ylabel("count", fontsize=5.5)
        ins.tick_params(axis="both", labelsize=4.5, length=2)
        for sp in ("top", "right"):
            ins.spines[sp].set_visible(False)

        out_pdf.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_pdf, format="pdf")
        plt.close(fig)

        return {
            "n_pairs": int(len(ref)),
            "max_rel_err": max_rel,
            "max_abs_diff": max_abs,
            "status": status,
            "out_pdf": str(out_pdf),
        }


def make_fig1d(csv_path: Path, out_pdf: Path) -> dict:
    """Fig 1d — unconditional branch-GRM rel-err scatter."""
    return _render_scatter_with_inset(
        csv_path, out_pdf,
        color=WONG["blue"],
        title="Unconditional branch GRM",
        panel_label="d",
    )


def make_fig1e(csv_path: Path, out_pdf: Path) -> dict:
    """Fig 1e — conditional (restrict_to_pathway) branch-GRM rel-err scatter."""
    return _render_scatter_with_inset(
        csv_path, out_pdf,
        color=WONG["vermillion"],
        title=r"Conditional branch GRM (restrict_to_pathway)",
        panel_label="e",
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

DEFAULT_REPO = Path("/mnt/data/GraphPop")
DEFAULT_FIG_DIR = DEFAULT_REPO / "paper/paper2_kinship_arg/figures"
DEFAULT_PANELS = (
    DEFAULT_REPO
    / "paper/paper2_kinship_arg/benchmarks/fig1de_out"
)


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument(
        "--fig1d-csv", type=Path,
        default=DEFAULT_PANELS / "fig1d_panel_data.csv")
    p.add_argument(
        "--fig1e-csv", type=Path,
        default=DEFAULT_PANELS / "fig1e_panel_data.csv")
    p.add_argument(
        "--fig1d-pdf", type=Path,
        default=DEFAULT_FIG_DIR / "fig1d.pdf")
    p.add_argument(
        "--fig1e-pdf", type=Path,
        default=DEFAULT_FIG_DIR / "fig1e.pdf")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    s1d = make_fig1d(args.fig1d_csv, args.fig1d_pdf)
    print(f"[fig1d] {s1d}")
    s1e = make_fig1e(args.fig1e_csv, args.fig1e_pdf)
    print(f"[fig1e] {s1e}")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
