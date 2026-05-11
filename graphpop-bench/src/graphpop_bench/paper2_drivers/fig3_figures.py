"""Paper 2 Fig 3 figure-gen — annotation-conditional GRMs (G2).

Renders three single-column vector PDFs from the panel CSVs:

- **fig3a.pdf** — h² recovery: unconditional / pathway /
                  anti-pathway predicate strip-plots overlaid
                  with the true h² reference line.
- **fig3b.pdf** — same for LoF-class / non-LoF.
- **fig3c.pdf** — composed-predicate rel-err scatter (mirrors
                  Fig 1d/1e visual language).

No msprime / egrm / Java dep — reads CSVs only.
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
    "grey":       "#999999",
}


# ---------------------------------------------------------------------------
# Fig 3a / 3b — h² recovery strip-plots
# ---------------------------------------------------------------------------

def _read_panel_csv(
    csv_path: Path,
) -> Tuple[Dict[str, List[float]], List[float]]:
    """Return (predicate → list of h2_estimate, list of true_h2)."""
    per_predicate: Dict[str, List[float]] = defaultdict(list)
    true_h2s: List[float] = []
    with open(csv_path) as fh:
        for row in csv.DictReader(fh):
            per_predicate[row["predicate"]].append(
                float(row["h2_estimate"]))
            true_h2s.append(float(row["true_h2"]))
    return per_predicate, true_h2s


def _render_h2_panel(
    csv_path: Path, out_pdf: Path,
    *, predicate_order: List[str], colors: Dict[str, str],
    title: str, panel_label: str, ylabel: str | None = None,
) -> dict:
    """Strip-plot of h²_HE across replicates, grouped by predicate."""
    per_predicate, true_h2s = _read_panel_csv(csv_path)
    if not per_predicate:
        raise ValueError(f"no rows in {csv_path}")
    mean_true_h2 = float(np.mean(true_h2s))

    with plt.rc_context(_RCPARAMS):
        fig, ax = plt.subplots(figsize=SINGLE_COL)
        rng = np.random.default_rng(0)
        for i, pred in enumerate(predicate_order):
            vals = np.asarray(per_predicate.get(pred, []))
            if vals.size == 0:
                continue
            jitter = rng.normal(scale=0.06, size=vals.size)
            ax.scatter(
                np.full_like(vals, i, dtype=float) + jitter,
                vals,
                s=10, alpha=0.7,
                facecolor=colors.get(pred, WONG["blue"]),
                edgecolor="black", linewidth=0.25,
                zorder=3,
            )
            # Mean marker.
            ax.scatter(
                [i], [vals.mean()],
                s=30, marker="D", facecolor="white",
                edgecolor=colors.get(pred, WONG["blue"]),
                linewidth=1.0, zorder=4,
            )

        ax.axhline(mean_true_h2, color=WONG["grey"], linewidth=0.5,
                   linestyle="--",
                   label=f"true h² ≈ {mean_true_h2:.2f}")
        ax.axhline(0, color="black", linewidth=0.4, alpha=0.4)
        ax.set_xticks(np.arange(len(predicate_order)))
        ax.set_xticklabels(predicate_order, fontsize=6)
        ax.set_xlabel("predicate")
        ax.set_ylabel(ylabel if ylabel else r"$h^2_{HE}$")
        ax.set_title(title)
        ax.legend(frameon=False, loc="best")
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
        "predicates": list(predicate_order),
        "true_h2": mean_true_h2,
        "out_pdf": str(out_pdf),
        "n_rows_total": sum(len(v) for v in per_predicate.values()),
    }


def make_fig3a(csv_path: Path, out_pdf: Path) -> dict:
    """Fig 3a: pathway-driven phenotype, h² across 3 predicates."""
    return _render_h2_panel(
        csv_path, out_pdf,
        predicate_order=["unconditional", "pathway", "anti_pathway"],
        colors={
            "unconditional": WONG["sky"],
            "pathway": WONG["blue"],
            "anti_pathway": WONG["vermillion"],
        },
        title="Pathway-restricted GRM isolates h²",
        panel_label="a",
    )


def make_fig3b(csv_path: Path, out_pdf: Path) -> dict:
    """Fig 3b: LoF-driven phenotype, h² across 3 predicates."""
    return _render_h2_panel(
        csv_path, out_pdf,
        predicate_order=["unconditional", "lof", "non_lof"],
        colors={
            "unconditional": WONG["sky"],
            "lof": WONG["green"],
            "non_lof": WONG["vermillion"],
        },
        title="LoF-class GRM isolates h²",
        panel_label="b",
    )


# ---------------------------------------------------------------------------
# Fig 3c — composed-predicate rel-err scatter
# ---------------------------------------------------------------------------

def _read_fig3c_csv(
    path: Path,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (ref, gp, rel_err, abs_diff). Same schema as Fig 1d/1e."""
    refs: List[float] = []
    gps: List[float] = []
    rels: List[float] = []
    abss: List[float] = []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            refs.append(float(row["egrm_ref"]))
            gps.append(float(row["graphpop_branch_grm"]))
            rels.append(float(row["rel_err"]))
            abss.append(float(row["abs_diff"]))
    return (np.asarray(refs), np.asarray(gps),
            np.asarray(rels), np.asarray(abss))


def make_fig3c(
    csv_path: Path, out_pdf: Path,
    *, rel_err_gate: float = 1e-6,
) -> dict:
    """Fig 3c — composed-predicate rel-err scatter + inset."""
    with plt.rc_context(_RCPARAMS):
        ref, gp, rel, abs_diff = _read_fig3c_csv(csv_path)
        max_rel = float(np.max(rel)) if len(rel) else 0.0
        max_abs = float(np.max(abs_diff)) if len(abs_diff) else 0.0

        fig, ax = plt.subplots(figsize=SINGLE_COL)
        ax.scatter(
            ref, gp,
            s=8, alpha=0.75,
            facecolor=WONG["purple"], edgecolor="black",
            linewidth=0.25, zorder=3,
        )
        lo = float(min(ref.min(), gp.min())) if len(ref) else 0.0
        hi = float(max(ref.max(), gp.max())) if len(ref) else 1.0
        pad = 0.04 * (hi - lo if hi > lo else 1.0)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad],
                color=WONG["grey"], linewidth=0.5, zorder=1,
                linestyle="--", label="y = x")
        ax.set_xlim(lo - pad, hi + pad)
        ax.set_ylim(lo - pad, hi + pad)
        ax.set_xlabel(r"$\mathrm{conditional\_egrm}$ reference")
        ax.set_ylabel(r"GraphPop predicate-composed branch_grm")
        ax.set_title(r"Predicate composition: pathway $\times$ time_window")
        ax.set_aspect("equal", adjustable="box")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.text(
            -0.18, 1.05, "c",
            transform=ax.transAxes,
            fontsize=9, fontweight="bold",
            va="top", ha="left",
        )
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

        # Rel-err inset.
        ins = fig.add_axes([0.61, 0.20, 0.30, 0.25])
        floor = max(1e-18, float(np.percentile(
            rel[rel > 0], 5)) if (rel > 0).any() else 1e-18)
        positive = np.clip(rel, floor, None)
        ins.hist(
            np.log10(positive), bins=12,
            color=WONG["purple"], edgecolor="black", linewidth=0.3,
        )
        ins.axvline(np.log10(rel_err_gate), color=WONG["grey"],
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


# ---------------------------------------------------------------------------
# Fig 3e — pathway-conditional GRM density: branch vs SNP
# ---------------------------------------------------------------------------

def make_fig3f(
    csv_path: Path, out_pdf: Path,
    *, tolerance: float = 0.05,
) -> dict:
    """Cross-population pathway-h² strip plot."""
    import csv as _csv
    from collections import defaultdict
    by_pop: dict = defaultdict(list)
    by_pop_truth: dict = defaultdict(list)
    with open(csv_path) as fh:
        for row in _csv.DictReader(fh):
            by_pop[row["population"]].append(float(row["h2_estimate"]))
            by_pop_truth[row["population"]].append(float(row["true_h2"]))
    if not by_pop:
        raise ValueError(f"no rows in {csv_path}")

    pops = sorted(by_pop.keys())
    means = [np.mean(by_pop[p]) for p in pops]
    cross_mean = float(np.mean(means))
    truth_mean = float(np.mean([
        v for vals in by_pop_truth.values() for v in vals]))

    colors_seq = [WONG["blue"], WONG["green"],
                   WONG["vermillion"], WONG["orange"],
                   WONG["purple"], WONG["sky"]]

    with plt.rc_context(_RCPARAMS):
        fig, ax = plt.subplots(figsize=SINGLE_COL)
        rng = np.random.default_rng(0)
        for i, pop in enumerate(pops):
            vals = np.asarray(by_pop[pop])
            jitter = rng.normal(scale=0.06, size=vals.size)
            ax.scatter(
                np.full_like(vals, i, dtype=float) + jitter, vals,
                s=10, alpha=0.6,
                facecolor=colors_seq[i % len(colors_seq)],
                edgecolor="black", linewidth=0.25,
                zorder=3,
            )
            ax.scatter(
                [i], [vals.mean()],
                s=40, marker="D", facecolor="white",
                edgecolor=colors_seq[i % len(colors_seq)],
                linewidth=1.0, zorder=4,
            )

        # ± tolerance band around the cross-pop mean.
        ax.axhspan(
            cross_mean - tolerance, cross_mean + tolerance,
            color=WONG["grey"], alpha=0.10, zorder=1,
            label=f"cross-pop ± {tolerance:.2f}",
        )
        ax.axhline(cross_mean, color=WONG["grey"], linewidth=0.5,
                    linestyle="--",
                    label=f"cross-pop mean ({cross_mean:.2f})")
        ax.axhline(truth_mean, color="black", linewidth=0.4,
                    linestyle=":", alpha=0.5,
                    label=f"true h² ({truth_mean:.2f})")

        ax.set_xticks(np.arange(len(pops)))
        ax.set_xticklabels(pops, fontsize=6)
        ax.set_xlabel("population")
        ax.set_ylabel(r"$h^2_{HE}$ (pathway-restricted)")
        ax.set_title("Cross-pop pathway-h² stability")
        ax.legend(frameon=False, loc="upper right", fontsize=5.5)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.text(
            -0.18, 1.05, "f",
            transform=ax.transAxes,
            fontsize=9, fontweight="bold",
            va="top", ha="left",
        )
        out_pdf.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_pdf, format="pdf")
        plt.close(fig)
    return {
        "populations": pops,
        "per_pop_mean": dict(zip(pops, means)),
        "cross_pop_mean": cross_mean,
        "truth_mean": truth_mean,
        "out_pdf": str(out_pdf),
    }


def make_fig3e(csv_path: Path, out_pdf: Path) -> dict:
    """Grouped-bar density panel + correlation line."""
    import csv as _csv
    from collections import defaultdict
    rows = list(_csv.DictReader(open(csv_path)))
    if not rows:
        raise ValueError(f"no rows in {csv_path}")
    by_size: dict = defaultdict(list)
    for r in rows:
        by_size[int(r["pathway_size"])].append(r)

    sizes = sorted(by_size.keys())
    branch_means = [np.mean([float(r["branch_nnz_frac"])
                              for r in by_size[s]]) for s in sizes]
    plink_means = [np.mean([float(r["plink_nnz_frac"])
                             for r in by_size[s]]) for s in sizes]
    branch_stds = [np.std([float(r["branch_nnz_frac"])
                            for r in by_size[s]], ddof=1)
                   if len(by_size[s]) > 1 else 0.0
                   for s in sizes]
    plink_stds = [np.std([float(r["plink_nnz_frac"])
                           for r in by_size[s]], ddof=1)
                  if len(by_size[s]) > 1 else 0.0
                  for s in sizes]
    pearson_means = [np.nanmean([float(r["corr_pearson"])
                                  for r in by_size[s]]) for s in sizes]

    with plt.rc_context(_RCPARAMS):
        fig, ax = plt.subplots(figsize=SINGLE_COL)
        x = np.arange(len(sizes))
        bar_w = 0.35
        ax.bar(x - bar_w/2, branch_means, width=bar_w,
               yerr=branch_stds, capsize=2,
               color=WONG["blue"], edgecolor="black",
               linewidth=0.4,
               error_kw=dict(linewidth=0.5),
               label="GraphPop branch GRM")
        ax.bar(x + bar_w/2, plink_means, width=bar_w,
               yerr=plink_stds, capsize=2,
               color=WONG["vermillion"], edgecolor="black",
               linewidth=0.4,
               error_kw=dict(linewidth=0.5),
               label="PLINK pathway-SNP GRM")
        ax.set_xticks(x)
        ax.set_xticklabels([str(s) for s in sizes])
        ax.set_xlabel("pathway size (# mutations)")
        ax.set_ylabel("off-diag entries with |G| > 1e-6")
        ax.set_ylim(0, 1.05)
        ax.set_title("Pathway-conditional GRM density")
        ax.legend(frameon=False, loc="lower right", fontsize=5.5)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.text(
            -0.18, 1.05, "e",
            transform=ax.transAxes,
            fontsize=9, fontweight="bold",
            va="top", ha="left",
        )

        # Pearson correlation on a twin axis (right).
        ax2 = ax.twinx()
        ax2.plot(x, pearson_means, "-o",
                 color=WONG["green"], markeredgecolor="black",
                 markeredgewidth=0.3, linewidth=0.8,
                 label="branch ↔ PLINK Pearson r")
        ax2.set_ylabel("Pearson r (branch vs PLINK)",
                        color=WONG["green"], fontsize=6)
        ax2.tick_params(axis="y", labelcolor=WONG["green"])
        ax2.set_ylim(-0.05, 1.05)
        ax2.spines["top"].set_visible(False)
        ax2.legend(frameon=False, loc="upper right", fontsize=5.5)

        out_pdf.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_pdf, format="pdf")
        plt.close(fig)
    return {
        "sizes": sizes,
        "branch_nnz_mean": branch_means,
        "plink_nnz_mean": plink_means,
        "pearson_mean": pearson_means,
        "out_pdf": str(out_pdf),
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

DEFAULT_PANELS = Path(
    "/mnt/data/GraphPop/paper/paper2_kinship_arg/benchmarks/fig3_out")
DEFAULT_FIG_DIR = Path(
    "/mnt/data/GraphPop/paper/paper2_kinship_arg/figures")


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--fig3a-csv", type=Path,
                   default=DEFAULT_PANELS / "fig3a_panel_data.csv")
    p.add_argument("--fig3b-csv", type=Path,
                   default=DEFAULT_PANELS / "fig3b_panel_data.csv")
    p.add_argument(
        "--fig3c-csv", type=Path,
        default=(DEFAULT_PANELS
                 / "fig3c_composed_pathway_time_panel_data.csv"))
    p.add_argument("--fig3e-csv", type=Path,
                   default=DEFAULT_PANELS / "fig3e_panel_data.csv")
    p.add_argument("--fig3f-csv", type=Path,
                   default=DEFAULT_PANELS / "fig3f_panel_data.csv")
    p.add_argument("--fig3a-pdf", type=Path,
                   default=DEFAULT_FIG_DIR / "fig3a.pdf")
    p.add_argument("--fig3b-pdf", type=Path,
                   default=DEFAULT_FIG_DIR / "fig3b.pdf")
    p.add_argument("--fig3c-pdf", type=Path,
                   default=DEFAULT_FIG_DIR / "fig3c.pdf")
    p.add_argument("--fig3e-pdf", type=Path,
                   default=DEFAULT_FIG_DIR / "fig3e.pdf")
    p.add_argument("--fig3f-pdf", type=Path,
                   default=DEFAULT_FIG_DIR / "fig3f.pdf")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    print(f"[fig3a] {make_fig3a(args.fig3a_csv, args.fig3a_pdf)}")
    print(f"[fig3b] {make_fig3b(args.fig3b_csv, args.fig3b_pdf)}")
    print(f"[fig3c] {make_fig3c(args.fig3c_csv, args.fig3c_pdf)}")
    if args.fig3e_csv.exists():
        print(f"[fig3e] {make_fig3e(args.fig3e_csv, args.fig3e_pdf)}")
    else:
        print(f"[fig3e] (skipped: {args.fig3e_csv} not found — run fig3e_panels first)")
    if args.fig3f_csv.exists():
        print(f"[fig3f] {make_fig3f(args.fig3f_csv, args.fig3f_pdf)}")
    else:
        print(f"[fig3f] (skipped: {args.fig3f_csv} not found — run fig3f_panels first)")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
