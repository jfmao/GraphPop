"""Paper 2 Fig 3a/3b driver — annotation-conditional h² recovery.

Demonstrates Paper 2's G2 claim: a predicate-restricted branch
GRM, when plugged into a heritability estimator, recovers the
true h² inside the predicate and ≈ 0 outside. Phase 1 sim-only,
Python-only.

**Strategy** (per `PLAN_fig3.md`):

1. Simulate `n_arg_draws` independent msprime cohorts.
2. Per ARG: assign the first `m_pathway` mutations to the
   pathway class, the next `m_lof` mutations to the LoF class.
   Compute 5 conditional eGRMs:
   - unconditional (egrm.varGRM_C)
   - pathway-restricted
   - anti-pathway-restricted
   - LoF-class-restricted
   - non-LoF-class-restricted  (= pathway + everything else)
3. Per ARG: draw `n_pheno_replicates` phenotypes for each
   panel:
   - Fig 3a: β supported on pathway mutations only;
              Var(genetic) = true_h²; ε ~ N(0, 1 - true_h²).
   - Fig 3b: β supported on LoF mutations only; same variance
              decomposition.
4. Run Haseman-Elston regression on each GRM × phenotype combo;
   record h²_HE per (ARG, replicate, predicate).
5. Emit `fig3a_panel_data.csv`, `fig3b_panel_data.csv`,
   `fig3_metadata.json`.

The headline observation the figure shows:
    h²_HE(predicate_aligned) ≈ true_h²
    h²_HE(predicate_anti)    ≈ 0

s-LDSC head-to-head is a Phase-2 follow-up; H5 wrapper ships
but the reference LD-score panel build is out-of-scope here.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Pure-logic helpers
# ---------------------------------------------------------------------------

def haseman_elston(
    grm: np.ndarray, y: np.ndarray,
) -> float:
    """Haseman-Elston regression on off-diagonal pairs.

    Solves   y_i y_j ≈ h² · G_ij   over (i < j).

    Returns the OLS slope (= h²_HE). The phenotype is centred
    internally; the GRM is used as-is.
    """
    if grm.ndim != 2 or grm.shape[0] != grm.shape[1]:
        raise ValueError(f"grm shape {grm.shape} must be (n, n)")
    if y.shape != (grm.shape[0],):
        raise ValueError(
            f"y shape {y.shape} must be ({grm.shape[0]},)")
    yc = y - y.mean()
    n = grm.shape[0]
    iu = np.triu_indices(n, k=1)
    g = grm[iu]
    z = np.outer(yc, yc)[iu]
    denom = float(np.sum(g * g))
    if denom == 0:
        return 0.0
    return float(np.sum(z * g) / denom)


@dataclass
class CohortParams:
    n_diploid: int = 50
    sequence_length: int = 20_000
    recomb_rate: float = 1e-5
    mut_rate: float = 1e-4
    population_size: int = 10_000

    @property
    def n_haploid(self) -> int:
        return self.n_diploid * 2


@dataclass
class PhenotypeReplicate:
    arg_idx: int
    replicate_idx: int
    panel: str         # "fig3a" or "fig3b"
    predicate: str     # "unconditional" | "pathway" | "anti_pathway" | ...
    true_h2: float
    h2_estimate: float


def simulate_arg_with_mutations(
    cohort: CohortParams, seed: int,
):
    """msprime sim + mutations; return (ts, genotype_matrix).

    Genotype matrix: (n_haploid, n_mutations) binary array.
    """
    import msprime  # noqa: F401  (local import keeps module importable)
    ts = msprime.sim_ancestry(
        samples=cohort.n_diploid,
        sequence_length=cohort.sequence_length,
        recombination_rate=cohort.recomb_rate,
        population_size=cohort.population_size,
        random_seed=seed,
    )
    ts = msprime.sim_mutations(
        ts, rate=cohort.mut_rate, random_seed=seed)
    G = ts.genotype_matrix().T.astype(np.float64)  # (haploid, mut)
    return ts, G


def child_nodes_of_mutations(ts, mutation_indices: Sequence[int]):
    """Resolve mutation indices → child-node ID set.

    Used by `lit_child_weight` to restrict the eGRM to branches
    bearing the requested mutations.
    """
    children = set()
    for mi in mutation_indices:
        mut = ts.mutation(int(mi))
        children.add(int(mut.node))
    return children


def compute_grm_unconditional(ts) -> np.ndarray:
    """Wrap `egrm.varGRM_C(ts, var=False)` → matrix."""
    from egrm import varGRM_C
    grm, _, _ = varGRM_C(ts, var=False)
    return np.asarray(grm, dtype=np.float64)


def compute_grm_conditional(
    ts, child_set,
) -> np.ndarray:
    """Predicate-restricted GRM via the M4.1 Python reference math.

    Inlines `conditional_egrm` + `lit_child_weight` from
    `graphpop-procedures/src/test/python/build_egrm_fixture.py`
    rather than importing the fixture-builder module — the
    builder's top-level code rewrites fixture files on import,
    which is a side-effect we must NOT trigger from the bench
    harness.
    """
    weight_fn = _lit_child_weight(child_set)
    mat, _ = _conditional_egrm(ts, weight_fn)
    return np.asarray(mat, dtype=np.float64)


def _lit_child_weight(lit_child_set):
    """Per-branch predicate: 1 iff the child node is in the set.

    Mirrors build_egrm_fixture.lit_child_weight verbatim.
    """
    s = set(int(x) for x in lit_child_set)
    def fn(p, c, pt, ct, s_int, e_int):
        return 1.0 if c in s else 0.0
    return fn


def _conditional_egrm(ts, branch_weight_fn):
    """eGRM with a per-branch weight multiplier in [0, 1].

    Mirrors build_egrm_fixture.conditional_egrm verbatim:
    gmap = identity, var=False, rlim=0, alim=inf, left=0,
    right=inf — so the unconditional limit (weight ≡ 1)
    reproduces egrm.varGRM exactly.
    """
    N = ts.num_samples
    mat = np.zeros([N, N])
    total_mu = 0.0
    for tree in ts.trees():
        if tree.total_branch_length == 0:
            continue
        interval_l = tree.interval[1] - tree.interval[0]
        if interval_l <= 0:
            continue
        for c in tree.nodes():
            descendants = list(tree.samples(c))
            n = len(descendants)
            if n == 0 or n == N:
                continue
            parent_c = tree.parent(c)
            if parent_c == -1:
                continue
            parent_time = tree.time(parent_c)
            child_time = tree.time(c)
            t = parent_time - child_time
            if t <= 0:
                continue
            w = branch_weight_fn(
                parent_c, c, parent_time, child_time,
                tree.interval[0], tree.interval[1])
            if w <= 0:
                continue
            mu = interval_l * t * w * 1e-8
            p = n / N
            mat[np.ix_(descendants, descendants)] += (
                mu / (p * (1.0 - p)))
            total_mu += mu
    if total_mu == 0:
        return np.zeros((N, N)), 0.0
    mat /= total_mu
    mat -= mat.mean(axis=0)
    mat -= mat.mean(axis=1, keepdims=True)
    return mat, total_mu


def simulate_phenotype(
    genotype_matrix: np.ndarray,
    causal_mutation_indices: Sequence[int],
    true_h2: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, float]:
    """y = X·β + ε with effects only on `causal_mutation_indices`.

    Returns (y, achieved_h2) where achieved_h2 is the
    empirically realised heritability on this finite sample
    (may differ from `true_h2` by sampling noise).

    Scaling: β_k ~ N(0, σ²_β) with σ²_β chosen so that the
    *empirical* variance of X@β equals `true_h2`. ε is then
    drawn N(0, 1 - true_h2) and y is left un-standardised (HE
    is scale-invariant for h²).
    """
    n_haploid, n_mut = genotype_matrix.shape
    if not causal_mutation_indices:
        raise ValueError("causal_mutation_indices must be non-empty")
    causal = np.asarray(list(causal_mutation_indices), dtype=int)
    if (causal >= n_mut).any() or (causal < 0).any():
        raise ValueError("causal index out of range")

    X_causal = genotype_matrix[:, causal]
    # Centre causal columns so β explains variance, not a mean shift.
    X_c = X_causal - X_causal.mean(axis=0, keepdims=True)

    # Initial unit-variance β; rescale below to hit empirical h².
    beta = rng.normal(size=causal.size)
    g_raw = X_c @ beta
    var_g = float(np.var(g_raw))
    if var_g <= 0:
        # All causal mutations were fixed → degenerate; return zero y.
        return np.zeros(n_haploid), 0.0
    g = g_raw * np.sqrt(true_h2 / var_g)
    eps = rng.normal(size=n_haploid) * np.sqrt(max(0.0, 1.0 - true_h2))
    y = g + eps
    achieved_h2 = float(np.var(g) / max(1e-12, np.var(y)))
    return y, achieved_h2


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

@dataclass
class Fig3ABResult:
    panel_3a_csv: Path
    panel_3b_csv: Path
    metadata_path: Path
    output_dir: Path


PANEL_3A_PREDICATES = ("unconditional", "pathway", "anti_pathway")
PANEL_3B_PREDICATES = ("unconditional", "lof", "non_lof")


def run_fig3ab(
    *,
    cohort: CohortParams,
    output_dir: Path,
    n_arg_draws: int = 5,
    n_pheno_replicates: int = 4,
    m_pathway: int = 30,
    m_lof: int = 30,
    true_h2: float = 0.5,
    seed: int = 2026,
    progress_callback=None,
) -> Fig3ABResult:
    """End-to-end Fig 3a + 3b panel data."""
    output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    arg_seeds = rng.integers(1, 2**31 - 1, size=n_arg_draws).tolist()

    rows_3a: List[PhenotypeReplicate] = []
    rows_3b: List[PhenotypeReplicate] = []

    for arg_idx, arg_seed in enumerate(arg_seeds):
        if progress_callback is not None:
            progress_callback("simulate", arg_idx, n_arg_draws)
        ts, G = simulate_arg_with_mutations(cohort, int(arg_seed))
        n_mut = G.shape[1]
        if n_mut < m_pathway + m_lof:
            raise ValueError(
                f"ARG #{arg_idx} only has {n_mut} mutations; "
                f"need ≥ {m_pathway + m_lof}")
        pathway_idx = list(range(0, m_pathway))
        lof_idx = list(range(m_pathway, m_pathway + m_lof))

        # Predicate child sets.
        pathway_children = child_nodes_of_mutations(ts, pathway_idx)
        lof_children = child_nodes_of_mutations(ts, lof_idx)
        # Complement child sets: every internal node NOT in the
        # positive set. We approximate "anti" as the full mutation
        # set minus the positive set; that's the natural complement
        # in the predicate algebra.
        all_mut_children = child_nodes_of_mutations(
            ts, range(n_mut))
        anti_pathway_children = all_mut_children - pathway_children
        non_lof_children = all_mut_children - lof_children

        # GRMs (compute once per ARG).
        if progress_callback is not None:
            progress_callback("grms", arg_idx, n_arg_draws)
        grm_unc = compute_grm_unconditional(ts)
        grm_pathway = compute_grm_conditional(ts, pathway_children)
        grm_anti_pathway = compute_grm_conditional(
            ts, anti_pathway_children)
        grm_lof = compute_grm_conditional(ts, lof_children)
        grm_non_lof = compute_grm_conditional(ts, non_lof_children)

        # Phenotype replicates for each panel.
        for r in range(n_pheno_replicates):
            # Fig 3a: pathway-driven phenotype.
            y_a, ach_a = simulate_phenotype(
                G, pathway_idx, true_h2, rng)
            for predicate, grm in (
                ("unconditional", grm_unc),
                ("pathway", grm_pathway),
                ("anti_pathway", grm_anti_pathway),
            ):
                h2 = haseman_elston(grm, y_a)
                rows_3a.append(PhenotypeReplicate(
                    arg_idx=arg_idx, replicate_idx=r,
                    panel="fig3a", predicate=predicate,
                    true_h2=ach_a, h2_estimate=h2))

            # Fig 3b: LoF-driven phenotype.
            y_b, ach_b = simulate_phenotype(
                G, lof_idx, true_h2, rng)
            for predicate, grm in (
                ("unconditional", grm_unc),
                ("lof", grm_lof),
                ("non_lof", grm_non_lof),
            ):
                h2 = haseman_elston(grm, y_b)
                rows_3b.append(PhenotypeReplicate(
                    arg_idx=arg_idx, replicate_idx=r,
                    panel="fig3b", predicate=predicate,
                    true_h2=ach_b, h2_estimate=h2))

    panel_3a_csv = output_dir / "fig3a_panel_data.csv"
    panel_3b_csv = output_dir / "fig3b_panel_data.csv"
    _write_panel(rows_3a, panel_3a_csv)
    _write_panel(rows_3b, panel_3b_csv)

    metadata = {
        "cohort": {
            "n_diploid": cohort.n_diploid,
            "n_haploid": cohort.n_haploid,
            "sequence_length": cohort.sequence_length,
            "recomb_rate": cohort.recomb_rate,
            "mut_rate": cohort.mut_rate,
            "population_size": cohort.population_size,
        },
        "n_arg_draws": n_arg_draws,
        "n_pheno_replicates": n_pheno_replicates,
        "m_pathway": m_pathway,
        "m_lof": m_lof,
        "true_h2_target": true_h2,
        "seed": seed,
        "h2_estimator": "haseman_elston",
        "grm_reference": "conditional_egrm (Python reference)",
        "notes": (
            "v1: Python-only HE on predicate-restricted GRMs; "
            "s-LDSC head-to-head comparator deferred to Phase 2."),
    }
    metadata_path = output_dir / "fig3_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2))

    return Fig3ABResult(
        panel_3a_csv=panel_3a_csv,
        panel_3b_csv=panel_3b_csv,
        metadata_path=metadata_path,
        output_dir=output_dir,
    )


def _write_panel(rows: Iterable[PhenotypeReplicate], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "arg_idx", "replicate_idx", "panel", "predicate",
            "true_h2", "h2_estimate",
        ])
        for r in rows:
            w.writerow([
                r.arg_idx, r.replicate_idx, r.panel, r.predicate,
                f"{r.true_h2:.6g}", f"{r.h2_estimate:.6g}",
            ])


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

DEFAULT_OUTPUT = Path(
    "/mnt/data/GraphPop/paper/paper2_kinship_arg/benchmarks/fig3_out")


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--n-diploid", type=int, default=50)
    p.add_argument("--sequence-length", type=int, default=20_000)
    p.add_argument("--n-arg-draws", type=int, default=5)
    p.add_argument("--n-pheno-replicates", type=int, default=4)
    p.add_argument("--m-pathway", type=int, default=30)
    p.add_argument("--m-lof", type=int, default=30)
    p.add_argument("--true-h2", type=float, default=0.5)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--quiet", action="store_true")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    cohort = CohortParams(
        n_diploid=args.n_diploid,
        sequence_length=args.sequence_length,
    )

    def _progress(stage, i, total):
        if not args.quiet:
            print(f"[fig3ab] {stage} {i + 1}/{total}", file=sys.stderr)

    result = run_fig3ab(
        cohort=cohort,
        output_dir=args.output_dir,
        n_arg_draws=args.n_arg_draws,
        n_pheno_replicates=args.n_pheno_replicates,
        m_pathway=args.m_pathway,
        m_lof=args.m_lof,
        true_h2=args.true_h2,
        seed=args.seed,
        progress_callback=None if args.quiet else _progress,
    )
    print(f"[fig3ab] DONE → {result.output_dir}", file=sys.stderr)
    print(f"  fig3a = {result.panel_3a_csv}", file=sys.stderr)
    print(f"  fig3b = {result.panel_3b_csv}", file=sys.stderr)
    print(f"  meta  = {result.metadata_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
