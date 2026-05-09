"""Regenerate the posterior-aggregator reference fixture for
PosteriorBranchGrmProcedureTest.

Generates M independent msprime simulations, ingestible ARG JSONs, and
two pairs of expected mean/SD matrices (unconditional + pathway-
restricted). Pinned seeds so the files are deterministic across hosts.

Run (after `pip install egrm tskit msprime` in the graphmana env):

    cd graphpop-procedures/src/test/python
    python build_egrm_fixture_posterior.py
"""
from __future__ import annotations

import json
from pathlib import Path

import msprime  # type: ignore
import numpy as np
import tskit  # type: ignore

from build_egrm_fixture import (
    conditional_egrm,
    lit_child_weight,
    product,
)

HERE = Path(__file__).resolve().parent
OUT_DIR = HERE.parent / "resources"
OUT_DIR.mkdir(parents=True, exist_ok=True)

SEEDS = [42, 43, 44, 45, 46]
N_DIPLOID = 10  # 20 haplotypes
SEQ_LEN = 50_000
RECOMB_RATE = 1e-5
MUT_RATE = 1e-4
MID = 20  # first MID mutations are "lit"

NODE_IS_SAMPLE = 1


def simulate(seed: int):
    ts = msprime.sim_ancestry(
        samples=N_DIPLOID,
        sequence_length=SEQ_LEN,
        recombination_rate=RECOMB_RATE,
        random_seed=seed,
    )
    ts = msprime.sim_mutations(ts, rate=MUT_RATE, random_seed=seed)
    return ts


def dump_arg_json(path: Path, ts, run_id: str):
    """Dump topology + mutations in the same shape as
    egrm_fixture_20samples_arg.json."""
    payload = {
        "schema_version": 1,
        "run_id": run_id,
        "source": "msprime",
        "n_samples": int(ts.num_samples),
        "n_trees": int(ts.num_trees),
        "n_edges": int(ts.num_edges),
        "n_mutations": int(ts.num_mutations),
        "sequence_length": int(ts.sequence_length),
        "nodes": [
            {
                "id": i,
                "time": float(node.time),
                "is_sample": bool(node.flags & NODE_IS_SAMPLE),
                "flags": int(node.flags),
            }
            for i, node in enumerate(ts.nodes())
        ],
        "edges": [
            {
                "parent": int(edge.parent),
                "child": int(edge.child),
                "start": int(edge.left),
                "end": int(edge.right),
            }
            for edge in ts.edges()
        ],
        "samples": list(map(int, ts.samples())),
        "mutations": [],
        "mutation_count": int(ts.num_mutations),
    }
    idx = 0
    for site in ts.sites():
        for mutation in site.mutations:
            payload["mutations"].append({
                "index": idx,
                "position": int(site.position),
                "child_node_id": int(mutation.node),
                "derived_state": str(mutation.derived_state),
            })
            idx += 1
    path.write_text(json.dumps(payload, indent=2))
    return payload


def lit_children_for_run(ts):
    """First MID mutations -> the tskit node id of their child branch."""
    lit = set()
    idx = 0
    for site in ts.sites():
        for mutation in site.mutations:
            if idx < MID:
                lit.add(int(mutation.node))
            idx += 1
    return lit


# -----------------------------------------------------------------------
# Simulate, ingest topology, compute per-run egrm (unconditional + pathway).
# -----------------------------------------------------------------------

uncond_matrices = []
pathway_matrices = []

for seed in SEEDS:
    ts = simulate(seed)
    arg_path = OUT_DIR / f"egrm_posterior_fixture_arg_seed{seed}.json"
    dump_arg_json(arg_path, ts, run_id=f"egrm_post_seed{seed}")

    # Unconditional egrm.
    mat_u, _ = conditional_egrm(ts, lambda *args: 1.0)
    uncond_matrices.append((mat_u + mat_u.T) / 2.0)

    # Pathway-restricted egrm using the same MID-mutations partition
    # as the single-run fixture.
    lit = lit_children_for_run(ts)
    pathway_w = lit_child_weight(lit)
    mat_p, _ = conditional_egrm(ts, product(pathway_w))
    pathway_matrices.append((mat_p + mat_p.T) / 2.0)

    print(f"seed={seed}  trees={ts.num_trees}  muts={ts.num_mutations}  "
          f"lit_branches={len(lit)}  -> {arg_path.name}")


def welford(matrices):
    """Welford-equivalent mean + SE from a list of matrices."""
    arr = np.stack(matrices, axis=0)
    mean = arr.mean(axis=0)
    sd = arr.std(axis=0, ddof=1)
    se = sd / np.sqrt(arr.shape[0])
    return mean, se


def dump_stat(path: Path, label: str, matrix: np.ndarray, kind: str, n_runs: int):
    payload = {
        "schema_version": 1,
        "label": label,
        "kind": kind,           # "mean" or "stderr"
        "n_runs": n_runs,
        "n_samples": int(matrix.shape[0]),
        "matrix": [[float(v) for v in row] for row in matrix],
    }
    path.write_text(json.dumps(payload, indent=2))
    print(f"  -> {path.name}")


mean_u, se_u = welford(uncond_matrices)
mean_p, se_p = welford(pathway_matrices)

dump_stat(OUT_DIR / "egrm_posterior_expected_mean.json",
          label="posterior unconditional", matrix=mean_u,
          kind="mean", n_runs=len(SEEDS))
dump_stat(OUT_DIR / "egrm_posterior_expected_sd.json",
          label="posterior unconditional", matrix=se_u,
          kind="stderr", n_runs=len(SEEDS))
dump_stat(OUT_DIR / "egrm_posterior_pathway_expected_mean.json",
          label="posterior restrict_to_pathway=P_test", matrix=mean_p,
          kind="mean", n_runs=len(SEEDS))
dump_stat(OUT_DIR / "egrm_posterior_pathway_expected_sd.json",
          label="posterior restrict_to_pathway=P_test", matrix=se_p,
          kind="stderr", n_runs=len(SEEDS))
