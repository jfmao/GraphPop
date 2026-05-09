"""Regenerate the eGRM reference fixture for BranchGrmProcedureTest.

One-shot script. Outputs go into ``../resources/`` and are checked into
git so ``mvn test`` is hermetic (no Python required at test time).

Outputs:
  egrm_fixture_20samples.trees       msprime tree-sequence file
  egrm_expected_20samples.json       {sample_ids, matrix} from egrm.varGRM

Run (after `pip install egrm tskit msprime` in the graphmana env):

    cd graphpop-procedures/src/test/python
    python build_egrm_fixture.py

The msprime seed is pinned so the file is deterministic across machines.
"""
from __future__ import annotations

import json
from pathlib import Path

import egrm  # type: ignore
import msprime  # type: ignore
import numpy as np

# Output dir is ../resources/ relative to this script.
HERE = Path(__file__).resolve().parent
OUT_DIR = HERE.parent / "resources"
OUT_DIR.mkdir(parents=True, exist_ok=True)

TREES_PATH = OUT_DIR / "egrm_fixture_20samples.trees"
JSON_PATH = OUT_DIR / "egrm_expected_20samples.json"
ARG_JSON_PATH = OUT_DIR / "egrm_fixture_20samples_arg.json"

# ---------------------------------------------------------------------------
# Simulate
# ---------------------------------------------------------------------------

N_DIPLOID = 10  # 20 haplotypes
SEQ_LEN = 50_000
RECOMB_RATE = 1e-5
MUT_RATE = 1e-4
SEED = 42

ts = msprime.sim_ancestry(
    samples=N_DIPLOID,
    sequence_length=SEQ_LEN,
    recombination_rate=RECOMB_RATE,
    random_seed=SEED,
)
ts = msprime.sim_mutations(ts, rate=MUT_RATE, random_seed=SEED)
ts.dump(str(TREES_PATH))

print(f"Simulated: {ts.num_samples} samples, {ts.num_trees} trees, "
      f"{ts.num_edges} edges, {ts.num_mutations} mutations.")
print(f"  -> {TREES_PATH}")

# ---------------------------------------------------------------------------
# Reference eGRM via egrm.varGRM
# ---------------------------------------------------------------------------

# var=False to skip variance estimation (we only validate the GRM matrix).
egrm_matrix, _vargrm, total_mu = egrm.varGRM(ts, var=False)

# Ensure exact symmetry (egrm produces a symmetric matrix; floating-point
# accumulation can introduce ~ulp differences that we don't care about).
egrm_sym = (egrm_matrix + egrm_matrix.T) / 2.0
asym = np.max(np.abs(egrm_matrix - egrm_matrix.T))
print(f"  asymmetry max abs: {asym:.3e}")

# ---------------------------------------------------------------------------
# Dump
# ---------------------------------------------------------------------------

sample_ids = list(map(int, ts.samples()))

payload = {
    "schema_version": 1,
    "n_samples": int(ts.num_samples),
    "n_trees": int(ts.num_trees),
    "n_edges": int(ts.num_edges),
    "n_mutations": int(ts.num_mutations),
    "sequence_length": int(ts.sequence_length),
    "msprime_seed": SEED,
    "recombination_rate": RECOMB_RATE,
    "mutation_rate": MUT_RATE,
    "egrm_total_mu": float(total_mu),
    "sample_ids": sample_ids,  # tskit local sample node IDs (haplotypes)
    "matrix": [[float(x) for x in row] for row in egrm_sym],
}

JSON_PATH.write_text(json.dumps(payload, indent=2))
print(f"  -> {JSON_PATH}")
print(f"matrix shape: {len(payload['matrix'])} x {len(payload['matrix'][0])}")
print(f"diagonal min/max: "
      f"{min(payload['matrix'][i][i] for i in range(N_DIPLOID*2)):.4f} / "
      f"{max(payload['matrix'][i][i] for i in range(N_DIPLOID*2)):.4f}")
print(f"off-diag min/max: "
      f"{min(payload['matrix'][i][j] for i in range(N_DIPLOID*2) for j in range(N_DIPLOID*2) if i != j):.4f} / "
      f"{max(payload['matrix'][i][j] for i in range(N_DIPLOID*2) for j in range(N_DIPLOID*2) if i != j):.4f}")

# ---------------------------------------------------------------------------
# Dump the ARG topology as JSON so the Java integration test can rebuild
# it in Neo4j without depending on the Python ARGIngester.
# ---------------------------------------------------------------------------

NODE_IS_SAMPLE = 1  # tskit flag

arg_payload = {
    "schema_version": 1,
    "run_id": "egrm_fixture_20samples",
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
}
ARG_JSON_PATH.write_text(json.dumps(arg_payload, indent=2))
print(f"  -> {ARG_JSON_PATH}")
