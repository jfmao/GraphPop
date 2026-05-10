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

TIME_WINDOW_PATH = OUT_DIR / "egrm_expected_time_window.json"
PATHWAY_PATH = OUT_DIR / "egrm_expected_pathway_half.json"
CONSEQUENCE_PATH = OUT_DIR / "egrm_expected_consequence_missense.json"
COMPOSED_PATH = OUT_DIR / "egrm_expected_composed_pathway_time.json"

# Deterministic mutation partitions for the conditional-predicate fixtures.
# Mutations 0..MID-1 are "in pathway P_test" / "missense"; the rest are not.
MID = 20  # fixture has 41 mutations; 20 are "lit", 21 are "dark"

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
    "mutations": [
        {
            "index": idx,
            "position": int(site.position),
            "child_node_id": int(mutation.node),
            "derived_state": str(mutation.derived_state),
        }
        for site in ts.sites()
        for idx, mutation in enumerate(site.mutations)
        # NB: enumerate restarts per-site so this is wrong for multi-mutation
        # sites; the fixture has at most one mutation per site (infinite-sites
        # default), so it's fine here.
    ],
    # Convenience: monotonic mutation index across the whole tree sequence.
    "mutation_count": int(ts.num_mutations),
}
# Replace the per-site index with a global mutation index for unambiguous
# downstream referencing.
mut_global_idx = 0
for m in arg_payload["mutations"]:
    m["index"] = mut_global_idx
    mut_global_idx += 1
ARG_JSON_PATH.write_text(json.dumps(arg_payload, indent=2))
print(f"  -> {ARG_JSON_PATH}")

# ---------------------------------------------------------------------------
# Conditional eGRM reference (mirrors egrm.varGRM source, plus per-branch
# weight). Used by BranchGrmProcedureTest to validate time_window /
# restrict_to_pathway / mutation_filter / composed predicates.
# ---------------------------------------------------------------------------


def conditional_egrm(ts, branch_weight_fn):
    """eGRM with a per-branch weight multiplier in [0, 1].

    Mirrors egrm.varGRM exactly when branch_weight_fn returns 1
    (gmap = identity, var=False, rlim=0, alim=inf, left=0, right=inf).
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
                continue  # root
            parent_time = tree.time(parent_c)
            child_time = tree.time(c)
            t = parent_time - child_time
            if t <= 0:
                continue
            w = branch_weight_fn(parent_c, c, parent_time, child_time,
                                 tree.interval[0], tree.interval[1])
            if w <= 0:
                continue
            mu = interval_l * t * w * 1e-8
            p = n / N
            mat[np.ix_(descendants, descendants)] += mu / (p * (1.0 - p))
            total_mu += mu

    if total_mu == 0:
        return np.zeros((N, N)), 0.0

    mat /= total_mu
    mat -= mat.mean(axis=0)
    mat -= mat.mean(axis=1, keepdims=True)
    return mat, total_mu


def time_window_weight(t_lo, t_hi):
    def fn(p, c, pt, ct, s, e):
        denom = pt - ct
        if denom <= 0:
            return 0.0
        overlap = min(pt, t_hi) - max(ct, t_lo)
        return max(0.0, overlap) / denom
    return fn


def lit_child_weight(lit_child_set):
    def fn(p, c, pt, ct, s, e):
        return 1.0 if c in lit_child_set else 0.0
    return fn


def product(*fns):
    def fn(*args):
        w = 1.0
        for f in fns:
            w *= f(*args)
            if w <= 0:
                return 0.0
        return w
    return fn


def dump_conditional(path, label, t_lo=None, t_hi=None,
                     pathway_lit_children=None,
                     consequence_lit_children=None,
                     extra_meta=None):
    fns = []
    if t_lo is not None and t_hi is not None:
        fns.append(time_window_weight(t_lo, t_hi))
    if pathway_lit_children is not None:
        fns.append(lit_child_weight(pathway_lit_children))
    if consequence_lit_children is not None:
        fns.append(lit_child_weight(consequence_lit_children))
    if not fns:
        fns.append(lambda *_: 1.0)
    weight_fn = product(*fns)

    mat, tmu = conditional_egrm(ts, weight_fn)
    sym = (mat + mat.T) / 2.0
    payload = {
        "schema_version": 1,
        "label": label,
        "n_samples": int(ts.num_samples),
        "egrm_total_mu": float(tmu),
        "matrix": [[float(x) for x in row] for row in sym],
    }
    if extra_meta:
        payload.update(extra_meta)
    path.write_text(json.dumps(payload, indent=2))
    print(f"  -> {path}  (total_mu={tmu:.4e})")


# Pick a time window that splits internal-node times roughly in half.
internal_times = sorted(
    [n.time for i, n in enumerate(ts.nodes()) if n.time > 0])
T_HI = float(internal_times[len(internal_times) // 2]) if internal_times else 0.5
T_LO = 0.0
print(f"\nConditional fixtures (time_window upper = median internal time = {T_HI:.4g})")

# Map mutation index -> child tskit node id (the branch the mutation lands on).
mut_to_child = []
for site in ts.sites():
    for mutation in site.mutations:
        mut_to_child.append(int(mutation.node))

pathway_lit = {mut_to_child[i] for i in range(min(MID, len(mut_to_child)))}
consequence_lit = pathway_lit  # same partition for the missense fixture

dump_conditional(
    TIME_WINDOW_PATH,
    label="time_window",
    t_lo=T_LO, t_hi=T_HI,
    extra_meta={"t_lo": T_LO, "t_hi": T_HI},
)
dump_conditional(
    PATHWAY_PATH,
    label="restrict_to_pathway=P_test",
    pathway_lit_children=pathway_lit,
    extra_meta={
        "pathway_id": "P_test",
        "lit_mutation_indices": list(range(MID)),
        "lit_child_node_ids": sorted(pathway_lit),
    },
)
dump_conditional(
    CONSEQUENCE_PATH,
    label="mutation_filter=missense_variant",
    consequence_lit_children=consequence_lit,
    extra_meta={
        "consequence": "missense_variant",
        "lit_mutation_indices": list(range(MID)),
        "lit_child_node_ids": sorted(consequence_lit),
    },
)
dump_conditional(
    COMPOSED_PATH,
    label="restrict_to_pathway=P_test & time_window",
    t_lo=T_LO, t_hi=T_HI,
    pathway_lit_children=pathway_lit,
    extra_meta={
        "pathway_id": "P_test",
        "t_lo": T_LO, "t_hi": T_HI,
        "lit_child_node_ids": sorted(pathway_lit),
    },
)

# ---------------------------------------------------------------------------
# M4.B + step 7: ancestry-decomposed eGRM reference.
# Synthetic painting: samples 0-9 = EUR, 10-19 = AFR. Internal nodes
# painted via majority-vote DFS over descendants.
# ---------------------------------------------------------------------------

PAINTING_PATH = OUT_DIR / "egrm_by_ancestry_fixture_painting.json"
DECOMP_EUR_PATH = OUT_DIR / "egrm_by_ancestry_expected_EUR.json"
DECOMP_AFR_PATH = OUT_DIR / "egrm_by_ancestry_expected_AFR.json"
DECOMP_PATHWAY_EUR_PATH = OUT_DIR / "egrm_by_ancestry_pathway_expected_EUR.json"
DECOMP_PATHWAY_AFR_PATH = OUT_DIR / "egrm_by_ancestry_pathway_expected_AFR.json"


def propagate_painting_majority(treeseq, sample_to_anc):
    """Majority-vote DFS propagation matching AncestryIngester."""
    from collections import Counter

    children = {}
    for edge in treeseq.edges():
        children.setdefault(int(edge.parent), []).append(int(edge.child))
    is_sample = {i: bool(node.flags & 1) for i, node in enumerate(treeseq.nodes())}

    desc_counts = {}

    def dfs(n):
        if n in desc_counts:
            return desc_counts[n]
        if is_sample.get(n, False):
            anc = sample_to_anc.get(n)
            counter = Counter()
            if anc is not None:
                counter[anc] = 1
            desc_counts[n] = counter
            return counter
        counter = Counter()
        for c in children.get(n, []):
            counter.update(dfs(c))
        desc_counts[n] = counter
        return counter

    # Emit ALL label probabilities (fractional painting) so the
    # decomposition partition Σ_a prob_a(node) = 1 is preserved.
    painting = {}
    for i in range(treeseq.num_nodes):
        dfs(i)
    for i in range(treeseq.num_nodes):
        counts = desc_counts[i]
        if not counts:
            continue
        total = sum(counts.values())
        painting[i] = [(label, counts[label] / total)
                        for label in sorted(counts)]
    return painting


def decompose_egrm_by_ancestry(treeseq, painting, branch_weight_fn):
    """Bucketed per-ancestry egrm; shared total_mu denominator.

    Mirrors the Java BranchGrmByAncestryComputer kernel.
    """
    N = treeseq.num_samples
    ancestries = sorted({a for entries in painting.values() for a, _ in entries})
    mats = {a: np.zeros([N, N]) for a in ancestries}
    total_mu = 0.0

    for tree in treeseq.trees():
        if tree.total_branch_length == 0:
            continue
        l = tree.interval[1] - tree.interval[0]
        if l <= 0:
            continue
        for c in tree.nodes():
            descendants = list(tree.samples(c))
            n = len(descendants)
            if n == 0 or n == N:
                continue
            parent_c = tree.parent(c)
            if parent_c == -1:
                continue
            t = tree.time(parent_c) - tree.time(c)
            if t <= 0:
                continue
            w = branch_weight_fn(parent_c, c, tree.time(parent_c), tree.time(c),
                                 tree.interval[0], tree.interval[1])
            if w <= 0:
                continue
            mu = l * t * w * 1e-8
            p = n / N
            update = mu / (p * (1 - p))

            child_paint = painting.get(c)
            if child_paint is not None:
                for anc, prob in child_paint:
                    if anc in mats and prob > 0:
                        mats[anc][np.ix_(descendants, descendants)] += prob * update

            total_mu += mu

    if total_mu == 0:
        return {a: np.zeros((N, N)) for a in ancestries}

    out = {}
    for a, m in mats.items():
        m /= total_mu
        m -= m.mean(axis=0)
        m -= m.mean(axis=1, keepdims=True)
        out[a] = (m + m.T) / 2.0
    return out


def dump_decomp(path, label, matrix, ancestry, n_runs=1):
    payload = {
        "schema_version": 1,
        "label": label,
        "ancestry": ancestry,
        "n_runs": n_runs,
        "n_samples": int(matrix.shape[0]),
        "matrix": [[float(v) for v in row] for row in matrix],
    }
    path.write_text(json.dumps(payload, indent=2))
    print(f"  -> {path.name}")


# Synthetic haplotype-ancestry partition (samples 0-9 = EUR, 10-19 = AFR).
sample_to_anc = {i: ("EUR" if i < 10 else "AFR") for i in range(20)}
painting = propagate_painting_majority(ts, sample_to_anc)

# Dump painting for the Java integration test setUp.
_painting_rows = [
    {"tskit_node_id": int(nid), "population_id": label, "posterior_prob": float(prob)}
    for nid in sorted(painting)
    for label, prob in painting[nid]
]
PAINTING_PATH.write_text(json.dumps({
    "schema_version": 1,
    "run_id": "egrm_fixture_20samples",
    "painter": "majority_vote",
    "rows": _painting_rows,
}, indent=2))
print(f"  -> {PAINTING_PATH.name}  ({len(painting)} nodes painted, "
      f"{len(_painting_rows)} edges)")

# Unconditional decomposition.
mats = decompose_egrm_by_ancestry(ts, painting, lambda *_: 1.0)
dump_decomp(DECOMP_EUR_PATH, "egrm_by_ancestry uncond", mats["EUR"], "EUR")
dump_decomp(DECOMP_AFR_PATH, "egrm_by_ancestry uncond", mats["AFR"], "AFR")

# Pathway-restricted decomposition.
pathway_w = lit_child_weight(pathway_lit)
mats_p = decompose_egrm_by_ancestry(ts, painting, pathway_w)
dump_decomp(DECOMP_PATHWAY_EUR_PATH, "egrm_by_ancestry pathway",
            mats_p["EUR"], "EUR")
dump_decomp(DECOMP_PATHWAY_AFR_PATH, "egrm_by_ancestry pathway",
            mats_p["AFR"], "AFR")

# ---------------------------------------------------------------------------
# M4.3 -- ARG-derived IBD reference via tskit.TreeSequence.ibd_segments.
# ---------------------------------------------------------------------------

IBD_PATH = OUT_DIR / "egrm_fixture_20samples_ibd.json"

# All-pairs, no TMRCA cap, no min span -- the full reference set.
ibd_result = ts.ibd_segments(within=ts.samples(),
                             store_pairs=True,
                             store_segments=True)

ibd_payload = {
    "schema_version": 1,
    "n_samples": int(ts.num_samples),
    "max_time": None,        # i.e. infinity
    "min_span": 0,
    "segments": [],
}
for pair, seglist in ibd_result.items():
    a, b = sorted((int(pair[0]), int(pair[1])))
    for seg in seglist:
        ibd_payload["segments"].append({
            "sample_a": a,
            "sample_b": b,
            "start": int(seg.left),
            "end": int(seg.right),
            "mrca_node_id": int(seg.node),
            "tmrca": float(ts.node(seg.node).time),
        })

# Sort deterministically.
ibd_payload["segments"].sort(key=lambda s: (s["sample_a"], s["sample_b"],
                                              s["start"]))
IBD_PATH.write_text(json.dumps(ibd_payload, indent=2))
print(f"  -> {IBD_PATH.name}  ({len(ibd_payload['segments'])} segments)")
