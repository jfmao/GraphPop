"""Unit tests for the R4 Fig 3f real-data driver."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from graphpop_bench.paper2_drivers import fig3f_real_panels as f3fr


def _msprime_available() -> bool:
    try:
        import msprime  # noqa: F401
        return True
    except ImportError:
        return False


@pytest.mark.skipif(
    not _msprime_available(),
    reason="msprime required",
)
def test_find_pathway_mutations_in_ts_filters_by_position():
    import msprime
    ts = msprime.sim_ancestry(
        samples=4, sequence_length=5_000,
        recombination_rate=1e-5, random_seed=2)
    ts = msprime.sim_mutations(
        ts, rate=1e-3, random_seed=2,
        model=msprime.BinaryMutationModel())
    all_positions = [int(s.position) for s in ts.sites()]
    pathway_positions = set(all_positions[:3])
    mut_indices, children = f3fr.find_pathway_mutations_in_ts(
        ts, pathway_positions)
    assert len(mut_indices) >= 1
    assert len(children) >= 1
    # Empty pathway set → empty results.
    mi, ch = f3fr.find_pathway_mutations_in_ts(ts, set())
    assert mi == [] and ch == set()


@pytest.mark.skipif(
    not _msprime_available(),
    reason="msprime required",
)
def test_run_fig3f_real_micro(tmp_path):
    """Tiny end-to-end on a synthetic 2-pop msprime tree-sequence."""
    import msprime
    demo = msprime.Demography()
    demo.add_population(name="AFR", initial_size=12_000)
    demo.add_population(name="EUR", initial_size=4_000)
    demo.add_population_split(time=2_000, derived=["EUR"],
                                ancestral="AFR")
    ts = msprime.sim_ancestry(
        samples={"AFR": 5, "EUR": 5}, demography=demo,
        sequence_length=5_000, recombination_rate=1e-5,
        random_seed=42)
    ts = msprime.sim_mutations(
        ts, rate=1e-3, random_seed=42,
        model=msprime.BinaryMutationModel())
    ts_path = tmp_path / "test.trees"
    ts.dump(str(ts_path))

    pop_labels = ["AFR"] * 10 + ["EUR"] * 10
    pathway_positions = set(
        int(s.position) for s in list(ts.sites())[:10])
    result = f3fr.run_fig3f_real(
        ts_path=ts_path,
        pop_labels_per_haploid=pop_labels,
        pathway_positions_in_region=pathway_positions,
        n_pheno_replicates=2,
        true_h2=0.5,
        tolerance=0.05,
        seed=2026,
        output_dir=tmp_path,
    )
    assert result.panel_csv.exists()
    assert result.summary_path.exists()
    with open(result.panel_csv) as fh:
        rows = list(csv.DictReader(fh))
    # 2 pops × 2 pheno replicates = 4 rows.
    assert len(rows) == 4
    pops = {r["population"] for r in rows}
    assert pops == {"AFR", "EUR"}
