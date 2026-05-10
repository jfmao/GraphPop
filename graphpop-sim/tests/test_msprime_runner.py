"""Tests for the msprime runner."""
from __future__ import annotations

import pytest

from graphpop_sim import MsprimeConfig, MsprimeRunner

pytest.importorskip("msprime")


def test_config_from_yaml(tmp_path):
    cfg_path = tmp_path / "cfg.yaml"
    cfg_path.write_text(
        "n_diploid: 5\n"
        "sequence_length: 1000\n"
        "recombination_rate: 1e-7\n"
        "mutation_rate: 1e-6\n"
        "seed: 7\n"
    )
    cfg = MsprimeConfig.from_yaml(cfg_path)
    assert cfg.n_diploid == 5
    assert cfg.sequence_length == 1000
    assert cfg.recombination_rate == 1e-7
    assert cfg.mutation_rate == 1e-6
    assert cfg.seed == 7


def test_config_defaults_when_yaml_empty(tmp_path):
    cfg_path = tmp_path / "empty.yaml"
    cfg_path.write_text("")
    cfg = MsprimeConfig.from_yaml(cfg_path)
    assert cfg.n_diploid == 100
    assert cfg.sequence_length == 1_000_000


def test_simulate_returns_tree_sequence_with_expected_sample_count():
    cfg = MsprimeConfig(
        n_diploid=4,
        sequence_length=10_000,
        recombination_rate=1e-7,
        mutation_rate=1e-6,
        seed=42,
    )
    ts = MsprimeRunner(cfg).simulate()
    # 4 diploid samples → 8 haplotypes.
    assert ts.num_samples == 8
    assert ts.sequence_length == 10_000


def test_simulate_and_dump_writes_a_trees_file(tmp_path):
    cfg = MsprimeConfig(
        n_diploid=4,
        sequence_length=5_000,
        recombination_rate=1e-7,
        mutation_rate=1e-6,
        seed=1,
    )
    out = tmp_path / "sim.trees"
    ts = MsprimeRunner(cfg).simulate_and_dump(out)
    assert out.exists()
    assert ts.num_samples == 8


def test_deterministic_under_same_seed():
    cfg = MsprimeConfig(
        n_diploid=4,
        sequence_length=5_000,
        recombination_rate=1e-7,
        mutation_rate=1e-6,
        seed=99,
    )
    ts1 = MsprimeRunner(cfg).simulate()
    ts2 = MsprimeRunner(cfg).simulate()
    assert ts1.num_mutations == ts2.num_mutations
    assert ts1.num_trees == ts2.num_trees


def test_piecewise_demography():
    cfg = MsprimeConfig(
        n_diploid=4,
        sequence_length=5_000,
        recombination_rate=1e-7,
        mutation_rate=1e-6,
        demography=[
            {"time": 0, "Ne": 1000},
            {"time": 100, "Ne": 10000},
        ],
        seed=42,
    )
    ts = MsprimeRunner(cfg).simulate()
    assert ts.num_samples == 8
