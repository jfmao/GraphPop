"""Tests for the SLiM runner.

Most tests skip when the SLiM binary isn't on PATH (CI-friendly).
Template loading + listing are exercised regardless.
"""
from __future__ import annotations

import pytest

from graphpop_sim.slim_runner import CURATED_TEMPLATES, SlimRunner


def test_list_templates_returns_curated_three():
    templates = SlimRunner.list_templates()
    assert set(templates) == set(CURATED_TEMPLATES)


def test_load_neutral_template_contents():
    src = SlimRunner.load_template("neutral")
    assert "initializeTreeSeq()" in src
    assert "initializeRecombinationRate(recomb_rate)" in src


def test_load_sweep_template_has_sweep_mutation():
    src = SlimRunner.load_template("sweep_genic")
    assert 'addNewDrawnMutation(m2' in src


def test_load_bottleneck_template_has_setSubpopulationSize():
    src = SlimRunner.load_template("bottleneck")
    assert "setSubpopulationSize" in src


def test_load_unknown_template_raises():
    with pytest.raises(ValueError):
        SlimRunner.load_template("nonsense")


def test_is_available_returns_bool():
    assert isinstance(SlimRunner.is_available(), bool)


@pytest.mark.skipif(
    not SlimRunner.is_available(),
    reason="SLiM binary not on PATH; integration test requires SLiM",
)
def test_neutral_simulation_completes(tmp_path):
    runner = SlimRunner()
    out = tmp_path / "neutral.trees"
    result = runner.run(
        template="neutral",
        params={
            "Ne": 100,
            "sequence_len": 5000,
            "recomb_rate": 1e-7,
            "mut_rate": 1e-7,
            "gens": 200,
            "seed": 1,
        },
        output_path=out,
    )
    assert result.return_code == 0
    assert out.exists()
