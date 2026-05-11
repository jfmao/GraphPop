"""Unit tests for the R3 Fig 2d/2e real-data driver."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from graphpop_bench.paper2_drivers import fig2de_real_panels as f2der
from graphpop_bench.paper2_drivers.subset_1000g import (
    DEFAULT_PANEL_PATH, PanelRow, load_panel,
)


def test_assign_haploid_pop_labels_two_per_diploid():
    panel = [
        PanelRow(sample="S0", pop="YRI", super_pop="AFR", gender="m"),
        PanelRow(sample="S1", pop="GBR", super_pop="EUR", gender="f"),
        PanelRow(sample="S2", pop="JPT", super_pop="EAS", gender="m"),
    ]
    out = f2der.assign_haploid_pop_labels(
        ["S0", "S1"], panel, use_sub_pop=False)
    assert out == ["AFR", "AFR", "EUR", "EUR"]


def test_assign_haploid_pop_labels_sub_pop():
    panel = [
        PanelRow(sample="S0", pop="YRI", super_pop="AFR", gender="m"),
        PanelRow(sample="S1", pop="GBR", super_pop="EUR", gender="f"),
    ]
    out = f2der.assign_haploid_pop_labels(
        ["S0", "S1"], panel, use_sub_pop=True)
    assert out == ["YRI", "YRI", "GBR", "GBR"]


def test_assign_haploid_pop_labels_unknown_sample_raises():
    panel = [
        PanelRow(sample="S0", pop="YRI", super_pop="AFR", gender="m"),
    ]
    with pytest.raises(KeyError, match="not in panel"):
        f2der.assign_haploid_pop_labels(["S_UNKNOWN"], panel)


# ---------------------------------------------------------------------------
# Live panel test (gated on the 1000G panel being on disk)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(
    not DEFAULT_PANEL_PATH.exists(),
    reason="1000G panel TSV not present",
)
def test_assign_haploid_pop_labels_eur_yri_count():
    """503 EUR + 108 YRI → 1006 + 216 = 1222 haploids."""
    panel = load_panel()
    from graphpop_bench.paper2_drivers.subset_1000g import (
        samples_for_super_pop, samples_for_sub_pop,
    )
    eur = samples_for_super_pop(panel, "EUR")
    yri = samples_for_sub_pop(panel, "YRI")
    diploids = list(eur) + list(yri)
    labels = f2der.assign_haploid_pop_labels(
        diploids, panel, use_sub_pop=False)
    assert len(labels) == 2 * (len(eur) + len(yri))
    assert labels.count("EUR") == 2 * len(eur)
    # YRI samples have super_pop "AFR".
    assert labels.count("AFR") == 2 * len(yri)
