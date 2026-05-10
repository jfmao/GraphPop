"""msprime orchestrator (M12.A).

Reads a YAML config, runs `msprime.sim_ancestry` +
`msprime.sim_mutations`, dumps the resulting tree sequence, and
optionally invokes the existing ``graphpop_import.ARGIngester``
to persist it as a new ``:ARGRun``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass
class MsprimeConfig:
    """Per-run simulation parameters."""

    n_diploid: int = 100
    sequence_length: int = 1_000_000
    recombination_rate: float = 1e-8
    mutation_rate: float = 1e-8
    demography: list[dict[str, Any]] = field(default_factory=list)
    seed: int = 42

    @classmethod
    def from_yaml(cls, path: str | Path) -> "MsprimeConfig":
        with open(path) as fh:
            data = yaml.safe_load(fh) or {}
        return cls(
            n_diploid=int(data.get("n_diploid", 100)),
            sequence_length=int(data.get("sequence_length", 1_000_000)),
            recombination_rate=float(data.get("recombination_rate", 1e-8)),
            mutation_rate=float(data.get("mutation_rate", 1e-8)),
            demography=list(data.get("demography") or []),
            seed=int(data.get("seed", 42)),
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "n_diploid": self.n_diploid,
            "sequence_length": self.sequence_length,
            "recombination_rate": self.recombination_rate,
            "mutation_rate": self.mutation_rate,
            "demography": list(self.demography),
            "seed": self.seed,
        }


class MsprimeRunner:
    """Run an msprime simulation given an ``MsprimeConfig``."""

    def __init__(self, config: MsprimeConfig):
        self.config = config

    def simulate(self) -> "Any":
        """Run sim_ancestry + sim_mutations; returns a tskit TreeSequence."""
        import msprime

        demography = self._build_demography()
        ts = msprime.sim_ancestry(
            samples=self.config.n_diploid,
            sequence_length=self.config.sequence_length,
            recombination_rate=self.config.recombination_rate,
            random_seed=self.config.seed,
            demography=demography,
        )
        ts = msprime.sim_mutations(
            ts, rate=self.config.mutation_rate,
            random_seed=self.config.seed,
        )
        return ts

    def simulate_and_dump(self, output_path: str | Path) -> "Any":
        """Run + dump tree sequence to ``output_path``. Returns the ts."""
        ts = self.simulate()
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        ts.dump(str(output_path))
        return ts

    def _build_demography(self):
        """Convert the YAML demography list into a msprime.Demography.

        Each entry is ``{time: t, Ne: n}`` — a piecewise-constant
        Ne trajectory. ``time: 0`` defines present-day Ne; later
        entries are population-size changes (older → larger time).
        v1 supports a single deme; multi-deme migration is deferred.
        """
        import msprime

        demography = msprime.Demography()
        epochs = sorted(self.config.demography, key=lambda d: d["time"])
        if not epochs:
            demography.add_population(name="A", initial_size=10_000)
            return demography
        present_ne = float(epochs[0]["Ne"])
        demography.add_population(name="A", initial_size=present_ne)
        for ep in epochs[1:]:
            demography.add_population_parameters_change(
                time=float(ep["time"]),
                initial_size=float(ep["Ne"]),
                population="A",
            )
        return demography
