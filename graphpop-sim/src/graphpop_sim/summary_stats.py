"""Per-tree-sequence summary statistics for ABC.

Computes π (nucleotide diversity), θ_W (Watterson's theta),
Tajima's D, and mean pairwise F_ST when populations are
defined. Operates directly on a tskit TreeSequence — no Neo4j
dependency.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np


@dataclass
class SummaryStats:
    """Bundle of per-tree-sequence summary statistics."""

    pi: float
    theta_w: float
    tajimas_d: float
    mean_fst: float | None
    n_samples: int
    n_segregating_sites: int
    sequence_length: int

    def as_vector(self, keys: list[str] | None = None) -> np.ndarray:
        """Return the named subset as a numpy vector (for ABC distance)."""
        keys = keys or ["pi", "theta_w", "tajimas_d"]
        d = asdict(self)
        return np.asarray(
            [float(d[k]) if d[k] is not None else 0.0 for k in keys],
            dtype=np.float64,
        )


def _harmonic(n: int) -> float:
    """Harmonic number H_{n-1} (used in Watterson's θ + Tajima's D)."""
    return float(np.sum(1.0 / np.arange(1, n)))


def _tajimas_d(pi: float, theta_w: float, n: int, S: int) -> float:
    """Tajima 1989: D = (π − θ_W) / sqrt(Var). Returns NaN if S=0 or n<4."""
    if n < 4 or S == 0:
        return float("nan")
    a1 = _harmonic(n)
    a2 = float(np.sum(1.0 / np.arange(1, n) ** 2))
    b1 = (n + 1) / (3 * (n - 1))
    b2 = 2 * (n * n + n + 3) / (9 * n * (n - 1))
    c1 = b1 - 1.0 / a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / (a1 * a1)
    e1 = c1 / a1
    e2 = c2 / (a1 * a1 + a2)
    var = e1 * S + e2 * S * (S - 1)
    if var <= 0:
        return float("nan")
    return (pi - theta_w) / float(np.sqrt(var))


def compute_summary_stats(
    ts: "Any",
    populations: list[list[int]] | None = None,
) -> SummaryStats:
    """Compute π / θ_W / Tajima's D + mean F_ST for a tree sequence.

    ts : tskit.TreeSequence
    populations : optional list of lists of sample-index lists; each
                  inner list defines one population. When provided
                  and there are ≥ 2 populations with ≥ 2 samples
                  each, the mean pairwise F_ST is reported.
    """
    n = int(ts.num_samples)
    S = int(ts.num_mutations)
    pi = float(ts.diversity())  # mean π across samples
    theta_w = S / _harmonic(n) if n >= 2 and S > 0 else 0.0
    tajimas_d = _tajimas_d(pi, theta_w, n, S)

    mean_fst: float | None = None
    if populations and len(populations) >= 2:
        valid = [p for p in populations if len(p) >= 2]
        if len(valid) >= 2:
            fst_values = []
            for i in range(len(valid)):
                for j in range(i + 1, len(valid)):
                    try:
                        fst = ts.Fst([valid[i], valid[j]])
                        if hasattr(fst, "__iter__"):
                            fst = float(np.asarray(fst).item())
                        fst_values.append(float(fst))
                    except (ValueError, ZeroDivisionError):
                        continue
            if fst_values:
                mean_fst = float(np.mean(fst_values))

    return SummaryStats(
        pi=pi,
        theta_w=theta_w,
        tajimas_d=tajimas_d,
        mean_fst=mean_fst,
        n_samples=n,
        n_segregating_sites=S,
        sequence_length=int(ts.sequence_length),
    )
