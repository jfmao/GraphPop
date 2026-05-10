"""Export the :Sample × :Variant bipartite graph from Neo4j to a
torch-free numpy / npz format.

Downstream consumers (PyTorch Geometric trainers, exploration
notebooks) build whatever HeteroData they need from this dump.
Keeping this module torch-free means it tests cleanly without GPU
deps.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


@dataclass
class ExportedGraph:
    """Plain-Python bipartite graph dump.

    Indexing convention: samples and variants each have packed
    indices in `[0, n)`. The CARRIES edge tensor stores
    `(sample_idx, variant_idx, gt)` triples where `gt ∈ {1, 2}`
    (heterozygous, homozygous-alt). The optional `populations`
    array stores per-sample population labels (string).
    """

    sample_ids: list[str] = field(default_factory=list)
    variant_ids: list[str] = field(default_factory=list)
    sample_idx: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int64))
    variant_idx: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int64))
    gt: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int8))
    populations: list[str] | None = None

    @property
    def n_samples(self) -> int:
        return len(self.sample_ids)

    @property
    def n_variants(self) -> int:
        return len(self.variant_ids)

    @property
    def n_edges(self) -> int:
        return int(self.sample_idx.shape[0])

    def save(self, path: str | Path) -> None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            path,
            sample_ids=np.asarray(self.sample_ids, dtype="U"),
            variant_ids=np.asarray(self.variant_ids, dtype="U"),
            sample_idx=self.sample_idx,
            variant_idx=self.variant_idx,
            gt=self.gt,
            populations=(np.asarray(self.populations, dtype="U")
                          if self.populations is not None
                          else np.zeros(0, dtype="U")),
        )

    @classmethod
    def load(cls, path: str | Path) -> "ExportedGraph":
        d = np.load(path, allow_pickle=False)
        pops = d["populations"]
        return cls(
            sample_ids=list(d["sample_ids"]),
            variant_ids=list(d["variant_ids"]),
            sample_idx=d["sample_idx"],
            variant_idx=d["variant_idx"],
            gt=d["gt"],
            populations=list(pops) if pops.size > 0 else None,
        )


class GraphExporter:
    """Pull samples / variants / CARRIES from a Neo4j driver.

    The driver protocol is duck-typed (anything with a `.session(database=...)`
    that returns a context manager exposing `.run(cypher, params)` works);
    tests pass an in-memory mock to avoid a live Neo4j.
    """

    def __init__(self, driver: Any, database: str = "neo4j"):
        self.driver = driver
        self.database = database

    def export(self,
                min_samples: int = 1,
                limit_variants: int | None = None) -> ExportedGraph:
        """Pull the full bipartite graph.

        min_samples : keep only variants with at least this many CARRIES
                       edges (drops singletons).
        limit_variants : cap the number of variants for development /
                          smoke runs (None = no cap).
        """
        with self.driver.session(database=self.database) as session:
            sample_rows = list(session.run(
                "MATCH (s:Sample) "
                "RETURN s.sampleId AS sid, s.population AS pop "
                "ORDER BY s.sampleId"
            ))
            sample_ids = [r["sid"] for r in sample_rows]
            populations = [r.get("pop") or "" for r in sample_rows]
            sid_to_idx = {sid: i for i, sid in enumerate(sample_ids)}

            limit_clause = f" LIMIT {limit_variants}" if limit_variants else ""
            variant_rows = list(session.run(
                "MATCH (v:Variant) "
                "WITH v, size((:Sample)-[:CARRIES]->(v)) AS n_carriers "
                "WHERE n_carriers >= $min "
                "RETURN v.variantId AS vid "
                "ORDER BY v.variantId" + limit_clause,
                {"min": min_samples},
            ))
            variant_ids = [r["vid"] for r in variant_rows]
            vid_to_idx = {vid: i for i, vid in enumerate(variant_ids)}

            edge_rows = list(session.run(
                "MATCH (s:Sample)-[r:CARRIES]->(v:Variant) "
                "WHERE v.variantId IN $vids "
                "RETURN s.sampleId AS sid, v.variantId AS vid, "
                "r.gt AS gt",
                {"vids": variant_ids},
            ))

        s_idx = np.empty(len(edge_rows), dtype=np.int64)
        v_idx = np.empty(len(edge_rows), dtype=np.int64)
        gt = np.empty(len(edge_rows), dtype=np.int8)
        keep = 0
        for r in edge_rows:
            si = sid_to_idx.get(r["sid"])
            vi = vid_to_idx.get(r["vid"])
            if si is None or vi is None:
                continue
            s_idx[keep] = si
            v_idx[keep] = vi
            gt[keep] = int(r["gt"]) if r["gt"] is not None else 1
            keep += 1
        s_idx = s_idx[:keep]
        v_idx = v_idx[:keep]
        gt = gt[:keep]

        return ExportedGraph(
            sample_ids=sample_ids,
            variant_ids=variant_ids,
            sample_idx=s_idx,
            variant_idx=v_idx,
            gt=gt,
            populations=populations if any(populations) else None,
        )
