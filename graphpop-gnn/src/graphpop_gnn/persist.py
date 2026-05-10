"""Persist trained sample embeddings back to Neo4j as
:Sample.embedding properties.

Keeps this module torch-free so it tests without GPU deps. Caller
passes in plain numpy / list-of-list arrays.
"""
from __future__ import annotations

from typing import Any, Iterable

import numpy as np


class EmbeddingPersister:
    """Write per-sample embedding vectors to :Sample.embedding."""

    def __init__(self, driver: Any, database: str = "neo4j"):
        self.driver = driver
        self.database = database

    def persist(
        self,
        sample_ids: Iterable[str],
        embeddings: np.ndarray,
        batch_size: int = 1000,
    ) -> int:
        """Write `len(sample_ids)` embedding vectors. Returns the
        number of :Sample nodes updated.
        """
        sids = list(sample_ids)
        if embeddings.ndim != 2:
            raise ValueError(
                "embeddings must be 2-D (n_samples, dim); "
                f"got shape {embeddings.shape}")
        if len(sids) != embeddings.shape[0]:
            raise ValueError(
                "sample_ids length must match embeddings.shape[0]; "
                f"got {len(sids)} ids and {embeddings.shape[0]} rows")
        n = 0
        with self.driver.session(database=self.database) as session:
            for i in range(0, len(sids), batch_size):
                batch = [
                    {
                        "sid": sids[j],
                        "emb": [float(x) for x in embeddings[j]],
                    }
                    for j in range(i, min(i + batch_size, len(sids)))
                ]
                result = session.run(
                    "UNWIND $batch AS r "
                    "MATCH (s:Sample {sampleId: r.sid}) "
                    "SET s.embedding = r.emb "
                    "RETURN count(s) AS c",
                    {"batch": batch},
                )
                rec = result.single()
                if rec is not None:
                    n += int(rec["c"])
        return n

    def clear(self) -> int:
        """Remove every :Sample.embedding property. Returns the number
        of :Sample nodes touched.
        """
        with self.driver.session(database=self.database) as session:
            rec = session.run(
                "MATCH (s:Sample) WHERE s.embedding IS NOT NULL "
                "REMOVE s.embedding RETURN count(s) AS c"
            ).single()
            return int(rec["c"]) if rec is not None else 0
