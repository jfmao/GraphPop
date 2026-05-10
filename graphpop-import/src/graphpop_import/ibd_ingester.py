"""External-tool IBD-segment ingest into the GraphPop ARG layer.

Populates ``:Sample-[:IBD_SEGMENT {chr, start, end, length_bp,
length_cM, source, created_at}]->:Sample`` from a hap-IBD-style TSV
or programmatic input. Mirrors the additive-ingest pattern from
:class:`graphpop_import.arg_importer.ARGIngester` and
:class:`graphpop_import.ancestry_ingester.AncestryIngester`.

Closes the matrix-side leg of M4.3: biobank pipelines (UK Biobank,
TOPMed, All-of-Us) ship hap-IBD output as the canonical biobank
relationship signal; this module wires that signal into Neo4j as
queryable :IBD_SEGMENT edges.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import TYPE_CHECKING, Any, Iterable, Iterator, Optional

if TYPE_CHECKING:  # pragma: no cover
    from neo4j import Driver, ManagedTransaction


@dataclass(frozen=True)
class IBDSegmentRow:
    """Single IBD segment between two samples."""

    sample_a: str
    sample_b: str
    chr: str
    start: int
    end: int
    length_cM: Optional[float] = None

    def length_bp(self) -> int:
        return int(self.end - self.start)

    def canonical(self) -> "IBDSegmentRow":
        """Return a copy with sample_a < sample_b lexicographically."""
        if self.sample_a <= self.sample_b:
            return self
        return IBDSegmentRow(
            sample_a=self.sample_b,
            sample_b=self.sample_a,
            chr=self.chr,
            start=self.start,
            end=self.end,
            length_cM=self.length_cM,
        )


@dataclass
class IBDSourceSummary:
    source: str
    n_edges: int
    n_pairs: int


@dataclass
class IBDIngestSummary:
    source: str
    n_segments: int
    n_pairs: int
    created_at: datetime = field(
        default_factory=lambda: datetime.now(timezone.utc))


class IBDIngester:
    """Additive Cypher ingest of ``:IBD_SEGMENT`` edges from external
    callers (hap-IBD, GERMLINE, iLASH, …)."""

    def __init__(
        self,
        driver: "Driver",
        *,
        database: str = "neo4j",
        batch_size: int = 10_000,
    ) -> None:
        self._driver = driver
        self._database = database
        self._batch_size = batch_size

    # ---- public API --------------------------------------------------

    def ingest(
        self,
        segments: Iterable[IBDSegmentRow],
        *,
        source: str = "hap_ibd",
        replace: bool = False,
    ) -> IBDIngestSummary:
        rows = [s.canonical() for s in segments]
        pairs = {(s.sample_a, s.sample_b) for s in rows}

        with self._driver.session(database=self._database) as session:
            if replace:
                session.execute_write(self._delete_source, source)
            else:
                session.execute_write(self._fail_if_source_exists, source)
            for batch in _chunked(rows, self._batch_size):
                session.execute_write(
                    self._create_ibd_batch, source, batch)

        return IBDIngestSummary(
            source=source,
            n_segments=len(rows),
            n_pairs=len(pairs),
        )

    def ingest_tsv(
        self,
        path: str,
        *,
        source: str = "hap_ibd",
        has_header: bool = False,
        replace: bool = False,
    ) -> IBDIngestSummary:
        """Parse a hap-IBD-style 6-column TSV.

        Column order (no header by default):

            sample_a<TAB>sample_b<TAB>chr<TAB>start_bp<TAB>end_bp<TAB>length_cM

        ``length_cM`` is optional; if omitted (5 cols), it is recorded
        as ``None`` and downstream IBD-kinship will use the ``ibd_bp``
        fallback path.
        """
        rows: list[IBDSegmentRow] = []
        with open(path) as fh:
            for lineno, raw in enumerate(fh, start=1):
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                if has_header and lineno == 1:
                    continue
                parts = line.split("\t")
                if len(parts) < 5:
                    raise ValueError(
                        f"line {lineno}: expected ≥5 TSV columns, got: {raw!r}"
                    )
                length_cM = (float(parts[5])
                             if len(parts) >= 6 and parts[5] != ""
                             else None)
                rows.append(IBDSegmentRow(
                    sample_a=parts[0],
                    sample_b=parts[1],
                    chr=parts[2],
                    start=int(parts[3]),
                    end=int(parts[4]),
                    length_cM=length_cM,
                ))
        return self.ingest(rows, source=source, replace=replace)

    def list_sources(self) -> list[IBDSourceSummary]:
        with self._driver.session(database=self._database) as session:
            return session.execute_read(self._list_sources)

    def delete_source(self, source: str) -> int:
        with self._driver.session(database=self._database) as session:
            return session.execute_write(self._delete_source, source)

    # ---- internal ----------------------------------------------------

    @staticmethod
    def _create_ibd_batch(
        tx: "ManagedTransaction", source: str,
        batch: list[IBDSegmentRow]
    ) -> None:
        rows = [
            {
                "sample_a": r.sample_a,
                "sample_b": r.sample_b,
                "chr": r.chr,
                "start": int(r.start),
                "end": int(r.end),
                "length_bp": int(r.end - r.start),
                "length_cM": float(r.length_cM) if r.length_cM is not None else None,
            }
            for r in batch
        ]
        tx.run(
            "UNWIND $rows AS r "
            "MATCH (a:Sample {sampleId: r.sample_a}), "
            "(b:Sample {sampleId: r.sample_b}) "
            "CREATE (a)-[:IBD_SEGMENT {"
            "  chr: r.chr, start: r.start, end: r.end, "
            "  length_bp: r.length_bp, length_cM: r.length_cM, "
            "  source: $source, created_at: datetime()"
            "}]->(b)",
            rows=rows, source=source)

    @staticmethod
    def _fail_if_source_exists(
        tx: "ManagedTransaction", source: str
    ) -> None:
        result = tx.run(
            "MATCH ()-[r:IBD_SEGMENT {source: $source}]->() "
            "RETURN count(r) AS c LIMIT 1",
            source=source)
        record = result.single()
        if record is not None and record["c"] > 0:
            raise RuntimeError(
                f"IBD segments already exist for source='{source}'. "
                "Pass replace=True to overwrite."
            )

    @staticmethod
    def _delete_source(
        tx: "ManagedTransaction", source: str
    ) -> int:
        result = tx.run(
            "MATCH ()-[r:IBD_SEGMENT {source: $source}]->() "
            "WITH count(r) AS n_deleted "
            "MATCH ()-[r:IBD_SEGMENT {source: $source}]->() "
            "DELETE r RETURN n_deleted",
            source=source)
        record = result.single()
        return int(record["n_deleted"]) if record is not None else 0

    @staticmethod
    def _list_sources(
        tx: "ManagedTransaction"
    ) -> list[IBDSourceSummary]:
        cypher = (
            "MATCH (a:Sample)-[r:IBD_SEGMENT]->(b:Sample) "
            "WITH r.source AS source, count(r) AS n_edges, "
            "count(DISTINCT [a.sampleId, b.sampleId]) AS n_pairs "
            "RETURN source, n_edges, n_pairs ORDER BY source"
        )
        out: list[IBDSourceSummary] = []
        for record in tx.run(cypher):
            out.append(IBDSourceSummary(
                source=record["source"],
                n_edges=int(record["n_edges"]),
                n_pairs=int(record["n_pairs"]),
            ))
        return out


# ---------------------------------------------------------------------------


def _chunked(iterable: Iterable[Any], size: int) -> Iterator[list[Any]]:
    batch: list[Any] = []
    for item in iterable:
        batch.append(item)
        if len(batch) >= size:
            yield batch
            batch = []
    if batch:
        yield batch
