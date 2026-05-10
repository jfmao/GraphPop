"""PLINK-compatible PED file exporter."""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, TextIO

from .reconstructor import PedigreeRow


class PedExporter:
    """Write a list of PedigreeRow to a PED file.

    PED format columns (whitespace-separated):
      FID  IID  FatherID  MotherID  Sex  Phenotype  [genotypes...]

    v1 emits the 6 mandatory columns only; downstream PLINK runs can
    join genotype columns via plink --merge.
    """

    @staticmethod
    def write(rows: Iterable[PedigreeRow], path: str | Path) -> int:
        rows_list = list(rows)
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as fh:
            return PedExporter._write_to(rows_list, fh)

    @staticmethod
    def to_string(rows: Iterable[PedigreeRow]) -> str:
        from io import StringIO
        buf = StringIO()
        PedExporter._write_to(list(rows), buf)
        return buf.getvalue()

    @staticmethod
    def _write_to(rows: list[PedigreeRow], fh: TextIO) -> int:
        for r in rows:
            fh.write(
                f"{r.family_id}\t{r.sample_id}\t{r.father_id}\t"
                f"{r.mother_id}\t{r.sex}\t{r.phenotype}\n"
            )
        return len(rows)
