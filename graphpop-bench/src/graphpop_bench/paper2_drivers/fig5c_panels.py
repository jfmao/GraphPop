"""Paper 2 Fig 5c driver — Cypher-one-liner vs Python-pipeline LOC.

Static introspection of the seven committed reference pipeline
files at `fig5c_pipeline_reference/` plus the Phase-1 Fig 5a
Cypher template — counts lines + files + reports
literature-estimated runtimes for each stage. No execution.

Runtime numbers carry `runtime_source` provenance so the figure
caption can disambiguate measured-vs-estimated.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import List, Sequence

# ---------------------------------------------------------------------------
# Pipeline-stage registry (matches the committed reference dir)
# ---------------------------------------------------------------------------

PIPELINE_REFERENCE_DIR = (
    Path(__file__).resolve().parent / "fig5c_pipeline_reference")


@dataclass
class StageSpec:
    stage: str
    tool: str
    file: str
    runtime_s: float
    runtime_source: str


# Order matters — driver emits rows in this order so the figure
# can show a stacked / grouped layout.
PIPELINE_STAGES: List[StageSpec] = [
    StageSpec("qc_filter", "PLINK 2.0", "01_qc_and_filter.sh",
              120, "PLINK docs (1000G chr22)"),
    StageSpec("compute_grm", "PLINK 2.0", "02_compute_grm.sh",
              300, "PLINK 2.0 benchmarks"),
    StageSpec("compute_kinship", "KING", "03_compute_kinship.sh",
              480, "Manichaikul 2010 §3"),
    StageSpec("admixture", "ADMIXTURE", "04_admixture.sh",
              3600, "Alexander 2009 §3"),
    StageSpec("local_ancestry", "RFMix", "05_local_ancestry.sh",
              1500, "Maples 2013 §3"),
    StageSpec("pathway_annot", "VEP + Reactome",
              "06_annotate_pathway.py",
              300, "Ensembl VEP estimate"),
    StageSpec("join_filter", "pandas",
              "07_join_and_filter.py",
              60, "pandas in-memory join"),
]

CYPHER_STAGE = StageSpec(
    stage="composed_query",
    tool="GraphPop Cypher",
    file="../../../../paper/paper2_kinship_arg/benchmarks/"
         "fig5_out/fig5a_cypher_template.txt",
    runtime_s=30,
    runtime_source="Phase-1 Fig 4 N=1000 + procedure benchmarks",
)


# ---------------------------------------------------------------------------
# LOC counters
# ---------------------------------------------------------------------------

_COMMENT_PATTERNS = (
    re.compile(r"^\s*#"),               # shell + Python
    re.compile(r"^\s*//"),              # Cypher
)


def count_lines(path: Path) -> tuple[int, int]:
    """Return ``(total_lines, code_lines)``.

    A "code line" is non-blank and not a comment (under the
    shell `#`, Python `#`, or Cypher `//` conventions). Docstring
    bodies are counted as code (we don't attempt to parse Python
    triple-strings).
    """
    if not path.exists():
        raise FileNotFoundError(f"file not found: {path}")
    total = 0
    code = 0
    for raw in path.read_text().splitlines():
        total += 1
        stripped = raw.strip()
        if not stripped:
            continue
        if any(p.match(stripped) for p in _COMMENT_PATTERNS):
            continue
        code += 1
    return total, code


# ---------------------------------------------------------------------------
# Introspection
# ---------------------------------------------------------------------------

@dataclass
class PanelRow:
    source: str                # "pipeline" or "graphpop"
    stage: str
    tool: str
    file: str
    total_lines: int
    code_lines: int
    runtime_s: float
    runtime_source: str

    def as_csv_row(self) -> list:
        return [self.source, self.stage, self.tool, self.file,
                self.total_lines, self.code_lines,
                f"{self.runtime_s:.6g}", self.runtime_source]


def introspect_pipeline_reference(
    reference_dir: Path = PIPELINE_REFERENCE_DIR,
) -> List[PanelRow]:
    """Walk the committed reference dir; emit one PanelRow per stage."""
    rows: List[PanelRow] = []
    for spec in PIPELINE_STAGES:
        path = reference_dir / spec.file
        total, code = count_lines(path)
        rows.append(PanelRow(
            source="pipeline",
            stage=spec.stage, tool=spec.tool, file=spec.file,
            total_lines=total, code_lines=code,
            runtime_s=spec.runtime_s,
            runtime_source=spec.runtime_source,
        ))
    return rows


def cypher_template_row(
    template_path: Path,
) -> PanelRow:
    total, code = count_lines(template_path)
    return PanelRow(
        source="graphpop",
        stage=CYPHER_STAGE.stage,
        tool=CYPHER_STAGE.tool,
        file=str(template_path),
        total_lines=total, code_lines=code,
        runtime_s=CYPHER_STAGE.runtime_s,
        runtime_source=CYPHER_STAGE.runtime_source,
    )


def write_panel_csv(
    rows: Sequence[PanelRow], path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "source", "stage", "tool", "file",
            "total_lines", "code_lines",
            "runtime_s", "runtime_source",
        ])
        for r in rows:
            w.writerow(r.as_csv_row())


def write_summary_json(
    rows: Sequence[PanelRow], path: Path,
) -> None:
    """Compact summary that the figure-caption can quote."""
    pipeline = [r for r in rows if r.source == "pipeline"]
    cypher = [r for r in rows if r.source == "graphpop"]
    summary = {
        "pipeline": {
            "n_files": len(pipeline),
            "total_lines": sum(r.total_lines for r in pipeline),
            "code_lines": sum(r.code_lines for r in pipeline),
            "runtime_s": sum(r.runtime_s for r in pipeline),
        },
        "graphpop": {
            "n_files": len(cypher),
            "total_lines": sum(r.total_lines for r in cypher),
            "code_lines": sum(r.code_lines for r in cypher),
            "runtime_s": sum(r.runtime_s for r in cypher),
        },
        "notes": (
            "Runtime numbers are literature estimates at 1000G-"
            "chr22 scale. Empirical real-data runtime is a "
            "Phase 2-bis deliverable when KING / ADMIXTURE / "
            "RFMix are installed and 1000G is ingested."),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(summary, indent=2))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

DEFAULT_OUTPUT = Path(
    "/mnt/data/GraphPop/paper/paper2_kinship_arg/benchmarks/fig5_out")
DEFAULT_CYPHER_TEMPLATE = (
    DEFAULT_OUTPUT / "fig5a_cypher_template.txt")


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--reference-dir", type=Path,
                   default=PIPELINE_REFERENCE_DIR)
    p.add_argument("--cypher-template", type=Path,
                   default=DEFAULT_CYPHER_TEMPLATE)
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    pipeline_rows = introspect_pipeline_reference(args.reference_dir)
    cypher_row = cypher_template_row(args.cypher_template)
    rows = pipeline_rows + [cypher_row]

    csv_path = args.output_dir / "fig5c_panel_data.csv"
    summary_path = args.output_dir / "fig5c_summary.json"
    write_panel_csv(rows, csv_path)
    write_summary_json(rows, summary_path)
    print(f"[fig5c] panel CSV → {csv_path}", file=sys.stderr)
    print(f"[fig5c] summary JSON → {summary_path}", file=sys.stderr)

    # Headline summary to stderr.
    pipeline = [r for r in rows if r.source == "pipeline"]
    cypher = [r for r in rows if r.source == "graphpop"]
    print(
        f"  pipeline: {len(pipeline)} files, "
        f"{sum(r.code_lines for r in pipeline)} code lines, "
        f"~{sum(r.runtime_s for r in pipeline) / 60:.1f} min",
        file=sys.stderr)
    print(
        f"  graphpop: {len(cypher)} file,  "
        f"{sum(r.code_lines for r in cypher)} code lines, "
        f"~{sum(r.runtime_s for r in cypher):.0f} s",
        file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
