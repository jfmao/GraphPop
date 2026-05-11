"""Paper 2 Fig 5 driver — graph-native query plane (G5).

Per `PLAN_fig5.md`: the canonical 5a/5b/5c panels are Phase 2
(need 1000G + ancestry + Reactome ingest). This v1 ships the
two pieces that don't depend on real data:

- **5a** — literal Cypher query template ("find cryptic 2nd-
  degree relatives in EUR-ancestry segments enriched for
  cardiovascular pathway"). Demonstrates the composition of
  branch_grm_by_ancestry × restrict_to_pathway × time_window
  × relate.classify in a SINGLE Cypher CALL chain.
- **5d** — ecosystem map: every shipped CLI subcommand +
  competitor wrapper + paper2 driver + Java `@Procedure`
  name, introspected from the codebase.

No live Cypher execution here — the embedded test Neo4j only
has the 20-sample fixture; the cryptic-relative + ancestry +
pathway ingest is Phase 2.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List

# ---------------------------------------------------------------------------
# 5a — Cypher template
# ---------------------------------------------------------------------------

CYPHER_TEMPLATE = """\
// Paper 2 Fig 5a — graph-native composed kinship query.
// Find cryptic 2nd-degree relatives in EUR-ancestry segments
// enriched for cardiovascular pathway (Reactome P_cardio_signaling).
CALL graphpop.kinship.branch_grm_by_ancestry(
    $runId,
    {
        restrict_to_pathway: "P_cardio_signaling",
        ancestry: "EUR",
        time_window: [0, 50]              // last ~50 generations
    }
) YIELD sample_a, sample_b, phi_ancestry
WHERE phi_ancestry > 0.0625                // 2nd-degree threshold (1/16)
WITH sample_a, sample_b, phi_ancestry
CALL graphpop.relate.classify(
    sample_a, sample_b, phi_ancestry
) YIELD degree, confidence
WHERE degree = 2 AND confidence > 0.8
RETURN sample_a, sample_b, phi_ancestry, degree, confidence
ORDER BY phi_ancestry DESC
LIMIT 50;
"""

CYPHER_ANNOTATIONS = [
    ("graphpop.kinship.branch_grm_by_ancestry",
     "M4.4: ancestry-stratified branch GRM"),
    ("restrict_to_pathway", "G2: pathway predicate"),
    ("ancestry: \"EUR\"", "G3: local-ancestry partition"),
    ("time_window: [0, 50]",
     "branch-time predicate (last ~50 generations)"),
    ("phi_ancestry > 0.0625",
     "2nd-degree IBD threshold (1/16)"),
    ("graphpop.relate.classify", "M5: degree-of-relatedness"),
]


def build_cypher_template(out_dir: Path) -> Path:
    """Write `fig5a_cypher_template.txt` + a JSON of annotations."""
    out_dir.mkdir(parents=True, exist_ok=True)
    text_path = out_dir / "fig5a_cypher_template.txt"
    text_path.write_text(CYPHER_TEMPLATE)
    json_path = out_dir / "fig5a_cypher_annotations.json"
    json_path.write_text(json.dumps(
        [{"token": t, "explanation": e}
         for t, e in CYPHER_ANNOTATIONS],
        indent=2))
    return text_path


# ---------------------------------------------------------------------------
# 5d — Ecosystem introspection
# ---------------------------------------------------------------------------

_PROCEDURE_RE = re.compile(
    r'@Procedure\(\s*name\s*=\s*"([^"]+)"\s*,?\s*mode\s*=\s*Mode\.(READ|WRITE)\s*\)')


def find_java_procedures(
    java_src_root: Path,
) -> List[Dict[str, str]]:
    """Scan `*.java` for `@Procedure(name="...", mode=Mode.X)`.

    Returns one dict per procedure: `{name, mode, file, module}`.
    `module` = the sub-folder under `org/graphpop/procedures/`
    (e.g. `pairwise`, `arg`, `selection`).
    """
    if not java_src_root.exists():
        return []
    out: List[Dict[str, str]] = []
    for path in sorted(java_src_root.rglob("*.java")):
        text = path.read_text(errors="ignore")
        for m in _PROCEDURE_RE.finditer(text):
            name = m.group(1)
            mode = m.group(2)
            try:
                rel = path.relative_to(java_src_root)
                parts = rel.parts
                # parts is something like
                # ("org","graphpop","procedures","pairwise","BranchGrmProcedure.java")
                idx = parts.index("procedures")
                module = (parts[idx + 1]
                          if idx + 1 < len(parts) - 1
                          else "root")
            except (ValueError, IndexError):
                module = "unknown"
            out.append({
                "name": name, "mode": mode,
                "file": str(path),
                "module": module,
            })
    return out


_KNOWN_WRAPPERS = [
    ("plink_grm", "PlinkGrmRunner", "plink2 / plink"),
    ("king", "KingRunner", "king"),
    ("tskit_branch_grm", "TskitBranchGrmRunner", "tskit (Python)"),
    ("egrm", "EgrmRunner", "egrm (Python)"),
    ("s_ldsc", "SLdscRunner", "ldsc.py"),
]


def list_competitor_wrappers() -> List[Dict[str, str]]:
    """Enumerate the wrappers shipped under `graphpop_bench.competitors`.

    Returns the five Phase-1 wrappers + their runner class +
    binary/library dependency. Static list (mirrors the README);
    the alternative (importlib introspection) would force-import
    every wrapper just to list them, which is overkill here.
    """
    return [
        {"label": label, "runner": runner, "dependency": dep}
        for label, runner, dep in _KNOWN_WRAPPERS
    ]


def list_paper2_drivers() -> List[Dict[str, str]]:
    """Enumerate the per-figure drivers shipped under
    `graphpop_bench.paper2_drivers`."""
    here = Path(__file__).resolve().parent
    out: List[Dict[str, str]] = []
    for path in sorted(here.glob("*.py")):
        if path.name in {"__init__.py", "fig5_panels.py",
                         "fig5_figures.py"}:
            continue
        stem = path.stem
        if stem.startswith("fig"):
            # Map fig1de_panels -> figure label "Fig 1d/1e"
            label = stem.split("_")[0]
            out.append({
                "module": f"graphpop_bench.paper2_drivers.{stem}",
                "figure": label,
            })
    return out


def list_cli_commands() -> List[str]:
    """List the top-level CLI surface area for graphpop-bench."""
    return [
        "graphpop-bench profile",
        "graphpop-bench run plink_grm",
        "graphpop-bench run king",
        "graphpop-bench run tskit_branch_grm",
        "graphpop-bench run egrm",
        "graphpop-bench run s_ldsc",
        "python -m graphpop_bench.paper2_drivers.fig1de_panels",
        "python -m graphpop_bench.paper2_drivers.fig1de_figures",
        "python -m graphpop_bench.paper2_drivers.fig2_panels",
        "python -m graphpop_bench.paper2_drivers.fig2_figures",
        "python -m graphpop_bench.paper2_drivers.fig3c_panels",
        "python -m graphpop_bench.paper2_drivers.fig3ab_panels",
        "python -m graphpop_bench.paper2_drivers.fig3_figures",
        "python -m graphpop_bench.paper2_drivers.fig4_panels",
        "python -m graphpop_bench.paper2_drivers.fig4_figures",
    ]


DEFAULT_JAVA_SRC_ROOT = Path(
    "/mnt/data/GraphPop/graphpop-procedures/src/main/java")


def introspect_ecosystem(
    java_src_root: Path = DEFAULT_JAVA_SRC_ROOT,
) -> Dict[str, list]:
    """Build the full ecosystem dict for `fig5d.pdf`."""
    return {
        "java_procedures": find_java_procedures(java_src_root),
        "competitor_wrappers": list_competitor_wrappers(),
        "paper2_drivers": list_paper2_drivers(),
        "cli_commands": list_cli_commands(),
    }


def write_ecosystem_json(eco: dict, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    p = out_dir / "fig5d_ecosystem.json"
    p.write_text(json.dumps(eco, indent=2))
    return p


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

DEFAULT_OUTPUT = Path(
    "/mnt/data/GraphPop/paper/paper2_kinship_arg/benchmarks/fig5_out")


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    p.add_argument("--java-src-root", type=Path,
                   default=DEFAULT_JAVA_SRC_ROOT)
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_arg_parser().parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cypher_path = build_cypher_template(args.output_dir)
    eco = introspect_ecosystem(args.java_src_root)
    eco_path = write_ecosystem_json(eco, args.output_dir)
    print(f"[fig5a] cypher template → {cypher_path}", file=sys.stderr)
    print(f"[fig5d] ecosystem JSON → {eco_path}", file=sys.stderr)
    print(
        f"  procedures={len(eco['java_procedures'])} "
        f"wrappers={len(eco['competitor_wrappers'])} "
        f"drivers={len(eco['paper2_drivers'])} "
        f"cli={len(eco['cli_commands'])}",
        file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
