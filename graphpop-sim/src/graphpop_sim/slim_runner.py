"""SLiM orchestrator (M12.B).

Renders a curated SLiM template with parameter overrides, runs
SLiM as a subprocess, and returns the path to the resulting
.trees file. Loading the trees + ingest happens at the CLI layer
(same path as msprime).

Tests skip when the SLiM binary is not on PATH (CI-friendly).
"""
from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Any


CURATED_TEMPLATES = ("neutral", "sweep_genic", "bottleneck")


@dataclass
class SlimRun:
    """Result of a SLiM run."""

    template: str
    params: dict[str, Any]
    trees_path: Path
    return_code: int


class SlimRunner:
    """Render + run a SLiM template."""

    def __init__(self, slim_binary: str | None = None):
        self.slim_binary = slim_binary or shutil.which("slim") or "slim"

    @staticmethod
    def list_templates() -> list[str]:
        return list(CURATED_TEMPLATES)

    @classmethod
    def load_template(cls, template: str) -> str:
        if template not in CURATED_TEMPLATES:
            raise ValueError(
                f"unknown template '{template}'; "
                f"available: {CURATED_TEMPLATES}")
        # Read packaged template source.
        package_files = resources.files("graphpop_sim.slim_templates")
        src = (package_files / f"{template}.slim").read_text()
        return src

    def run(
        self,
        template: str,
        params: dict[str, Any],
        output_path: str | Path,
        verbose: bool = False,
    ) -> SlimRun:
        """Render template + invoke SLiM. Returns the run's metadata.

        Each parameter is passed to SLiM via -d 'key=value' so the
        template's `initialize()` and event blocks see it as a
        global variable.
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        params = {**params, "output_path": f"'{output_path}'"}

        # Run the on-disk template; SLiM resolves variables via -d.
        template_path = (resources.files("graphpop_sim.slim_templates")
                          / f"{template}.slim")
        cmd = [self.slim_binary]
        for k, v in params.items():
            cmd += ["-d", f"{k}={v}"]
        cmd.append(str(template_path))

        result = subprocess.run(
            cmd,
            capture_output=not verbose,
            text=True,
            check=False,
        )
        return SlimRun(
            template=template,
            params=params,
            trees_path=output_path,
            return_code=result.returncode,
        )

    @classmethod
    def is_available(cls) -> bool:
        """Return True if the SLiM binary is on PATH."""
        return shutil.which("slim") is not None
