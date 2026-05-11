"""Competitor-tool wrappers. Populated incrementally in Step H1+.

Each wrapper emits a GraphPop-schema-compatible output so paired
comparison is mechanical.

Phase 1 wrappers (Paper 2 sim-only panels):
- plink_grm         (PLINK 2.0 `--make-grm-bin`)             [H1, shipped]
- king              (KING-robust `--kinship`)                [H2, shipped]
- tskit_branch_grm  (tskit `genetic_relatedness_matrix`)     [H3, shipped]
- egrm              (Fan, Mancuso & Chiang 2022 reference)   [H4, shipped]
- s_ldsc            (Finucane 2015 partitioned-h²)           [H5, pending]

Phase 2 wrappers (Paper 2 real-data panels, deferred):
- arg_rhe    (Zhu et al. 2025)
- threads    (Gunnarsson et al. 2024)
- admixture
- singer     (Deng et al. 2025)
"""
from .egrm import EgrmResult, EgrmRunner
from .king import KingResult, KingRunner, parse_kin_files
from .plink_grm import PlinkGrmResult, PlinkGrmRunner, parse_grm_bin
from .tskit_branch_grm import (
    TskitBranchGrmResult,
    TskitBranchGrmRunner,
)

__all__ = [
    "EgrmResult", "EgrmRunner",
    "KingResult", "KingRunner", "parse_kin_files",
    "PlinkGrmResult", "PlinkGrmRunner", "parse_grm_bin",
    "TskitBranchGrmResult", "TskitBranchGrmRunner",
]
