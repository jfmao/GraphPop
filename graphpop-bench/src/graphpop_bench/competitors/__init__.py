"""Competitor-tool wrappers. Populated incrementally in Step H1+.

Each wrapper emits a GraphPop-schema-compatible output so paired
comparison is mechanical.

Phase 1 wrappers (Paper 2 sim-only panels):
- plink_grm  (PLINK 2.0 `--make-grm-bin`)            [H1, shipped]
- king       (KING-robust `--related`)               [H2, pending]
- tskit_branch_grm                                    [H3, pending]
- egrm       (Fan 2022 reference)                     [H4, pending]
- s_ldsc     (Finucane 2015 partitioned-h²)           [H5, pending]

Phase 2 wrappers (Paper 2 real-data panels, deferred):
- arg_rhe    (Zhu et al. 2025)
- threads    (Gunnarsson et al. 2024)
- admixture
- singer     (Deng et al. 2025)
"""
from .plink_grm import PlinkGrmResult, PlinkGrmRunner, parse_grm_bin

__all__ = ["PlinkGrmResult", "PlinkGrmRunner", "parse_grm_bin"]
