"""Competitor-tool wrappers. Populated in Step H1+ (one module per
external tool: plink_grm, king, tskit_branch_grm, egrm, s_ldsc,
admixture, threads, singer, arg_rhe). Each wrapper emits a
GraphPop-schema-compatible output so paired comparison is
mechanical.

Phase 1 wrappers (immediate): plink_grm, king, tskit_branch_grm,
egrm, s_ldsc.

Phase 2 wrappers (deferred): admixture, threads, singer, arg_rhe.
"""
