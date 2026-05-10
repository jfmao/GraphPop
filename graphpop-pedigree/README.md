# graphpop-pedigree

PRIMUS-style pedigree reconstruction from GraphPop's `:RELATIVE` labels
(M9, deferred from M5).

## What it does

Given a cohort with per-pair relationship labels written by
`graphpop.relate.classify` (parent_child / full_sibling /
second_degree / etc.) and family clustering by
`graphpop.relate.families`, this package reconstructs the most-likely
pedigree for each family and exports it as a PLINK-compatible PED
file or as `:PEDIGREE_PARENT_OF` edges back into the graph.

V1 covers:

- Parent-child trios (1 parent + 1 parent + 1 child via two
  parent_child edges + one full_sibling between non-parent samples
  inferred as siblings of the trio child)
- Full-sibling cohorts with or without parents
- 3-generation lineages when sex metadata
  (`:Sample.sex ∈ {"M","F"}`) is available

V1 defers: consanguinity / inbreeding inference, pedigrees beyond
3 generations, half-siblings.
