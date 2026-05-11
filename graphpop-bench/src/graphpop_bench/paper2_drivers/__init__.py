"""Per-figure driver scripts for Paper 2 (kinship + ARG).

Each module orchestrates one figure's data pipeline: runs the
relevant competitor wrappers from `graphpop_bench.competitors`,
extracts GraphPop's procedure output via Maven-driven Java test
dumps, joins by sample pair, and emits tidy panel CSVs.

Figure-rendering modules (matplotlib) are co-located so the
panel-CSV → PDF pipeline is reproducible from this package.

Shipped:
- fig1de_panels   — Fig 1d/1e branch-GRM rel-err panel data
- fig1de_figures  — Fig 1d/1e PDF rendering

Pending:
- fig2 (posterior branch GRM)
- fig3 (annotation-conditional GRMs)
- fig4 (biobank-scale validation)
- fig5 (graph-native query plane)

The submodules are NOT auto-imported from this `__init__` so
that ``python -m graphpop_bench.paper2_drivers.fig1de_panels``
runs cleanly without the "module-already-in-sys.modules"
RuntimeWarning. Import the helpers you need explicitly:

    from graphpop_bench.paper2_drivers.fig1de_panels import (
        join_pairs, load_pair_tsv,
    )
"""
