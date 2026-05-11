# Fig 5c reference Python pipeline

The seven stages of a representative pop-gen workflow that
the Phase-1 Fig 5a Cypher template (one query, 20 lines)
replaces. These files are **committed reference artifacts**,
not executable from `fig5c_panels.py` — the driver only
**reads** them to count lines + files.

| File | Tool | What it does |
|---|---|---|
| `01_qc_and_filter.sh` | PLINK 2.0 | Variant QC + filter |
| `02_compute_grm.sh` | PLINK 2.0 | Genotype GRM `--make-grm-bin` |
| `03_compute_kinship.sh` | KING | Robust kinship `--kinship` |
| `04_admixture.sh` | ADMIXTURE | Global ancestry K-means |
| `05_local_ancestry.sh` | RFMix | Local-ancestry painting |
| `06_annotate_pathway.py` | VEP + Reactome | Variant → pathway mapping |
| `07_join_and_filter.py` | pandas | Cross-table join + filter |

To actually *run* this pipeline you would also need:

- 1000G chr22 VCF (`--vcf 1000G_chr22.vcf.gz`)
- 1000G reference panel + sample-map (`05_local_ancestry.sh`)
- Reactome pathway JSON (`06_annotate_pathway.py`)
- VEP-annotated variant table (`06_annotate_pathway.py`)

None of those are downloaded here; this is purely the
*pipeline-surface-area* demonstration that Paper 2 Fig 5c
contrasts against the GraphPop one-liner.

Runtime estimates from published per-tool benchmarks (see
`PLAN_fig5c.md` § "Runtime estimates"): the full pipeline
adds up to **≈ 106 minutes** on a 1000G-scale cohort. The
GraphPop equivalent (Phase-1 Fig 5a Cypher) is estimated at
**≈ 30 seconds** on the same data, once the cohort is
ingested.
