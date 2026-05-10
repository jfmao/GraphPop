# GraphPop Publication & Benchmarking Roadmap

**Status**: Draft v1 (2026-05-10). Companion to the personal planning
file at `~/.claude/plans/for-next-stage-of-peppy-noodle.md` (which
carries auto-mode execution details that don't belong in the repo).

---

## 1. Context

Through commit `3e04ef3` on `develop`, the GraphPop codebase covers
phases 1–6 plus the recombination v2 follow-ups (M13):

- **Phase 4 — pairwise statistics + ARG** (closed at M6):
  `kinship.king`, `kinship.branch_grm` and its conditional /
  posterior / ancestry-decomposed / PCA / HE-regression variants,
  `kinship.ibs`, ARG ingest, ancestry painting, IBD segments,
  `arg.tmrca`, `arg.coalescence_rate`, `arg.branch_diversity`,
  `arg.allele_age`.
- **M5**: `relate.classify`, `relate.families`.
- **Phase 5 — inference engines** (closed at M11):
  `demography.ne_trajectory`, `selection.allele_age_scan`,
  `selection.branch_outlier_scan`,
  `recombination.arg_breakpoints`, `recombination.ld_decay`.
- **Phase 4 deferred** (closed at M9–M10): `community.louvain`,
  `graphpop-pedigree` (PRIMUS-lite), `embedding.knn`,
  `embedding.cluster`, `graphpop-gnn` (GraphSAGE + InfoNCE
  training pipeline).
- **Phase 6 — simulation integration** (closed at M12):
  `graphpop-sim` (msprime + SLiM templates + rejection-ABC).
- **M13 recombination v2**: `recombination.ldhat_mcmc`,
  `recombination.hmm_smooth`, `recombination.stratified_ld_decay`,
  `recombination.hotspots`.

464 `mvn test` + 220 `pytest` green; 63 MCP tools; 70+ Cypher
procedures.

**Publication context**:

- **Paper 1** (platform): manuscript v5 on bioRxiv, formally
  submitted to a top journal, awaiting decision.
- **Target cadence**: 1 strong paper per year at top venues
  (Nature Methods, Genome Research, Molecular Biology and
  Evolution, PNAS).

This document defines Papers 2 → 5 (Y+1 → Y+4) and the shared
benchmarking infrastructure each one inherits.

---

## 2. Why publication plan before benchmarking?

Benchmarking is expensive (compute + writing) and benchmarks built
for the wrong target are wasted. The publication plan dictates:

- **Simulation design** — sweep recovery vs known-ρ recovery vs
  piecewise-Ne vs admixed cohort. Each paper has a different
  ground truth.
- **Competitor set** — top journals demand head-to-head against the
  strongest available baseline (PLINK / KING vs LDhat / pyrho vs
  selscan / SweepFinder2 vs PCA / ADMIXTURE).
- **Compute envelope** — 1k vs 50k vs biobank-scale 500k. The
  envelope is a paper-positioning choice.
- **Case study** — 1000G ancestry vs HapMap recombination vs UK
  Biobank kinship vs rice 3K agronomy. Each implies different
  data-use agreements + preprocessing.

Once each paper's benchmark spec is locked, the benchmarks are
mechanical execution.

---

## 3. Discipline rules

Two mandatory constraints apply to every paper-related task:

### 3.1 Plan-before-act

Before any non-trivial task (anything beyond a one-line edit),
produce a short written plan and commit it to disk **before**
executing. The plan is short — a section header, the goal, the
inputs, the outputs, the success criterion — but it must exist as
a diff'able artifact before the task starts.

### 3.2 Literature-survey-first for every paper

Before any manuscript-structuring or benchmark-design work for a
paper, conduct a thorough literature survey using bioRxiv,
PubMed, and web-search MCP tools. Output: a `literature_survey.md`
committed under that paper's folder. The survey defines:

- the most-recent state-of-the-art (no claiming novelty against
  5-year-old baselines);
- the precise gap our paper claims to close;
- the strongest competing methods to benchmark against;
- the mathematical results we'll cite vs derive ourselves.

These rules apply even when working under autonomous-execution
("go") approval.

---

## 4. Paper sequence

Each paper targets ONE sharp scientific question with clearly
demarcated methodological novelty. All papers ride on the same
graph-native substrate but each stands on its own scientific
contribution.

### Paper 2 (Y+1) — *Branch-GRM and conditional kinship from biobank-scale ARGs*

**Target**: Nature Methods.

**Research gap**: ARG-derived branch GRMs (Fan, Mancuso & Chiang
2022; Tang & Chiang 2025) demonstrate the principle on small
cohorts but the literature has three open gaps as of our initial
sketch (to be re-validated by the literature survey):

- **G1 (uncertainty)** — no published posterior over the branch
  GRM that propagates ARG-inference uncertainty into downstream
  kinship.
- **G2 (annotation co-residence)** — no published GRM
  conditioning on functional categories (pathway, consequence)
  propagated via the ARG.
- **G3 (local ancestry)** — no published GRM decomposed by
  inferred ancestry without external software pipelines.

**Headline novelty**:

- Closed-form `branch_grm` with three composable predicates
  (`restrict_to_pathway`, `mutation_filter`, `time_window`).
  Validated to < 10⁻⁶ rel-err vs `egrm.varGRM`.
- Posterior `branch_grm_posterior` via Welford aggregation across
  SINGER posterior samples (G1).
- `branch_grm_by_ancestry` decomposes by `:HAS_ANCESTRY`
  partition; exact: Σ_a B_ij^a == B_ij^uncond (G3).
- Algorithm V matvec form + Lanczos PCA + HE-regression
  heritability — all on the same packed ARG primitives.
- Biobank scale: never materialises the full n×n matrix.

**Benchmarking**:

| Axis | Plan |
|---|---|
| Simulation validation | msprime piecewise-Ne sims + tsdate-inferred ARGs; recover known kinship within 5 % rel-err |
| Comparison | PLINK GRM, KING-robust, REAP, RKM, GCTA on raw VCFs |
| Compute scale | 1k / 10k / 50k / 100k diploid samples; CPU-h + RAM peak |
| Conditional uniqueness | G2 path-specific GRM no competitor exposes |
| Application case study | 1000G + UK Biobank; recover known pedigrees via M5 |
| Math proofs | branch-GRM variance; partition exactness; Lanczos eigenpair convergence |

**Code status**: shipped (M4.1, M4.B, M4.4). Remaining: literature
survey, manuscript planning, benchmark execution, writing.

---

### Paper 3 (Y+2) — *Cross-validating ARG-derived and LD-derived recombination maps*

**Target**: Genome Research (alt: MBE).

**Research gap**: Recombination maps are inferred via either LD-
decay (LDhat, LDpop, pyrho) or ARG-derived (Hapne, iSMC, Speidel
et al.) methods. Each class has known failure modes — LD methods
fail in admixed cohorts; ARG methods fail when the inferred ARG
has biased breakpoints. **No published tool ships both side-by-
side**, and **no study has characterised the regions where they
systematically disagree**.

**Sharp question**: For a given cohort, where do the ARG-derived
ρ-map and the LD-derived ρ-map disagree, and what does the
disagreement tell us about admixture, selection, or ARG-inference
error?

**Headline novelty**:

- First side-by-side ARG-breakpoint + LD-decay ρ-maps in one
  query plane.
- Per-window posterior via LDhat-style MCMC (`ldhat_mcmc`).
- Per-window smoothing via pyrho-style HMM (`hmm_smooth`).
- Stratified ρ-maps by `:Sample.population` / `:Sample.sex`
  quantify admixture-induced discordance.
- Hotspot FDR via Benjamini-Hochberg.
- A novel **agreement score** flagging windows where ARG and LD
  diverge beyond a bootstrap-defined threshold.

**Benchmarking**:

| Axis | Plan |
|---|---|
| Simulation validation | msprime piecewise-constant ρ-maps; both estimators recover the truth within 25 % |
| Hotspot recovery | k synthetic hotspots; recall ≥ 0.8 at FDR ≤ 0.05 |
| Comparison | LDhat, LDpop, pyrho, iSMC |
| Compute scale | chr22 end-to-end; biobank-scale benchmark on UK Biobank chr20 |
| Disagreement diagnostic | Admixed cohort (CEU+YRI); ARG-LD disagreement tracks admixture-LD regions |
| Application case study | HapMap recombination hotspots + ARG-LD disagreement diagnostic |
| Math proofs | Hudson 1985 closed-form + bisection; Metropolis-Hastings detailed balance; HMM forward-backward correctness |

**Code status**: shipped (M11 + M13). Writing + benchmark-runs.

---

### Paper 4 (Y+3) — *Self-supervised GNN sample embeddings unify ancestry + relatedness summarisation*

**Target**: Nature Methods.

**Research gap**: PCA + ADMIXTURE remain the standard summaries
of cohort structure but require *two separate algorithms*, are
sensitive to LD filtering, and don't expose downstream queryable
embeddings. GNN-based pop-gen embeddings exist for variant-level
representations (e.g. PopVAE, DeepGenome) but not for sample-
level embeddings at biobank GPU scale, and none ship graph-
native KNN + clustering query procedures.

**Sharp question**: Can a self-supervised GraphSAGE-style GNN on
the `:Sample × :Variant` bipartite produce sample embeddings that
match PCA + ADMIXTURE on ancestry recovery, while also revealing
relatedness structure without an explicit kinship pass?

**Headline novelty**:

- 2-layer GraphSAGE + InfoNCE contrastive loss on co-carrier
  positive pairs; trained on RTX 4090.
- Embeddings persist as `:Sample.embedding` and feed
  `embedding.knn` / `embedding.cluster` query procedures.
- One pass recovers both ancestry (super-pop clustering) and
  relatedness (top-k neighbours align with M5 `:RELATIVE`
  edges).
- GPU training pipeline decoupled from Java query procedures.

**Benchmarking**:

| Axis | Plan |
|---|---|
| Ancestry recovery | 1000G super-population clustering; ARI ≥ 0.9 vs panel labels |
| Relatedness recovery | Synthetic trios + cousin pairs; top-5 neighbours ≥ 0.8 recall |
| Comparison | PCA (scikit-allel), ADMIXTURE, UMAP, popVAE |
| Compute scale | 1000G, HGDP, UK Biobank subset |
| GPU efficiency | RTX 4090 wall-clock per epoch + memory; A100 scaling |
| Robustness | LD-filtering / MAF / missingness ablations |
| Application case study | Cryptic relatedness in a synthetic admixed cohort |
| Math proofs | InfoNCE gradient; GraphSAGE Lipschitz bound; k-means++ variance |

**Code status**: shipped (M10). Writing + benchmark-runs.

---

### Paper 5 (Y+4) — *ARG-aware selection scans on biobank cohorts*

**Target**: Genome Research (alt: MBE).

**Research gap**: Classical selection scans (iHS, XP-EHH, nSL,
SweepFinder2, OmegaPlus) operate on haplotype patterns; none
use inferred ARGs directly. ARG-aware selection statistics have
been proposed (Speidel et al. 2019; Stern et al. 2019) but no
tool ships them as composable graph procedures with proper FDR
control.

**Sharp question**: Do ARG-derived selection signals (branch-
length outlier + allele-age z-score) recover known sweeps where
classical haplotype-based scans fail, or vice versa?

**Headline novelty**:

- `selection.allele_age_scan` — per-variant z-score of allele
  age conditional on derived-allele frequency.
- `selection.branch_outlier_scan` — Speidel-style per-window
  branch length restricted to focal samples.
- Cross-validation with classical methods (iHS, XP-EHH, nSL,
  SweepFinder2).
- BH-FDR on selection windows.

**Benchmarking**:

| Axis | Plan |
|---|---|
| Sweep recovery | msprime + SLiM sweep sims (genic / soft / polygenic); recover position within ±10 kb at z < −3 |
| FPR | Neutral sims; FPR ≤ 0.05 at q ≤ 0.05 |
| Comparison | selscan iHS, XP-EHH, nSL, SweepFinder2, OmegaPlus, RAiSD, iSAFE |
| Compute scale | 1000G chr2 (LCT), chr15 (SLC24A5) |
| Application case study | LCT / SLC24A5 / EDAR known sweeps |
| Math proofs | ARG conditional likelihood for allele age; branch-length null distribution; BH-FDR correctness |

**Code status**: shipped (M8). Writing + benchmark-runs.

---

## 5. Cross-paper infrastructure

A new package `graphpop-bench/` (sister to `graphpop-sim`,
`graphpop-gnn`, `graphpop-pedigree`) carries the shared
benchmarking machinery:

- **Simulation orchestrator**: extended from M12; adds piecewise-
  Ne, piecewise-ρ, sweep templates as YAML configs.
- **Comparison harness**: Python wrappers around competitor
  tools — one module per tool. Each wrapper emits a
  GraphPop-schema-compatible output so paired comparison is
  mechanical.
- **Compute profiling**: `/usr/bin/time -v` for CPU + RAM,
  `nvidia-smi` for GPU, ingest-stage timing for end-to-end.
- **Application data ETL**: 1000G chr22, HapMap, HGDP, rice 3K
  pipelines.

Per-paper folder layout under `paper/paper<N>_<topic>/`:

```
paper/paper<N>_<topic>/
├── literature_survey.md   # mandatory before any other work
├── manuscript_plan.md     # depends on literature_survey
├── benchmark_plan.md      # depends on manuscript_plan
├── benchmarks/            # scripts + intermediate outputs
├── figures/               # final figure-generation scripts
├── manuscript/            # LaTeX or markdown
└── data/                  # data agreements; never raw data
```

### Shared math-proof lemmas

A `paper/math/` folder hosts cross-paper lemma statements +
proofs:

- **Lemma A**: Hudson 1985 closed-form `E[r²|n, ρ·d]` + bisection
  convergence (Paper 3).
- **Lemma B**: Metropolis-Hastings detailed balance for the
  log-uniform-prior Gaussian-likelihood ρ-sampler (Paper 3).
- **Lemma C**: HMM forward-backward correctness for a banded-
  walk transition matrix on a discrete state space (Paper 3).
- **Lemma D**: Branch-GRM partition correctness:
  `Σ_ancestry B^a = B^unconditional` (Paper 2).
- **Lemma E**: InfoNCE gradient + convergence to the density-
  ratio estimator (Paper 4).
- **Lemma F**: Benjamini-Hochberg FDR control under positive-
  regression dependence (Papers 3 + 5).

---

## 6. Parallel-compute readiness

Audit of the 19 M6–M13 procedures relevant to all 4 papers:

- **Mode.READ** declared on all 18 read procedures; **Mode.WRITE**
  on `community.louvain` (persistence). Concurrent Cypher-level
  invocations across sessions are supported by default.
- **Within-procedure parallelism** (parallelStream on per-window
  loops, ForkJoinPool on per-pair sweeps) is NOT yet enabled.
  Hotspots are likely: LDhat MCMC, HMM smoothing, ARG-window
  procedures, GNN training (already GPU). v2 enhancement: profile
  → identify the slowest 3 → add `parallelStream()` after the
  data-load step. The Cypher transaction isn't thread-safe, so
  the parallel section must operate on pre-loaded in-memory data
  (the existing `LdPairLoader` / `ARGTraversal.load` pattern).

The papers themselves don't block on internal parallelism — the
benchmarks will reveal whether it's needed for paper-scale
runtimes. Add it surgically when the data demands it.

---

## 7. Strategic execution order (over 4 years)

1. **Paper 2 (kinship + ARG)** — kicks off when Paper 1's review
   decision returns. Strongest methodological story (3 closed
   gaps) and most-complete benchmark spec.
2. **Comparison harness for Paper 2** (PLINK, KING, REAP, RKM)
   seeds the cross-paper infrastructure.
3. **Paper 3 (recombination)** — starts after Paper 2 submission.
   Harness extended with LDhat / pyrho / iSMC wrappers.
4. **Paper 4 (GNN)** — starts as Paper 3 enters review. Harness
   extended with PCA / ADMIXTURE / popVAE wrappers.
5. **Paper 5 (selection)** — final paper. Harness extended with
   selscan / SweepFinder2 etc.

Each paper carries:

- one sharp research question,
- three or more methodological novelties,
- a head-to-head comparison against the strongest baselines,
- a compute-scale envelope demonstrated up to the largest
  publicly available cohort,
- a real-data application case study,
- a formal mathematical appendix.

Total: 4 top-tier papers over 4 years on the already-shipped
GraphPop substrate.

---

## 8. Out of scope

- **Phase 7 (GraphRAG + agents)** — explicitly deferred per
  current direction.
- **Paper 1 platform-paper revisions** — handled when the
  journal decision returns.
- **v2 follow-ups** to individual procedures (RJ-MCMC,
  ancestry-block stratification, continuous-state HMM, etc.)
  unless reviewer requests demand them.
- **Side projects** (GraphMana brief comm, educational slides)
  — tracked separately.
