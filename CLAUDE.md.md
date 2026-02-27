# GraphPop — Graph-Native Population Genomics Platform

## Project Identity
GraphPop is a graph-native compute engine for population genomic inference.
Core thesis: every population genomics computation is a traversal, aggregation,
or diffusion on a population variation graph. We replace VCF-centric workflows
with a unified Neo4j graph database.

## Current Phase: Phase 1 — Data Foundation + Summary Statistics
See `docs/GraphPop_Compiled_Holistic_Design.md` for full design.
See `docs/GraphPop_Revised_Implementation_Roadmap.html` for timeline.
Active work tracked in `tasks/todo.md`.

## Technology Stack
- **Graph DB:** Neo4j 5.x Community (WSL2), upgrade to Enterprise+Infinigraph later
- **Compute Engine:** Java 21+ plugin (graphpop-procedures.jar), SIMD via Vector API
- **Import Pipeline:** Python (cyvcf2, pandas) → CSV → neo4j-admin bulk import
- **Python API:** graphpop-py (neo4j-driver, plotting)
- **GNN/Embeddings:** PyTorch Geometric on RTX 4090 (Phase 3+)
- **Simulation:** msprime, SLiM (Phase 3B+)
- **Build:** Maven (Java), pyproject.toml (Python)

## Environment
- **OS:** WSL2 Ubuntu on Windows, project at `/mnt/e/GraphPop`
- **GPUs:** RTX 4090 (24GB VRAM) + Intel UHD 770
- **RAM:** 64 GB (plan upgrade to 128 GB)
- **Storage:** 1 TB + 4 TB NVMe
- **Neo4j data:** `/mnt/e/GraphPop/data/neo4j/`
- **Neo4j config:** page_cache=20GB, heap.max=16GB

## Architecture (must understand)

### Graph Schema v3 — Stratified Aggregation Model
- Variant nodes carry population allele count arrays (ac[], an[], af[] per pop)
- Sample→Variant via sparse CARRIES relationships (only non-ref genotypes)
- Variant→Variant via NEXT (positional chain) and LD (above threshold)
- Functional: Variant→Gene→Pathway→GOTerm

### Dual-Path Computation
- **FAST PATH:** Population stats from allele count arrays on Variant nodes.
  O(V×K) independent of sample size. Used for: π, θ_W, Tajima's D, F_ST, SFS, H_e, H_o
- **FULL PATH:** Individual stats via CARRIES traversal. Sparse, efficient for rare variants.
  Used for: LD (r², D'), iHS, XP-EHH, kinship, PCA/GRM, IBD, ROH

### Seven Layers (call downward only)
1. Data Foundation (schema, import, VectorOps SIMD)
2. Summary Statistics (SFS, diversity, divergence — fast + full paths)
3. Pairwise Statistics (kinship, IBS/IBD, ROH)
4. Graph Algorithms & Embeddings (community detection, GNN, PCA)
5. Inference Engines (demography, selection, recombination)
6. Simulation Integration (msprime/SLiM, ABC)
7. AI-Augmented Inference (GraphRAG, agents)

## Benchmark Datasets
1. **1000 Genomes chr22** — 2,504 samples, ~1.1M variants (human validation)
2. **Rice 3K chr1** — 3,000 samples, ~2M variants (crop validation)
3. **msprime simulated** — exact ground truth
Validate ALL statistics against vcftools, scikit-allel, PLINK, selscan, pixy.

## Key Naming Conventions
- Procedures: `graphpop.{module}.{method}` (e.g., `graphpop.diversity()`)
- Node IDs: `"chr1:12345:A:T"` for variants
- Run IDs: `"{method}_{pop}_{params_hash}_{timestamp}"`
- Package: `org.graphpop` (Java), `graphpop` (Python)

## Workflow Rules

### Plan First
- Enter plan mode for ANY non-trivial task (3+ steps or architectural decisions)
- Write plan to `tasks/todo.md` with checkable items
- Check in before starting implementation
- If something goes sideways, STOP and re-plan — don't keep pushing

### Subagent Strategy
- Use subagents to keep main context clean
- One task per subagent (e.g., "research LD computation algorithms")
- Offload: benchmarking, test writing, documentation, data exploration

### Self-Improvement
- After ANY correction: update `tasks/lessons.md`
- Write rules to prevent the same mistake
- Review lessons at session start

### Verification Standards
- Never mark a task complete without proving it works
- For statistics: validate against classical tools to <0.01% relative error
- For Java: unit tests with known-answer vectors
- For Python: pytest with fixtures from benchmark datasets
- Run `mvn test` for Java, `pytest` for Python before declaring done

### Code Quality
- Simplicity first: minimal code changes
- No temporary fixes — find root causes
- Java: follow Neo4j stored procedure conventions
- Python: type hints, docstrings, PEP 8
- All VectorOps have BOTH dense and sparse variants

### Domain Rules (Population Genomics)
- Homozygous reference = NO CARRIES edge (implicit, saves ~90% storage)
- CARRIES.gt: 1=HET, 2=HOM_ALT
- Missing data: -1 in dosage contexts; absent CARRIES = ref homozygote
- Allele count arrays indexed by population: pop_ids[], ac[], an[], af[]
- Ancestral state needed for unfolded SFS, Fay & Wu's H
- Always track polarization status (is_polarized boolean)

## Task Management
1. Write plan to `tasks/todo.md` with checkable items
2. Mark items complete as you go: `- [x] task`
3. Log architectural decisions in `tasks/decisions.md`
4. After corrections: update `tasks/lessons.md`
5. After completion: add review section to `tasks/todo.md`
