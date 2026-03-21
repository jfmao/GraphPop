# Novel Graph-Native Population Genomics Investigations

## Context

With GraphPop's Neo4j database (human 1000 Genomes) now containing:
- ~84M Variant nodes with ac/an/af[26 pops], XP-EHH peak scores (1,800 variants)
- 1.46M GenomicWindow nodes with π, θ_W, Tajima's D, PBS, H1/H12/sweep_type
- 3,202 Sample nodes with ROH per chromosome (roh_{chr}_froh/n_roh/total_kb)
- 91,973 Gene nodes → Pathway → GOTerm (functional annotation)
- CARRIES relationships (Sample → Variant, sparse, non-ref only)

## Five Novel Investigations

### 1. Pathway-Level Population Differentiation (highest priority)
**Question:** Which biological pathways are most divergently selected across human populations?
**Method:** Traverse Pathway→Gene→Variant, aggregate per-population af[] arrays,
compute mean Hudson F_ST per pathway, rank pathways.
**Why novel:** No VCF tool computes FST at pathway level natively. One Cypher query
replaces three tools + manual joins.
**Output:** Table of pathways ranked by population differentiation.

### 2. Convergent Positive Selection Across Independent Population Pairs
**Question:** Which genes/pathways show XP-EHH selection signals in African, European,
AND East Asian comparisons simultaneously?
**Method:** MATCH Variant nodes where xpehh_YRI_vs_CEU, xpehh_CHB_vs_CEU, AND
xpehh_YRI_vs_CHB are all non-null → traverse to Gene → aggregate by Pathway.
**Why novel:** Convergent adaptation detected as a graph pattern in one query.
**Output:** Gene/pathway list with convergent selection evidence.

### 3. ROH × Sweep Intersection — Inbreeding vs. Adaptation
**Question:** Do inbred samples (high FROH) preferentially carry adaptive alleles in
selection peak regions?
**Method:** Sample.roh_{chr}_froh > threshold → CARRIES → Variant with xpehh_* score →
Gene. Cross-reference by population.
**Why novel:** Connects individual inbreeding (sample-level) to population-level
adaptive evolution (variant-level) in one traversal.
**Output:** Per-population correlation of FROH vs. adaptive allele count.

### 4. Hard Sweep LD Network Topology
**Question:** Are hard sweep peak variants located in highly connected LD network hubs?
**Method:** Variant nodes with xpehh > threshold AND h12 > 0.9 → count LD edges →
compare to genome-wide LD degree distribution.
**Why novel:** Tests "haplotype surfing" theory via graph topology; never directly
measured this way in the literature.
**Output:** LD degree distribution for sweep peaks vs. neutral variants.

### 5. Multi-Hop Selection → GO Term Enrichment (4-hop traversal)
**Question:** Which GO terms are enriched in hard sweep genes per population?
**Method:** GenomicWindow(sweep_type='hard', population=P, h12>0.9)
→ Variant(SNP) → Gene → Pathway → GOTerm. Count distinct genes per GO term.
**Why novel:** 4-hop traversal across the full schema in one Cypher query.
VCF-centric workflows need 3+ separate tools.
**Output:** GO enrichment table per population, comparable across all 26 pops.

## Recommended First Paper

**Title candidate:** "Graph-native integration reveals pathway-level convergent positive
selection in human populations"
**Core result:** Investigations 1 + 2 combined: pathway FST ranking cross-validated
by XP-EHH convergence signals. No existing pipeline can produce this in one step.

## Prerequisites / Caveats
- Variant→Gene (IN_GENE) edge completeness needs verification before pathway aggregation
- XP-EHH write-back covers only top 1,800 peaks (not genome-wide) — investigations
  2–4 are limited to top signals; genome-wide XP-EHH scores would strengthen #4
- Pathway-level FST requires correction for pathway size and LD between member variants
- Rice data needed in Neo4j for cross-species comparative investigation (future)
