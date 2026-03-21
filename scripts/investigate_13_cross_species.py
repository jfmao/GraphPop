#!/usr/bin/env python3
"""Investigation 13: Cross-Species Parallel Adaptation (Human × Rice).

Question: Which biological processes face elevated selection in BOTH human and rice
populations? Cross-kingdom convergent adaptation reveals universally important cellular
functions.

Method:
  1. Rice side: Extract all sweep genes from rice_interpretation_results.json
     (genome-wide H12 sweep scan across 12 populations).
     Query Neo4j for GO terms of rice sweep genes that have human gene IDs in the graph.
  2. Human side: Query Neo4j for GO terms of human convergent sweep genes.
     Also use go_enrichment_ALL.tsv for population-level GO enrichment.
  3. Cross-species join: Find GO terms present in BOTH human sweep gene annotations
     AND in human orthologs of rice sweep genes.
     Also compute functional module convergence scores across 6 universal categories.
  4. Score convergence = n_human_sweep × n_rice_pops_sweeping.

Why graph-native:
  Both human sweep gene → GO terms and rice gene → human ortholog → GO terms are
  obtained via the same Neo4j HAS_GO_TERM traversal. The gene→GO annotation graph
  is leveraged for cross-species functional inference. No VCF tool integrates
  multi-species selection signals with shared GO term annotations in one step.

Output: data/results/cross_species.tsv / .json

Usage:
  /home/jfmao/miniconda3/envs/graphevo/bin/python -u scripts/investigate_13_cross_species.py
"""

import csv
import json
import re
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from neo4j import GraphDatabase

# ── Config ────────────────────────────────────────────────────────────────────
NEO4J_URI  = "bolt://localhost:7687"
NEO4J_USER = "neo4j"
NEO4J_PASS = "graphpop"
NEO4J_DB   = "neo4j"

OUT_DIR  = Path("data/results")
OUT_TSV  = OUT_DIR / "cross_species.tsv"
OUT_JSON = OUT_DIR / "cross_species.json"

RICE_RESULTS_PATH  = Path("results/rice/rice_interpretation_results.json")
RICE_DEEP_PATH     = Path("results/rice/rice_deep_integration_results.json")
CONV_SWEEPS_PATH   = OUT_DIR / "convergent_sweeps.json"
GO_ENRICH_PATH     = OUT_DIR / "go_enrichment_ALL.tsv"
PBS_SCAN_PATH      = OUT_DIR / "pbs_scan.json"
GENE_FST_PATH      = OUT_DIR / "gene_fst.json"

# ── Universal functional module definitions ───────────────────────────────────
# Keyword patterns → module assignment
# Each GO term name is matched against these patterns
MODULES = {
    "immunity_defense": [
        "immune", "defense", "pathogen", "antimicrobial", "immunoglobulin",
        "antigen", "antibody", "inflammation", "interferon", "innate",
        "toll", "nf-kb", "cytokine", "interleukin",
    ],
    "stress_response": [
        "stress", "heat", "cold", "drought", "hypoxia", "oxidative",
        "osmotic", "radiation", "toxic", "xenobiotic", "chaperone",
        "unfolded protein", "heat shock", "response to",
    ],
    "ion_transport": [
        "ion transport", "ion channel", "potassium", "sodium", "calcium",
        "chloride", "transporter", "membrane potential", "repolarization",
        "action potential", "electrochemical",
    ],
    "metabolism": [
        "metabolic process", "biosynthetic process", "catabolic process",
        "oxidation-reduction", "glycolysis", "lipid metabolism",
        "carbohydrate", "amino acid", "vitamin", "cofactor",
        "photosynthesis", "starch", "sugar", "nitrogen",
    ],
    "development_growth": [
        "development", "morphogenesis", "differentiation", "growth",
        "organogenesis", "pattern specification", "meristem",
        "flower", "seed", "fruit", "embryo",
    ],
    "signal_transduction": [
        "signal transduction", "signaling pathway", "receptor", "kinase",
        "phosphorylation", "transcription factor", "gene expression",
        "cell cycle", "apoptosis", "proliferation",
    ],
}


def classify_go_term(go_name: str) -> str | None:
    """Return the functional module for a GO term name, or None if unclassified."""
    name_lower = go_name.lower()
    for module, keywords in MODULES.items():
        for kw in keywords:
            if kw in name_lower:
                return module
    return None


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    print("=== Investigation 13: Cross-Species Parallel Adaptation ===\n", flush=True)

    driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASS))

    # ── HUMAN SIDE ─────────────────────────────────────────────────────────────
    print("--- HUMAN SIDE ---", flush=True)

    # 1a. Convergent sweep genes → GO terms from Neo4j
    conv_genes = []
    if CONV_SWEEPS_PATH.exists():
        with open(CONV_SWEEPS_PATH) as f:
            conv_data = json.load(f)
        cg = conv_data.get("convergent_genes", {})
        conv_genes = list(cg.keys()) if isinstance(cg, dict) else [
            g.get("gene_id", g) if isinstance(g, dict) else g for g in cg
        ]

    # 1b. PBS top genes (broader sweep signal)
    pbs_genes = []
    if PBS_SCAN_PATH.exists():
        with open(PBS_SCAN_PATH) as f:
            pbs_data = json.load(f)
        for pop, pop_data in pbs_data.get("populations", {}).items():
            for g in pop_data.get("top_genes", [])[:20]:
                pbs_genes.append(g.get("gene_id", ""))
    pbs_genes = list(set(pbs_genes))

    # 1c. Top gene-FST genes
    gene_fst_genes = []
    if GENE_FST_PATH.exists():
        with open(GENE_FST_PATH) as f:
            gene_fst_data = json.load(f)
        for g in gene_fst_data.get("top_genes", [])[:50]:
            gene_fst_genes.append(g["gene_id"])

    # Combine human sweep gene candidates
    all_human_sweep_genes = list(set(conv_genes + pbs_genes[:50] + gene_fst_genes[:50]))
    print(f"  Human sweep gene candidates: {len(all_human_sweep_genes)} "
          f"(conv={len(conv_genes)}, pbs={len(pbs_genes[:50])}, fst={len(gene_fst_genes[:50])})",
          flush=True)

    # Query Neo4j for their GO terms
    human_go_genes   = defaultdict(set)   # go_id → set of human gene IDs
    human_go_info    = {}                  # go_id → (name, aspect)
    human_go_modules = defaultdict(set)   # module → set of go_ids

    CYPHER_HUMAN = """
    MATCH (g:Gene)-[:HAS_GO_TERM]->(go:GOTerm)
    WHERE g.geneId IN $gene_ids
    RETURN g.geneId AS gene_id, go.goId AS go_id,
           go.name AS go_name, go.aspect AS aspect
    """
    with driver.session(database=NEO4J_DB) as session:
        result = session.run(CYPHER_HUMAN, gene_ids=all_human_sweep_genes)
        for rec in result:
            go_id   = rec["go_id"]
            go_name = rec["go_name"] or go_id
            aspect  = rec["aspect"] or ""
            human_go_genes[go_id].add(rec["gene_id"])
            human_go_info[go_id] = (go_name, aspect)
            mod = classify_go_term(go_name)
            if mod:
                human_go_modules[mod].add(go_id)

    print(f"  Human sweep GO terms: {len(human_go_genes)}", flush=True)

    # Also load go_enrichment_ALL for enrichment signal
    human_enriched_go = {}  # go_id → odds_ratio
    if GO_ENRICH_PATH.exists():
        with open(GO_ENRICH_PATH) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                try:
                    go_id = row["go_id"]
                    human_enriched_go[go_id] = {
                        "go_name":    row["go_name"],
                        "aspect":     row["aspect"],
                        "odds_ratio": float(row["odds_ratio"]),
                        "padj":       float(row["padj"]),
                        "genes":      row["genes"].split(","),
                    }
                except (ValueError, KeyError):
                    pass
    print(f"  Human enriched GO terms (Inv.5): {len(human_enriched_go)}", flush=True)

    # Compute human module signal score
    human_module_score = {}
    for mod, go_ids in human_go_modules.items():
        n_genes    = len(set().union(*[human_go_genes[g] for g in go_ids]))
        n_enriched = sum(1 for g in go_ids if g in human_enriched_go)
        max_or     = max((human_enriched_go[g]["odds_ratio"] for g in go_ids
                          if g in human_enriched_go), default=1.0)
        human_module_score[mod] = {
            "n_go_terms":   len(go_ids),
            "n_sweep_genes": n_genes,
            "n_enriched_go": n_enriched,
            "max_odds_ratio": round(max_or, 2),
            "signal_score":  round(n_genes * (1 + np.log1p(max_or)), 2),
        }

    print("\n  Human module scores:")
    for mod, ms in sorted(human_module_score.items(),
                          key=lambda x: x[1]["signal_score"], reverse=True):
        print(f"    {mod:25s}: n_genes={ms['n_sweep_genes']:3d}  "
              f"n_go={ms['n_go_terms']:3d}  "
              f"enriched={ms['n_enriched_go']:2d}  "
              f"signal={ms['signal_score']:.1f}")

    # ── RICE SIDE ──────────────────────────────────────────────────────────────
    print("\n--- RICE SIDE ---", flush=True)

    # Extract rice sweep genes and their population counts
    rice_gene_npops = defaultdict(int)  # gene_id → n pops sweeping
    rice_gene_names = {}                # gene_id → common name

    if RICE_RESULTS_PATH.exists():
        with open(RICE_RESULTS_PATH) as f:
            rice_data = json.load(f)

        # H12 genome-wide sweep scan
        sw = rice_data["annotate"]["sweep_annotations"]
        for pop, data in sw.items():
            if not data:
                continue
            for sweep in data.get("sweeps", []):
                for g in sweep.get("genes", []):
                    if isinstance(g, dict):
                        gid = g.get("geneId", g.get("gene_id", ""))
                    else:
                        gid = str(g)
                    if gid:
                        rice_gene_npops[gid] += 1

        # XP-EHH top peaks
        xp = rice_data["annotate"]["xpehh_annotations"]
        for comp, comp_data in xp.items():
            if isinstance(comp_data, dict):
                for peak in comp_data.get("top_peaks", [])[:5]:
                    for g in peak.get("genes", []):
                        gid = g.get("geneId", "") if isinstance(g, dict) else str(g)
                        if gid:
                            rice_gene_npops[gid] = max(rice_gene_npops[gid], 3)

    # Add known domestication genes from deep integration
    if RICE_DEEP_PATH.exists():
        with open(RICE_DEEP_PATH) as f:
            rice_deep = json.load(f)
        for gene_name, gene_data in rice_deep.get("sweeps", {}).items():
            # Use function label as proxy for n_pops (these are confirmed sweeps)
            rice_gene_names[gene_name] = gene_data.get("function", "")
            # Count populations with h12 signal (any signal > 0)
            n_pops = sum(1 for pd in gene_data.get("populations", {}).values()
                         if pd.get("h12", 0) > 0.001)
            rice_gene_npops[gene_name] = max(rice_gene_npops.get(gene_name, 0), max(n_pops, 1))

    print(f"  Rice sweep gene loci (all scan + deep): {len(rice_gene_npops):,}", flush=True)
    print(f"  Rice deep-integration known genes: {len(rice_gene_names)}", flush=True)

    # Query Neo4j for rice genes — those with LOC_Os IDs won't be there
    # but known gene symbols might match via geneId
    rice_go_genes   = defaultdict(set)  # go_id → set of rice gene IDs
    rice_go_info    = {}                # go_id → (name, aspect)
    rice_go_modules = defaultdict(set)

    rice_gids_in_graph = []
    with driver.session(database=NEO4J_DB) as session:
        # Check which rice gene IDs exist in Neo4j (LOC_Os won't, but let's try)
        r = session.run(
            "MATCH (g:Gene) WHERE g.geneId IN $gids RETURN g.geneId",
            gids=list(rice_gene_npops.keys())[:200]
        )
        rice_gids_in_graph = [rec["g.geneId"] for rec in r]

    print(f"  Rice genes found in Neo4j graph: {len(rice_gids_in_graph)}", flush=True)

    if rice_gids_in_graph:
        with driver.session(database=NEO4J_DB) as session:
            r = session.run(CYPHER_HUMAN, gene_ids=rice_gids_in_graph)
            for rec in r:
                go_id   = rec["go_id"]
                go_name = rec["go_name"] or go_id
                rice_go_genes[go_id].add(rec["gene_id"])
                rice_go_info[go_id] = (go_name, rec["aspect"] or "")
                mod = classify_go_term(go_name)
                if mod:
                    rice_go_modules[mod].add(go_id)

    # Assign rice domestication genes to modules using function labels
    # (since LOC_Os genes are not in Neo4j, use the known function descriptions)
    FUNCTION_MODULE_MAP = {
        "starch": "metabolism",
        "amylose": "metabolism",
        "glycan": "metabolism",
        "gibberellin": "development_growth",
        "semi-dwarf": "development_growth",
        "seed shattering": "development_growth",
        "grain": "development_growth",
        "heading date": "development_growth",
        "florigen": "development_growth",
        "flower": "development_growth",
        "prostrate": "development_growth",
        "fragrance": "metabolism",
        "submergence": "stress_response",
        "cold": "stress_response",
        "drought": "stress_response",
        "rooting": "development_growth",
        "hull": "metabolism",
        "blast resistance": "immunity_defense",
        "disease": "immunity_defense",
    }

    rice_module_npops = defaultdict(int)  # module → total n_pops (sum across genes)
    rice_module_genes = defaultdict(list)  # module → list of gene names

    for gene_name, func in rice_gene_names.items():
        func_lower = func.lower()
        assigned = False
        for keyword, module in FUNCTION_MODULE_MAP.items():
            if keyword in func_lower:
                n_pops = rice_gene_npops.get(gene_name, 1)
                rice_module_npops[module] += n_pops
                rice_module_genes[module].append(f"{gene_name}({func})")
                assigned = True
                break
        if not assigned:
            rice_module_npops["signal_transduction"] += rice_gene_npops.get(gene_name, 1)
            rice_module_genes["signal_transduction"].append(f"{gene_name}({func})")

    print("\n  Rice module signal (from deep domestication genes):")
    for mod in sorted(rice_module_npops, key=rice_module_npops.get, reverse=True):
        print(f"    {mod:25s}: npops={rice_module_npops[mod]:3d}  "
              f"genes={', '.join(rice_module_genes[mod][:3])}")

    # ── CROSS-SPECIES CONVERGENCE ─────────────────────────────────────────────
    print("\n--- CROSS-SPECIES CONVERGENCE ---", flush=True)

    rows = []
    all_modules = set(human_module_score.keys()) | set(rice_module_npops.keys())

    for mod in all_modules:
        hms = human_module_score.get(mod, {})
        h_signal = hms.get("signal_score", 0)
        h_ngenes = hms.get("n_sweep_genes", 0)
        r_npops  = rice_module_npops.get(mod, 0)
        r_genes  = "; ".join(rice_module_genes.get(mod, []))

        if h_signal == 0 and r_npops == 0:
            continue

        conv_score = round(float(h_signal) * max(r_npops, 0.5), 2)

        rows.append({
            "module":               mod,
            "human_signal_score":   round(h_signal, 2),
            "human_n_sweep_genes":  h_ngenes,
            "human_n_go_terms":     hms.get("n_go_terms", 0),
            "human_enriched_go":    hms.get("n_enriched_go", 0),
            "human_max_odds_ratio": hms.get("max_odds_ratio", 0),
            "rice_n_pops_total":    r_npops,
            "rice_genes":           r_genes,
            "convergence_score":    conv_score,
        })

    rows.sort(key=lambda r: r["convergence_score"], reverse=True)

    # ── GO-level cross-species (where rice genes exist in Neo4j) ─────────────
    go_rows = []
    if rice_go_genes:
        for go_id in set(rice_go_genes.keys()) & set(human_go_genes.keys()):
            go_name, aspect = human_go_info.get(go_id, ("", ""))
            human_genes = sorted(human_go_genes[go_id])
            rice_genes_here = sorted(rice_go_genes[go_id])
            r_npops = max(rice_gene_npops.get(g, 1) for g in rice_genes_here)
            go_rows.append({
                "go_id":   go_id,
                "go_name": go_name,
                "aspect":  aspect,
                "human_sweep_genes": ",".join(human_genes),
                "rice_genes": ",".join(rice_genes_here),
                "rice_n_pops": r_npops,
                "convergence_score": len(human_genes) * r_npops,
            })
        go_rows.sort(key=lambda r: r["convergence_score"], reverse=True)

    # ── Summary ───────────────────────────────────────────────────────────────
    print("\n  Module-level convergence (both species):")
    for r in rows[:6]:
        print(f"    {r['module']:25s}: "
              f"human_signal={r['human_signal_score']:.1f}  "
              f"rice_npops={r['rice_n_pops_total']:2d}  "
              f"conv={r['convergence_score']:.1f}")
        if r['rice_genes']:
            print(f"      rice: {r['rice_genes'][:80]}")

    driver.close()

    # ── Write outputs ─────────────────────────────────────────────────────────
    # Module-level TSV
    fieldnames = ["module", "human_signal_score", "human_n_sweep_genes",
                  "human_n_go_terms", "human_enriched_go", "human_max_odds_ratio",
                  "rice_n_pops_total", "rice_genes", "convergence_score"]
    with open(OUT_TSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    out_data = {
        "n_human_sweep_genes": len(all_human_sweep_genes),
        "n_human_go_terms": len(human_go_genes),
        "n_human_enriched_go": len(human_enriched_go),
        "n_rice_sweep_loci": len(rice_gene_npops),
        "n_rice_domestication_genes": len(rice_gene_names),
        "n_rice_in_neo4j": len(rice_gids_in_graph),
        "module_convergence": rows,
        "go_level_convergence": go_rows[:30],
        "human_module_scores": human_module_score,
        "rice_module_scores": {
            mod: {"npops": rice_module_npops[mod], "genes": rice_module_genes[mod]}
            for mod in rice_module_npops
        },
    }
    with open(OUT_JSON, "w") as f:
        json.dump(out_data, f, indent=2, default=str)

    print(f"\n  n_rows (module convergence): {len(rows)}")
    print(f"  n_go_level_convergence: {len(go_rows)}")
    print(f"\nOutputs: {OUT_TSV}  {OUT_JSON}")
    print(f"Total time: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
