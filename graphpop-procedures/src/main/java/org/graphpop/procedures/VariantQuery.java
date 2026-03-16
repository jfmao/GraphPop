package org.graphpop.procedures;

import java.util.HashMap;
import java.util.Map;

/**
 * Shared Cypher query builder for variant retrieval with optional
 * consequence and pathway filters.
 *
 * <p>Extracted from DiversityProcedure's private buildVariantQuery method
 * so all FAST PATH procedures can share the same query construction logic.</p>
 *
 * <p>When {@code consequence} or {@code pathway} options are provided, the query
 * joins through HAS_CONSEQUENCE and/or IN_PATHWAY edges. Otherwise, a simple
 * Variant-only query is returned.</p>
 *
 * <p>Consequence and pathway values are passed as Cypher parameters ($consequence,
 * $pathway) to prevent injection. Use {@link #params(Map, Map)} to merge them
 * into the base parameter map before executing the query.</p>
 */
final class VariantQuery {

    private VariantQuery() {}

    /**
     * Build a Cypher query for fetching variants in a region, with optional
     * consequence and pathway filters from the options map.
     *
     * <p>Consequence and pathway are referenced as $consequence and $pathway
     * parameters. Use {@link #params(Map, Map)} to build the parameter map.</p>
     *
     * @param options user-supplied options (may contain "consequence" and/or "pathway")
     * @return Cypher query string with parameters $chr, $start, $end
     *         (and optionally $consequence, $pathway)
     */
    static String build(Map<String, Object> options) {
        String consequence = options != null ? (String) options.get("consequence") : null;
        String pathway = options != null ? (String) options.get("pathway") : null;
        return build(consequence, pathway);
    }

    /**
     * Build a Cypher query for fetching variants in a region.
     *
     * @param consequence VEP consequence filter (e.g., "missense_variant"), or null
     * @param pathway     Reactome pathway name filter, or null
     * @return Cypher query string with parameters $chr, $start, $end
     *         (and optionally $consequence, $pathway)
     */
    static String build(String consequence, String pathway) {
        StringBuilder sb = new StringBuilder();

        if (pathway != null && consequence != null) {
            sb.append("MATCH (v:Variant)-[c:HAS_CONSEQUENCE]->(g:Gene)-[:IN_PATHWAY]->(p:Pathway) ");
            sb.append("WHERE v.chr = $chr AND v.pos >= $start AND v.pos <= $end ");
            sb.append("AND c.consequence = $consequence ");
            sb.append("AND p.name = $pathway ");
        } else if (pathway != null) {
            sb.append("MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene)-[:IN_PATHWAY]->(p:Pathway) ");
            sb.append("WHERE v.chr = $chr AND v.pos >= $start AND v.pos <= $end ");
            sb.append("AND p.name = $pathway ");
        } else if (consequence != null) {
            sb.append("MATCH (v:Variant)-[c:HAS_CONSEQUENCE]->(g:Gene) ");
            sb.append("WHERE v.chr = $chr AND v.pos >= $start AND v.pos <= $end ");
            sb.append("AND c.consequence = $consequence ");
        } else {
            sb.append("MATCH (v:Variant) ");
            sb.append("WHERE v.chr = $chr AND v.pos >= $start AND v.pos <= $end ");
        }

        sb.append("RETURN DISTINCT v");
        return sb.toString();
    }

    /**
     * Build a Cypher query for fetching all variants on a chromosome (no region bounds),
     * sorted by position. Used by GenomeScanProcedure.
     */
    static String buildChromosome(Map<String, Object> options) {
        String consequence = options != null ? (String) options.get("consequence") : null;
        String pathway = options != null ? (String) options.get("pathway") : null;

        StringBuilder sb = new StringBuilder();

        if (pathway != null && consequence != null) {
            sb.append("MATCH (v:Variant)-[c:HAS_CONSEQUENCE]->(g:Gene)-[:IN_PATHWAY]->(p:Pathway) ");
            sb.append("WHERE v.chr = $chr ");
            sb.append("AND c.consequence = $consequence ");
            sb.append("AND p.name = $pathway ");
        } else if (pathway != null) {
            sb.append("MATCH (v:Variant)-[:HAS_CONSEQUENCE]->(g:Gene)-[:IN_PATHWAY]->(p:Pathway) ");
            sb.append("WHERE v.chr = $chr ");
            sb.append("AND p.name = $pathway ");
        } else if (consequence != null) {
            sb.append("MATCH (v:Variant)-[c:HAS_CONSEQUENCE]->(g:Gene) ");
            sb.append("WHERE v.chr = $chr ");
            sb.append("AND c.consequence = $consequence ");
        } else {
            sb.append("MATCH (v:Variant) ");
            sb.append("WHERE v.chr = $chr ");
        }

        sb.append("RETURN DISTINCT v ORDER BY v.pos");
        return sb.toString();
    }

    /**
     * Merge consequence/pathway from the user options into a base parameter map.
     * Returns the base map unchanged if no consequence/pathway is set.
     *
     * @param options    user-supplied options (may contain "consequence" and/or "pathway")
     * @param baseParams base Cypher parameters (e.g., chr, start, end)
     * @return merged parameter map (may be the same instance if no merging needed)
     */
    static Map<String, Object> params(Map<String, Object> options,
                                       Map<String, Object> baseParams) {
        String consequence = options != null ? (String) options.get("consequence") : null;
        String pathway = options != null ? (String) options.get("pathway") : null;

        if (consequence == null && pathway == null) return baseParams;

        Map<String, Object> merged = new HashMap<>(baseParams);
        if (consequence != null) merged.put("consequence", consequence);
        if (pathway != null) merged.put("pathway", pathway);
        return merged;
    }
}
