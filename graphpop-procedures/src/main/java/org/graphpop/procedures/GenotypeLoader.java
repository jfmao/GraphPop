package org.graphpop.procedures;

import org.neo4j.graphdb.Label;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Transaction;

import java.util.*;

/**
 * Loads CARRIES edge data into dosage vectors for individual-level statistics.
 *
 * <p>Dosage encoding: 0 = ref/ref, 1 = het, 2 = hom_alt.
 * Non-carriers are implicitly 0 (sparse storage in CARRIES edges means only
 * non-ref genotypes are stored).</p>
 */
final class GenotypeLoader {

    private GenotypeLoader() {}

    /**
     * Build sample index: maps sampleId → array position for a given population.
     * The order is deterministic (sorted by sampleId).
     */
    static Map<String, Integer> buildSampleIndex(Transaction tx, String population) {
        var result = tx.execute(
                "MATCH (s:Sample)-[:IN_POPULATION]->(:Population {populationId: $pop}) " +
                "RETURN s.sampleId AS sid ORDER BY sid",
                Map.of("pop", population)
        );

        Map<String, Integer> index = new LinkedHashMap<>();
        try {
            while (result.hasNext()) {
                String sid = (String) result.next().get("sid");
                index.put(sid, index.size());
            }
        } finally {
            result.close();
        }
        return index;
    }

    /**
     * Load all carriers for a single variant into a dosage array.
     * dosage[sampleIndex] = gt value from CARRIES edge (1 for het, 2 for hom_alt).
     * Non-carriers remain 0 (ref/ref).
     */
    static double[] loadDosage(Transaction tx, String variantId,
                               Map<String, Integer> sampleIndex, int nSamples) {
        double[] dosage = new double[nSamples];

        var result = tx.execute(
                "MATCH (s:Sample)-[c:CARRIES]->(v:Variant {variantId: $varId}) " +
                "WHERE s.sampleId IN $sids " +
                "RETURN s.sampleId AS sid, c.gt AS gt",
                Map.of("varId", variantId, "sids", new ArrayList<>(sampleIndex.keySet()))
        );

        try {
            while (result.hasNext()) {
                Map<String, Object> row = result.next();
                String sid = (String) row.get("sid");
                Integer idx = sampleIndex.get(sid);
                if (idx != null) {
                    dosage[idx] = ((Number) row.get("gt")).doubleValue();
                }
            }
        } finally {
            result.close();
        }
        return dosage;
    }

    /**
     * Batch-load carriers for multiple variants in a single query.
     * Returns a map from variantId to dosage array.
     *
     * <p>This is much more efficient than calling loadDosage per variant,
     * as it issues a single Cypher query for the entire batch.</p>
     */
    static Map<String, double[]> loadDosageBatch(Transaction tx, List<String> variantIds,
                                                  Map<String, Integer> sampleIndex, int nSamples) {
        // Pre-allocate dosage arrays (all zeros = ref/ref)
        Map<String, double[]> dosageMap = new HashMap<>();
        for (String vid : variantIds) {
            dosageMap.put(vid, new double[nSamples]);
        }

        if (variantIds.isEmpty()) return dosageMap;

        var result = tx.execute(
                "MATCH (s:Sample)-[c:CARRIES]->(v:Variant) " +
                "WHERE v.variantId IN $varIds AND s.sampleId IN $sids " +
                "RETURN v.variantId AS vid, s.sampleId AS sid, c.gt AS gt",
                Map.of("varIds", variantIds,
                       "sids", new ArrayList<>(sampleIndex.keySet()))
        );

        try {
            while (result.hasNext()) {
                Map<String, Object> row = result.next();
                String vid = (String) row.get("vid");
                String sid = (String) row.get("sid");
                Integer idx = sampleIndex.get(sid);
                if (idx != null) {
                    double[] dosage = dosageMap.get(vid);
                    if (dosage != null) {
                        dosage[idx] = ((Number) row.get("gt")).doubleValue();
                    }
                }
            }
        } finally {
            result.close();
        }
        return dosageMap;
    }

    /**
     * Load carrier haplotypes for a single variant, organized by phase.
     *
     * <p>For phased data, each diploid sample contributes two haplotypes.
     * For a HET (gt=1), the phase property indicates which haplotype carries alt:
     * phase=0 means haplotype 0 carries alt, phase=1 means haplotype 1 carries alt.
     * For a HOM_ALT (gt=2), both haplotypes carry alt.</p>
     *
     * @return two int arrays [hap0, hap1] each of length nSamples.
     *         hap0[i] = allele on haplotype 0 for sample i (0 or 1).
     *         hap1[i] = allele on haplotype 1 for sample i (0 or 1).
     */
    static int[][] loadHaplotypes(Transaction tx, String variantId,
                                  Map<String, Integer> sampleIndex, int nSamples) {
        int[] hap0 = new int[nSamples];
        int[] hap1 = new int[nSamples];

        var result = tx.execute(
                "MATCH (s:Sample)-[c:CARRIES]->(v:Variant {variantId: $varId}) " +
                "WHERE s.sampleId IN $sids " +
                "RETURN s.sampleId AS sid, c.gt AS gt, c.phase AS phase",
                Map.of("varId", variantId, "sids", new ArrayList<>(sampleIndex.keySet()))
        );

        try {
            while (result.hasNext()) {
                Map<String, Object> row = result.next();
                String sid = (String) row.get("sid");
                Integer idx = sampleIndex.get(sid);
                if (idx == null) continue;

                int gt = ((Number) row.get("gt")).intValue();
                if (gt == 2) {
                    // Hom alt: both haplotypes carry alt
                    hap0[idx] = 1;
                    hap1[idx] = 1;
                } else if (gt == 1) {
                    // Het: phase tells us which haplotype has alt
                    Object phaseObj = row.get("phase");
                    int phase = phaseObj != null ? ((Number) phaseObj).intValue() : 0;
                    if (phase == 0) {
                        hap0[idx] = 1;
                    } else {
                        hap1[idx] = 1;
                    }
                }
            }
        } finally {
            result.close();
        }
        return new int[][]{hap0, hap1};
    }
}
