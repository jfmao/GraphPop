package org.graphpop.procedures;

import org.neo4j.graphdb.*;

import java.util.*;

/**
 * Dense in-memory haplotype matrix for a genomic region and population.
 *
 * <p>Loads all variant × haplotype data in a single pass via the graph API,
 * enabling parallel computation with zero DB access. The matrix layout is:
 * {@code haplotypes[variantIndex][haplotypeIndex]} where
 * {@code haplotypeIndex = 2 * sampleIndex + phase} (0=ref, 1=alt).</p>
 *
 * <p>Memory usage: nVariants × 2 × nSamples bytes. For EUR on chr22:
 * ~1M × 1006 ≈ 1 GB (fits in 4 GB heap).</p>
 */
final class HaplotypeMatrix {

    private static final RelationshipType CARRIES_REL = RelationshipType.withName("CARRIES");

    /** Variant IDs in position order. */
    final String[] variantIds;

    /** Genomic positions in sorted order. */
    final long[] positions;

    /** Node references for writing results back. */
    final Node[] nodes;

    /** Allele frequency for the target population at each variant. */
    final double[] afs;

    /** Allele counts for the target population. */
    final int[] acs;

    /** Allele numbers for the target population. */
    final int[] ans;

    /**
     * Dense haplotype matrix: {@code haplotypes[variantIdx][hapIdx]}.
     * hapIdx = 2 * sampleIdx + phase. Values: 0 (ref) or 1 (alt).
     */
    final byte[][] haplotypes;

    final int nVariants;
    final int nHaplotypes;
    final int nSamples;

    /** Sample index used during loading (sampleId → array position). */
    final Map<String, Integer> sampleIndex;

    /** Lazy dosage cache: populated on first dosage(variantIdx) call. */
    private volatile double[][] dosageCache;

    private HaplotypeMatrix(String[] variantIds, long[] positions, Node[] nodes,
                            double[] afs, int[] acs, int[] ans,
                            byte[][] haplotypes, int nSamples,
                            Map<String, Integer> sampleIndex) {
        this.variantIds = variantIds;
        this.positions = positions;
        this.nodes = nodes;
        this.afs = afs;
        this.acs = acs;
        this.ans = ans;
        this.haplotypes = haplotypes;
        this.nVariants = variantIds.length;
        this.nHaplotypes = nSamples * 2;
        this.nSamples = nSamples;
        this.sampleIndex = sampleIndex;
    }

    /**
     * Build a HaplotypeMatrix from raw arrays for unit testing (no DB required).
     */
    static HaplotypeMatrix forTest(String[] variantIds, long[] positions,
                                   double[] afs, int[] acs, int[] ans,
                                   byte[][] haplotypes, int nSamples) {
        return new HaplotypeMatrix(variantIds, positions, null, afs, acs, ans,
                haplotypes, nSamples, Map.of());
    }

    /**
     * Get the dosage vector for a variant (0/1/2 per sample).
     * dosage[s] = haplotypes[v][2*s] + haplotypes[v][2*s+1].
     *
     * <p>Uses a lazy cache: first call computes and stores the result,
     * subsequent calls return the cached array.</p>
     */
    double[] dosage(int variantIdx) {
        double[][] cache = dosageCache;
        if (cache != null && cache[variantIdx] != null) {
            return cache[variantIdx];
        }
        double[] dos = computeDosage(variantIdx);
        if (cache == null) {
            synchronized (this) {
                if (dosageCache == null) {
                    dosageCache = new double[nVariants][];
                }
                cache = dosageCache;
            }
        }
        cache[variantIdx] = dos;
        return dos;
    }

    private double[] computeDosage(int variantIdx) {
        double[] dos = new double[nSamples];
        byte[] h = haplotypes[variantIdx];
        for (int s = 0; s < nSamples; s++) {
            dos[s] = h[2 * s] + h[2 * s + 1];
        }
        return dos;
    }

    /**
     * Get a single haplotype value.
     */
    int haplotype(int variantIdx, int sampleIdx, int phase) {
        return haplotypes[variantIdx][2 * sampleIdx + phase];
    }

    /**
     * Load a haplotype matrix for a chromosome region using a pre-built sample index.
     * This overload supports custom sample subsets — the caller builds the sample index
     * from an explicit list of sample IDs.
     *
     * @param tx          active transaction
     * @param chr         chromosome
     * @param sampleIndex pre-built sample index (sampleId → array position)
     * @param start       region start (use {@code null} for whole chromosome)
     * @param end         region end (use {@code null} for whole chromosome)
     * @return loaded matrix, or null if no variants/samples found
     */
    static HaplotypeMatrix load(Transaction tx, String chr,
                                Map<String, Integer> sampleIndex,
                                Long start, Long end) {
        return loadInternal(tx, chr, sampleIndex, start, end);
    }

    /**
     * Load a haplotype matrix for a chromosome region and population.
     *
     * @param tx    active transaction
     * @param chr   chromosome
     * @param pop   population ID
     * @param start region start (use {@code null} for whole chromosome)
     * @param end   region end (use {@code null} for whole chromosome)
     * @return loaded matrix, or null if no variants/samples found
     */
    static HaplotypeMatrix load(Transaction tx, String chr, String pop,
                                Long start, Long end) {
        // 1. Build sample index for the population
        Map<String, Integer> sampleIndex = GenotypeLoader.buildSampleIndex(tx, pop);
        return loadInternal(tx, chr, sampleIndex, start, end);
    }

    private static HaplotypeMatrix loadInternal(Transaction tx, String chr,
                                                Map<String, Integer> sampleIndex,
                                                Long start, Long end) {
        int nSamples = sampleIndex.size();
        if (nSamples == 0) return null;
        int nHaplotypes = nSamples * 2;

        // 2. Query all variants in region sorted by position
        String query;
        Map<String, Object> params;
        if (start != null && end != null) {
            query = "MATCH (v:Variant) WHERE v.chr = $chr AND v.pos >= $start AND v.pos <= $end " +
                    "RETURN v ORDER BY v.pos";
            params = Map.of("chr", chr, "start", start, "end", end);
        } else {
            query = "MATCH (v:Variant) WHERE v.chr = $chr RETURN v ORDER BY v.pos";
            params = Map.of("chr", chr);
        }

        var result = tx.execute(query, params);

        // 3. Collect variant metadata
        List<String> vidList = new ArrayList<>();
        List<Long> posList = new ArrayList<>();
        List<Node> nodeList = new ArrayList<>();

        try {
            while (result.hasNext()) {
                Map<String, Object> row = result.next();
                Node node = (Node) row.get("v");

                vidList.add((String) node.getProperty("variantId"));
                posList.add(((Number) node.getProperty("pos")).longValue());
                nodeList.add(node);
            }
        } finally {
            result.close();
        }

        int nVariants = vidList.size();
        if (nVariants == 0) return null;

        // 4. Convert to arrays
        String[] variantIds = vidList.toArray(new String[0]);
        long[] positions = new long[nVariants];
        Node[] nodes = nodeList.toArray(new Node[0]);
        for (int i = 0; i < nVariants; i++) {
            positions[i] = posList.get(i);
        }

        // 5. Load haplotype data via graph API (one pass over CARRIES edges)
        //    Compute ac/an/af from the loaded haplotypes
        byte[][] haplotypes = new byte[nVariants][nHaplotypes];
        int[] acs = new int[nVariants];
        int[] ans = new int[nVariants];
        double[] afs = new double[nVariants];

        for (int v = 0; v < nVariants; v++) {
            Node varNode = nodes[v];
            byte[] hapRow = haplotypes[v];

            for (Relationship rel : varNode.getRelationships(Direction.INCOMING, CARRIES_REL)) {
                Node sampleNode = rel.getStartNode();
                String sid = (String) sampleNode.getProperty("sampleId");
                Integer sIdx = sampleIndex.get(sid);
                if (sIdx == null) continue; // not in target sample set

                int gt = ((Number) rel.getProperty("gt")).intValue();
                if (gt == 2) {
                    hapRow[2 * sIdx] = 1;
                    hapRow[2 * sIdx + 1] = 1;
                } else if (gt == 1) {
                    Object phaseObj = rel.getProperty("phase", null);
                    int phase = phaseObj != null ? ((Number) phaseObj).intValue() : 0;
                    hapRow[2 * sIdx + phase] = 1;
                }
            }

            // Compute AC from haplotypes
            int ac = 0;
            for (int h = 0; h < nHaplotypes; h++) ac += hapRow[h];
            acs[v] = ac;
            ans[v] = nHaplotypes;
            afs[v] = (double) ac / nHaplotypes;
        }

        return new HaplotypeMatrix(variantIds, positions, nodes, afs, acs, ans,
                haplotypes, nSamples, sampleIndex);
    }

    /**
     * Load two population matrices from the same variant set.
     * More efficient than loading separately since variants are queried once.
     */
    static HaplotypeMatrix[] loadPair(Transaction tx, String chr,
                                      String pop1, String pop2,
                                      Long start, Long end) {
        Map<String, Integer> sampleIndex1 = GenotypeLoader.buildSampleIndex(tx, pop1);
        Map<String, Integer> sampleIndex2 = GenotypeLoader.buildSampleIndex(tx, pop2);
        int nSamples1 = sampleIndex1.size();
        int nSamples2 = sampleIndex2.size();
        if (nSamples1 == 0 || nSamples2 == 0) return null;

        // Query variants once
        String query;
        Map<String, Object> params;
        if (start != null && end != null) {
            query = "MATCH (v:Variant) WHERE v.chr = $chr AND v.pos >= $start AND v.pos <= $end " +
                    "RETURN v ORDER BY v.pos";
            params = Map.of("chr", chr, "start", start, "end", end);
        } else {
            query = "MATCH (v:Variant) WHERE v.chr = $chr RETURN v ORDER BY v.pos";
            params = Map.of("chr", chr);
        }

        var result = tx.execute(query, params);

        List<String> vidList = new ArrayList<>();
        List<Long> posList = new ArrayList<>();
        List<Node> nodeList = new ArrayList<>();
        List<double[]> afArrList = new ArrayList<>();
        List<int[]> acArrList = new ArrayList<>();
        List<int[]> anArrList = new ArrayList<>();

        PopulationContext ctx1 = null;
        PopulationContext ctx2 = null;

        try {
            while (result.hasNext()) {
                Map<String, Object> row = result.next();
                Node node = (Node) row.get("v");

                if (ctx1 == null) {
                    ctx1 = PopulationContext.resolve(tx, node, pop1);
                    ctx2 = PopulationContext.resolve(tx, node, pop2);
                }

                vidList.add((String) node.getProperty("variantId"));
                posList.add(((Number) node.getProperty("pos")).longValue());
                nodeList.add(node);
                afArrList.add(ArrayUtil.toDoubleArray(node.getProperty("af")));
                acArrList.add(ArrayUtil.toIntArray(node.getProperty("ac")));
                anArrList.add(ArrayUtil.toIntArray(node.getProperty("an")));
            }
        } finally {
            result.close();
        }

        int nVariants = vidList.size();
        if (nVariants == 0) return null;

        String[] variantIds = vidList.toArray(new String[0]);
        long[] positions = new long[nVariants];
        Node[] nodes = nodeList.toArray(new Node[0]);
        double[] afs1 = new double[nVariants], afs2 = new double[nVariants];
        int[] acs1 = new int[nVariants], acs2 = new int[nVariants];
        int[] ans1 = new int[nVariants], ans2 = new int[nVariants];

        for (int i = 0; i < nVariants; i++) {
            positions[i] = posList.get(i);
            double[] afArr = afArrList.get(i);
            int[] acArr = acArrList.get(i);
            int[] anArr = anArrList.get(i);
            afs1[i] = afArr[ctx1.index]; afs2[i] = afArr[ctx2.index];
            acs1[i] = acArr[ctx1.index]; acs2[i] = acArr[ctx2.index];
            ans1[i] = anArr[ctx1.index]; ans2[i] = anArr[ctx2.index];
        }

        // Load haplotypes for both populations from the same variant nodes
        byte[][] hap1 = new byte[nVariants][nSamples1 * 2];
        byte[][] hap2 = new byte[nVariants][nSamples2 * 2];

        for (int v = 0; v < nVariants; v++) {
            Node varNode = nodes[v];
            for (Relationship rel : varNode.getRelationships(Direction.INCOMING, CARRIES_REL)) {
                Node sampleNode = rel.getStartNode();
                String sid = (String) sampleNode.getProperty("sampleId");

                int gt = ((Number) rel.getProperty("gt")).intValue();
                Object phaseObj = rel.getProperty("phase", null);
                int phase = phaseObj != null ? ((Number) phaseObj).intValue() : 0;

                Integer sIdx1 = sampleIndex1.get(sid);
                if (sIdx1 != null) {
                    if (gt == 2) {
                        hap1[v][2 * sIdx1] = 1;
                        hap1[v][2 * sIdx1 + 1] = 1;
                    } else if (gt == 1) {
                        hap1[v][2 * sIdx1 + phase] = 1;
                    }
                }

                Integer sIdx2 = sampleIndex2.get(sid);
                if (sIdx2 != null) {
                    if (gt == 2) {
                        hap2[v][2 * sIdx2] = 1;
                        hap2[v][2 * sIdx2 + 1] = 1;
                    } else if (gt == 1) {
                        hap2[v][2 * sIdx2 + phase] = 1;
                    }
                }
            }
        }

        HaplotypeMatrix m1 = new HaplotypeMatrix(variantIds, positions, nodes,
                afs1, acs1, ans1, hap1, nSamples1, sampleIndex1);
        HaplotypeMatrix m2 = new HaplotypeMatrix(variantIds, positions, nodes,
                afs2, acs2, ans2, hap2, nSamples2, sampleIndex2);
        return new HaplotypeMatrix[]{m1, m2};
    }
}
