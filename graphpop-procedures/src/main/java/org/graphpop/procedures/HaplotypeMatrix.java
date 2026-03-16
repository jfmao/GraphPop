package org.graphpop.procedures;

import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Transaction;

import java.util.*;

/**
 * Dense in-memory haplotype matrix for a genomic region and population.
 *
 * <p>Loads all variant &times; haplotype data in a single pass via the graph API,
 * enabling parallel computation with zero DB access. The matrix layout is
 * bit-packed: {@code haplotypes[variantIndex][byteIndex]} where each bit
 * represents one haplotype (0=ref, 1=alt). Haplotype index
 * {@code h = 2 * sampleIndex + phase}; the bit is at
 * {@code haplotypes[v][h >> 3] & (1 << (h & 7))}.</p>
 *
 * <p>Memory usage: nVariants &times; ceil(2 &times; nSamples / 8) bytes.
 * For EUR on chr22: ~1M &times; 252 bytes &asymp; 252 MB (vs 2 GB unpacked).</p>
 */
final class HaplotypeMatrix {

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
     * Bit-packed haplotype matrix: {@code haplotypes[variantIdx][byteIdx]}.
     * Bit {@code h & 7} of byte {@code h >> 3} stores haplotype h.
     * h = 2 * sampleIdx + phase. Values: 0 (ref) or 1 (alt).
     */
    final byte[][] haplotypes;

    final int nVariants;
    final int nHaplotypes;
    final int nSamples;

    /** Sample index used during loading (sampleId &rarr; array position). */
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
     * Set a single bit in a packed byte array.
     */
    private static void setBit(byte[] packed, int idx) {
        packed[idx >> 3] |= (byte)(1 << (idx & 7));
    }

    /**
     * Pack an unpacked byte[] row (one byte per haplotype, values 0/1)
     * into a bit-packed byte[] row.
     */
    static byte[] packRow(byte[] unpacked) {
        int packedLen = (unpacked.length + 7) >> 3;
        byte[] packed = new byte[packedLen];
        for (int i = 0; i < unpacked.length; i++) {
            if (unpacked[i] == 1) {
                setBit(packed, i);
            }
        }
        return packed;
    }

    /**
     * Get a single haplotype value (0 or 1) from the bit-packed matrix.
     */
    int hap(int variantIdx, int hapIdx) {
        return (haplotypes[variantIdx][hapIdx >> 3] >> (hapIdx & 7)) & 1;
    }

    /**
     * Build a HaplotypeMatrix from raw unpacked arrays for unit testing (no DB required).
     * The unpacked byte[][] (one byte per haplotype) is auto-packed to bit format.
     */
    static HaplotypeMatrix forTest(String[] variantIds, long[] positions,
                                   double[] afs, int[] acs, int[] ans,
                                   byte[][] unpackedHaplotypes, int nSamples) {
        byte[][] packed = new byte[unpackedHaplotypes.length][];
        for (int v = 0; v < unpackedHaplotypes.length; v++) {
            packed[v] = packRow(unpackedHaplotypes[v]);
        }
        return new HaplotypeMatrix(variantIds, positions, null, afs, acs, ans,
                packed, nSamples, Map.of());
    }

    /**
     * Filter two independently-loaded matrices to their shared polymorphic
     * intersection. Both output matrices will contain exactly the same positions
     * in the same order, with only variants that are polymorphic (0 &lt; AC &lt; AN)
     * in both populations retained.
     *
     * <p>This matches scikit-allel/selscan XP-EHH behavior where input is the
     * shared intersection of per-population polymorphic variants.</p>
     */
    static HaplotypeMatrix[] intersectPolymorphic(HaplotypeMatrix a, HaplotypeMatrix b) {
        // Build position → index maps
        Map<Long, Integer> posToIdxA = new HashMap<>();
        for (int i = 0; i < a.nVariants; i++) {
            if (a.acs[i] > 0 && a.acs[i] < a.ans[i]) {
                posToIdxA.put(a.positions[i], i);
            }
        }

        // Find shared polymorphic positions (in sorted order from b)
        List<int[]> pairs = new ArrayList<>();  // [idxA, idxB]
        for (int j = 0; j < b.nVariants; j++) {
            if (b.acs[j] > 0 && b.acs[j] < b.ans[j]) {
                Integer idxA = posToIdxA.get(b.positions[j]);
                if (idxA != null) {
                    pairs.add(new int[]{idxA, j});
                }
            }
        }

        int n = pairs.size();
        String[] vidsA = new String[n], vidsB = new String[n];
        long[] posArr = new long[n];
        Node[] nodesA = a.nodes != null ? new Node[n] : null;
        Node[] nodesB = b.nodes != null ? new Node[n] : null;
        double[] afsA = new double[n], afsB = new double[n];
        int[] acsA = new int[n], acsB = new int[n];
        int[] ansA = new int[n], ansB = new int[n];
        byte[][] hapA = new byte[n][], hapB = new byte[n][];

        for (int k = 0; k < n; k++) {
            int iA = pairs.get(k)[0], iB = pairs.get(k)[1];
            vidsA[k] = a.variantIds[iA]; vidsB[k] = b.variantIds[iB];
            posArr[k] = a.positions[iA];
            if (nodesA != null) nodesA[k] = a.nodes[iA];
            if (nodesB != null) nodesB[k] = b.nodes[iB];
            afsA[k] = a.afs[iA]; afsB[k] = b.afs[iB];
            acsA[k] = a.acs[iA]; acsB[k] = b.acs[iB];
            ansA[k] = a.ans[iA]; ansB[k] = b.ans[iB];
            hapA[k] = a.haplotypes[iA]; hapB[k] = b.haplotypes[iB];
        }

        HaplotypeMatrix mA = new HaplotypeMatrix(vidsA, posArr, nodesA,
                afsA, acsA, ansA, hapA, a.nSamples, a.sampleIndex);
        HaplotypeMatrix mB = new HaplotypeMatrix(vidsB, posArr, nodesB,
                afsB, acsB, ansB, hapB, b.nSamples, b.sampleIndex);
        return new HaplotypeMatrix[]{mA, mB};
    }

    /**
     * Get the dosage vector for a variant (0/1/2 per sample).
     * dosage[s] = hap(v, 2*s) + hap(v, 2*s+1).
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
        for (int s = 0; s < nSamples; s++) {
            dos[s] = hap(variantIdx, 2 * s) + hap(variantIdx, 2 * s + 1);
        }
        return dos;
    }

    /**
     * Get a single haplotype value (legacy accessor).
     */
    int haplotype(int variantIdx, int sampleIdx, int phase) {
        return hap(variantIdx, 2 * sampleIdx + phase);
    }

    /**
     * Load a haplotype matrix for a chromosome region using a pre-built sample index.
     * This overload supports custom sample subsets — the caller builds the sample index
     * from an explicit list of sample IDs.
     *
     * @param tx          active transaction
     * @param chr         chromosome
     * @param sampleIndex pre-built sample index (sampleId &rarr; array position)
     * @param start       region start (use {@code null} for whole chromosome)
     * @param end         region end (use {@code null} for whole chromosome)
     * @return loaded matrix, or null if no variants/samples found
     */
    static HaplotypeMatrix load(Transaction tx, String chr,
                                Map<String, Integer> sampleIndex,
                                Long start, Long end) {
        return loadInternal(tx, chr, sampleIndex, start, end, null);
    }

    /**
     * Load a haplotype matrix filtered by variant type.
     *
     * @param variantType variant type filter (e.g. "SNP"), or null for all types
     */
    static HaplotypeMatrix load(Transaction tx, String chr,
                                Map<String, Integer> sampleIndex,
                                Long start, Long end, String variantType) {
        return loadInternal(tx, chr, sampleIndex, start, end, variantType);
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
        Map<String, Integer> sampleIndex = GenotypeLoader.buildSampleIndex(tx, pop);
        return loadInternal(tx, chr, sampleIndex, start, end, null);
    }

    /**
     * Load a haplotype matrix for a population, filtered by variant type.
     *
     * @param variantType variant type filter (e.g. "SNP"), or null for all types
     */
    static HaplotypeMatrix load(Transaction tx, String chr, String pop,
                                Long start, Long end, String variantType) {
        Map<String, Integer> sampleIndex = GenotypeLoader.buildSampleIndex(tx, pop);
        return loadInternal(tx, chr, sampleIndex, start, end, variantType);
    }

    private static HaplotypeMatrix loadInternal(Transaction tx, String chr,
                                                Map<String, Integer> sampleIndex,
                                                Long start, Long end,
                                                String variantType) {
        int nSamples = sampleIndex.size();
        if (nSamples == 0) return null;
        int nHaplotypes = nSamples * 2;
        int packedRowLen = (nHaplotypes + 7) >> 3;

        // 1. Build packed_index lookup: matrixPos → VCF column index
        int[] packedIndices = GenotypeLoader.buildPackedIndices(tx, sampleIndex);

        // 2. Query variants in region sorted by position, optionally filtered by type
        StringBuilder qb = new StringBuilder("MATCH (v:Variant) WHERE v.chr = $chr");
        Map<String, Object> params = new HashMap<>();
        params.put("chr", chr);
        if (start != null && end != null) {
            qb.append(" AND v.pos >= $start AND v.pos <= $end");
            params.put("start", start);
            params.put("end", end);
        }
        if (variantType != null) {
            qb.append(" AND v.variant_type = $variantType");
            params.put("variantType", variantType);
        }
        qb.append(" RETURN v ORDER BY v.pos");
        String query = qb.toString();

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

        // 5. Load haplotype data from packed arrays on Variant nodes
        //    Compute ac/an/af from packed genotypes (ploidy-aware)
        //    Bit-packed: one bit per haplotype
        byte[][] haplotypes = new byte[nVariants][packedRowLen];
        int[] acs = new int[nVariants];
        int[] ans = new int[nVariants];
        double[] afs = new double[nVariants];

        // Read ploidy_packed once from first variant (constant per chromosome)
        byte[] ploidyPacked = null;
        {
            Object obj = nodes[0].getProperty("ploidy_packed", null);
            if (obj instanceof byte[] pb && pb.length > 0) ploidyPacked = pb;
        }

        for (int v = 0; v < nVariants; v++) {
            Node varNode = nodes[v];
            byte[] hapRow = haplotypes[v];

            byte[] gtPacked = (byte[]) varNode.getProperty("gt_packed");
            byte[] phasePacked = (byte[]) varNode.getProperty("phase_packed");

            int ac = 0, an = 0;
            for (int s = 0; s < nSamples; s++) {
                int pi = packedIndices[s];
                if (pi < 0) continue;
                int gt = PackedGenotypeReader.genotype(gtPacked, pi);

                if (gt == PackedGenotypeReader.GT_HOM_ALT) {
                    setBit(hapRow, 2 * s);
                    setBit(hapRow, 2 * s + 1);
                } else if (gt == PackedGenotypeReader.GT_HET) {
                    int phase = PackedGenotypeReader.phase(phasePacked, pi);
                    setBit(hapRow, 2 * s + phase);
                }
                // GT_HOM_REF → both 0 (already zero-initialized)

                if (gt == PackedGenotypeReader.GT_MISSING) continue;
                boolean isHaploid = PackedGenotypeReader.ploidy(ploidyPacked, pi) == 1;
                int ploidy = isHaploid ? 1 : 2;
                an += ploidy;
                if (gt == PackedGenotypeReader.GT_HOM_ALT) ac += ploidy;
                else if (gt == PackedGenotypeReader.GT_HET) ac += 1;
            }

            acs[v] = ac;
            ans[v] = an;
            afs[v] = an > 0 ? (double) ac / an : 0.0;
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
        return loadPair(tx, chr, pop1, pop2, start, end, null);
    }

    /**
     * Load two population matrices, optionally filtered by variant type.
     *
     * @param variantType variant type filter (e.g. "SNP"), or null for all types
     */
    static HaplotypeMatrix[] loadPair(Transaction tx, String chr,
                                      String pop1, String pop2,
                                      Long start, Long end,
                                      String variantType) {
        return loadPair(tx, chr, pop1, pop2, start, end, variantType, false);
    }

    /**
     * Load paired haplotype matrices from a single variant pass.
     *
     * @param sharedPolyOnly if true, only include variants polymorphic in BOTH
     *                       populations (matching scikit-allel/selscan XP-EHH behavior
     *                       where input is the shared intersection)
     */
    static HaplotypeMatrix[] loadPair(Transaction tx, String chr,
                                      String pop1, String pop2,
                                      Long start, Long end,
                                      String variantType,
                                      boolean sharedPolyOnly) {
        Map<String, Integer> sampleIndex1 = GenotypeLoader.buildSampleIndex(tx, pop1);
        Map<String, Integer> sampleIndex2 = GenotypeLoader.buildSampleIndex(tx, pop2);
        int nSamples1 = sampleIndex1.size();
        int nSamples2 = sampleIndex2.size();
        if (nSamples1 == 0 || nSamples2 == 0) return null;

        int nHaplotypes1 = nSamples1 * 2;
        int nHaplotypes2 = nSamples2 * 2;
        int packedRowLen1 = (nHaplotypes1 + 7) >> 3;
        int packedRowLen2 = (nHaplotypes2 + 7) >> 3;

        // Build packed_index lookups for both populations
        int[] packedIndices1 = GenotypeLoader.buildPackedIndices(tx, sampleIndex1);
        int[] packedIndices2 = GenotypeLoader.buildPackedIndices(tx, sampleIndex2);

        // Query variants once, optionally filtered by type
        StringBuilder qb = new StringBuilder("MATCH (v:Variant) WHERE v.chr = $chr");
        Map<String, Object> params = new HashMap<>();
        params.put("chr", chr);
        if (start != null && end != null) {
            qb.append(" AND v.pos >= $start AND v.pos <= $end");
            params.put("start", start);
            params.put("end", end);
        }
        if (variantType != null) {
            qb.append(" AND v.variant_type = $variantType");
            params.put("variantType", variantType);
        }
        qb.append(" RETURN v ORDER BY v.pos");
        String query = qb.toString();

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

                // When sharedPolyOnly, skip variants not polymorphic in both pops.
                if (sharedPolyOnly) {
                    int[] acArr = ArrayUtil.toIntArray(node.getProperty("ac"));
                    int[] anArr = ArrayUtil.toIntArray(node.getProperty("an"));
                    int ac1 = acArr[ctx1.index], an1 = anArr[ctx1.index];
                    int ac2 = acArr[ctx2.index], an2 = anArr[ctx2.index];
                    boolean poly1 = ac1 > 0 && ac1 < an1;
                    boolean poly2 = ac2 > 0 && ac2 < an2;
                    if (!poly1 || !poly2) continue;
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

        // Load haplotypes for both populations from packed arrays (bit-packed)
        byte[][] hap1 = new byte[nVariants][packedRowLen1];
        byte[][] hap2 = new byte[nVariants][packedRowLen2];

        for (int v = 0; v < nVariants; v++) {
            Node varNode = nodes[v];
            byte[] gtPacked = (byte[]) varNode.getProperty("gt_packed");
            byte[] phasePacked = (byte[]) varNode.getProperty("phase_packed");

            // Population 1
            for (int s = 0; s < nSamples1; s++) {
                int pi = packedIndices1[s];
                if (pi < 0) continue;
                int gt = PackedGenotypeReader.genotype(gtPacked, pi);
                if (gt == PackedGenotypeReader.GT_HOM_ALT) {
                    setBit(hap1[v], 2 * s);
                    setBit(hap1[v], 2 * s + 1);
                } else if (gt == PackedGenotypeReader.GT_HET) {
                    int phase = PackedGenotypeReader.phase(phasePacked, pi);
                    setBit(hap1[v], 2 * s + phase);
                }
            }

            // Population 2
            for (int s = 0; s < nSamples2; s++) {
                int pi = packedIndices2[s];
                if (pi < 0) continue;
                int gt = PackedGenotypeReader.genotype(gtPacked, pi);
                if (gt == PackedGenotypeReader.GT_HOM_ALT) {
                    setBit(hap2[v], 2 * s);
                    setBit(hap2[v], 2 * s + 1);
                } else if (gt == PackedGenotypeReader.GT_HET) {
                    int phase = PackedGenotypeReader.phase(phasePacked, pi);
                    setBit(hap2[v], 2 * s + phase);
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
