package org.graphpop.procedures;

import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Transaction;
import org.neo4j.procedure.Context;
import org.neo4j.procedure.Description;
import org.neo4j.procedure.Mode;
import org.neo4j.procedure.Name;
import org.neo4j.procedure.Procedure;

import java.util.*;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Computes Cross-Population Extended Haplotype Homozygosity (XP-EHH).
 *
 * <p>Architecture: chunked load → parallel compute → genome-wide standardize → sequential write.
 * <ol>
 *   <li>Determine chromosome bounds</li>
 *   <li>Iterate in genomic windows (default 5 Mb core + 2 Mb margin each side):
 *     <ol>
 *       <li>Load paired haplotype matrices for window + margin via
 *           {@link HaplotypeMatrix#loadPair} (single-threaded)</li>
 *       <li>Compute iHH per population in parallel for focal variants in the
 *           core window only (ForkJoinPool)</li>
 *       <li>Collect unstandardized XP-EHH values; release matrices for GC</li>
 *     </ol>
 *   </li>
 *   <li>Standardize genome-wide (mean/std across all collected values)</li>
 *   <li>Write XP-EHH properties on Variant nodes (single-threaded)</li>
 * </ol>
 * </p>
 *
 * <p>Chunking reduces peak memory from O(nVariants_chr × nSamples × 2 pops) to
 * O(nVariants_chunk × nSamples × 2 pops). For chr1 with ~700k variants, this
 * drops peak usage from ~3–5 GB to ~100–200 MB per chunk.</p>
 *
 * <p>Usage:
 * <pre>
 * CALL graphpop.xpehh('chr22', 'EUR', 'AFR', {min_af: 0.05})
 * YIELD variantId, pos, af_pop1, af_pop2, xpehh_unstd, xpehh
 *
 * // Custom chunk size for memory-constrained environments:
 * CALL graphpop.xpehh('chr1', 'EUR', 'AFR', {chunk_size: 2000000, ehh_margin: 1000000})
 * </pre>
 * </p>
 */
public class XPEHHProcedure {

    // XP-EHH uses ALL haplotypes (not carrier-only like iHS), so there is
    // no statistical reason to require a minimum MAF. Reference tools (selscan,
    // scikit-allel) process every shared polymorphic variant. Default to 0.0.
    private static final double DEFAULT_MIN_AF = 0.0;

    /** Default core window size (bp) for chunked processing. */
    private static final long DEFAULT_CHUNK_SIZE = 5_000_000L;

    /** Default margin (bp) on each side of the core window for EHH walk context. */
    private static final long DEFAULT_EHH_MARGIN = 2_000_000L;

    /**
     * Maximum variants to process before switching to unstd_only mode automatically.
     * Beyond this threshold, result streaming exceeds Neo4j transaction memory.
     */
    private static final int AUTO_UNSTD_THRESHOLD = 500_000;

    @Context
    public Transaction tx;

    public static class XPEHHResult {
        public String variantId;
        public long pos;
        public double af_pop1;
        public double af_pop2;
        public double xpehh_unstd;
        public double xpehh;

        public XPEHHResult() {}
    }

    /** Lightweight per-variant result holder, collected across chunks. */
    private static class FocalResult {
        final String variantId;
        final long pos;
        final double afPop1, afPop2;
        final double xpehhUnstd;
        final Node node;

        FocalResult(String variantId, long pos, double afPop1, double afPop2,
                    double xpehhUnstd, Node node) {
            this.variantId = variantId;
            this.pos = pos;
            this.afPop1 = afPop1;
            this.afPop2 = afPop2;
            this.xpehhUnstd = xpehhUnstd;
            this.node = node;
        }
    }

    @Procedure(name = "graphpop.xpehh", mode = Mode.WRITE)
    @Description("Compute cross-population extended haplotype homozygosity (XP-EHH) between two populations.")
    @SuppressWarnings("unchecked")
    public Stream<XPEHHResult> xpehh(
            @Name("chr") String chr,
            @Name("pop1") String pop1,
            @Name("pop2") String pop2,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        // Default to SNPs only — XP-EHH EHH walk path must match reference tools
        // (selscan, rehh, scikit-allel) which all receive pre-filtered SNP input.
        // Indels in the walk path cause extra haplotype group splits that distort
        // the EHH decay curve. Users can override with variant_type: 'ALL'.
        if (options == null) {
            options = Map.of("variant_type", "SNP");
        } else if (!options.containsKey("variant_type")) {
            options = new HashMap<>(options);
            options.put("variant_type", "SNP");
        } else if ("ALL".equals(options.get("variant_type"))) {
            options = new HashMap<>(options);
            options.remove("variant_type");  // no type filter
        }
        VariantFilter filter = VariantFilter.fromOptions(options);
        double minAf = filter.minAf > 0 ? filter.minAf : DEFAULT_MIN_AF;
        Long regionStart = null;
        Long regionEnd = null;
        if (options != null) {
            if (options.containsKey("start"))
                regionStart = ((Number) options.get("start")).longValue();
            if (options.containsKey("end"))
                regionEnd = ((Number) options.get("end")).longValue();
        }

        // Parse sample lists and variant type
        List<String> samples1 = options != null ? (List<String>) options.get("samples1") : null;
        List<String> samples2 = options != null ? (List<String>) options.get("samples2") : null;
        String variantType = filter.variantType;
        EHHComputer ehhComputer = EHHComputer.fromOptions(options);

        boolean useCustomSamples = (samples1 != null && !samples1.isEmpty())
                || (samples2 != null && !samples2.isEmpty());

        // unstd_only mode: write only xpehh_unstd, return summary (no standardization,
        // no result streaming). Use for large chromosomes where per-region calls are
        // followed by external global standardization.
        boolean unstdOnly = options != null && Boolean.TRUE.equals(options.get("unstd_only"));

        // Chunking parameters
        long chunkSize = DEFAULT_CHUNK_SIZE;
        long ehhMargin = DEFAULT_EHH_MARGIN;
        if (options != null) {
            if (options.containsKey("chunk_size"))
                chunkSize = ((Number) options.get("chunk_size")).longValue();
            if (options.containsKey("ehh_margin"))
                ehhMargin = ((Number) options.get("ehh_margin")).longValue();
        }

        // Determine effective bounds for windowing
        long chrStart, chrEnd;
        if (regionStart != null && regionEnd != null) {
            chrStart = regionStart;
            chrEnd = regionEnd;
        } else {
            long[] bounds = queryChromosomeBounds(chr, variantType);
            if (bounds == null) return Stream.empty();
            chrStart = regionStart != null ? regionStart : bounds[0];
            chrEnd = regionEnd != null ? regionEnd : bounds[1];
        }

        // Pre-build sample indices for custom samples path (reuse across chunks)
        Map<String, Integer> customIdx1 = null, customIdx2 = null;
        if (useCustomSamples) {
            customIdx1 = (samples1 != null && !samples1.isEmpty())
                    ? GenotypeLoader.buildSampleIndex(tx, samples1)
                    : GenotypeLoader.buildSampleIndex(tx, pop1);
            customIdx2 = (samples2 != null && !samples2.isEmpty())
                    ? GenotypeLoader.buildSampleIndex(tx, samples2)
                    : GenotypeLoader.buildSampleIndex(tx, pop2);
        }

        // Collect results across all chunks
        List<FocalResult> allFocal = new ArrayList<>();

        for (long winStart = chrStart; winStart <= chrEnd; winStart += chunkSize) {
            long winCoreEnd = Math.min(winStart + chunkSize - 1, chrEnd);
            long loadStart = Math.max(winStart - ehhMargin, 0);
            long loadEnd = winCoreEnd + ehhMargin;

            // Load pair for this window (with margin for EHH walk context)
            HaplotypeMatrix[] matPair;
            if (useCustomSamples) {
                HaplotypeMatrix raw1 = HaplotypeMatrix.load(tx, chr, customIdx1,
                        loadStart, loadEnd, variantType);
                HaplotypeMatrix raw2 = HaplotypeMatrix.load(tx, chr, customIdx2,
                        loadStart, loadEnd, variantType);
                if (raw1 == null || raw2 == null) continue;
                matPair = HaplotypeMatrix.intersectPolymorphic(raw1, raw2);
                if (matPair[0].nVariants == 0) continue;
            } else {
                matPair = HaplotypeMatrix.loadPair(tx, chr, pop1, pop2,
                        loadStart, loadEnd, variantType, true);
                if (matPair == null) continue;
            }

            final HaplotypeMatrix m1 = matPair[0];
            final HaplotypeMatrix m2 = matPair[1];

            // Only process focal variants in core window [winStart, winCoreEnd]
            // Margin variants are present for EHH walk context but not focal
            int[] focalIndices = identifyFocalInWindow(m1, m2, minAf, filter,
                    winStart, winCoreEnd);
            if (focalIndices.length == 0) continue;

            // Parallel iHH computation per variant (no DB access)
            double[] xpehhUnstd = new double[focalIndices.length];
            boolean[] valid = new boolean[focalIndices.length];

            IntStream.range(0, focalIndices.length).parallel().forEach(i -> {
                int vidx = focalIndices[i];

                int[] allHaps1 = allHaplotypeIndices(m1.nHaplotypes);
                int[] allHaps2 = allHaplotypeIndices(m2.nHaplotypes);

                double iHH1 = ehhComputer.iHH(m1, vidx, allHaps1, true);
                double iHH2 = ehhComputer.iHH(m2, vidx, allHaps2, true);

                if (iHH1 > 0 && iHH2 > 0) {
                    xpehhUnstd[i] = Math.log(iHH1 / iHH2);
                    valid[i] = true;
                }
            });

            // Collect valid results (lightweight — matrices can be GC'd after this)
            for (int i = 0; i < focalIndices.length; i++) {
                if (!valid[i]) continue;
                int vidx = focalIndices[i];
                allFocal.add(new FocalResult(
                        m1.variantIds[vidx], m1.positions[vidx],
                        m1.afs[vidx], m2.afs[vidx],
                        xpehhUnstd[i], m1.nodes[vidx]));
            }
        }

        if (allFocal.isEmpty()) return Stream.empty();

        // Auto-detect: switch to unstd_only if too many variants for result streaming
        if (!unstdOnly && allFocal.size() > AUTO_UNSTD_THRESHOLD) {
            unstdOnly = true;
        }

        String propNameUnstd = "xpehh_unstd_" + pop1 + "_" + pop2;

        if (unstdOnly) {
            // Write only unstandardized values — caller handles global standardization.
            // Return a single summary row to avoid result streaming OOM.
            for (FocalResult fr : allFocal) {
                fr.node.setProperty(propNameUnstd, fr.xpehhUnstd);
            }

            XPEHHResult summary = new XPEHHResult();
            summary.variantId = "SUMMARY";
            summary.pos = allFocal.size();
            summary.af_pop1 = 0;
            summary.af_pop2 = 0;
            summary.xpehh_unstd = 0;
            summary.xpehh = 0;
            return Stream.of(summary);
        }

        // Genome-wide standardization (for small enough result sets)
        double sum = 0.0;
        for (FocalResult fr : allFocal) sum += fr.xpehhUnstd;
        double mean = sum / allFocal.size();

        double varSum = 0.0;
        for (FocalResult fr : allFocal) {
            double d = fr.xpehhUnstd - mean;
            varSum += d * d;
        }
        double std = allFocal.size() > 1 ? Math.sqrt(varSum / (allFocal.size() - 1)) : 1.0;

        // Sequential write to graph
        String propName = "xpehh_" + pop1 + "_" + pop2;

        List<XPEHHResult> results = new ArrayList<>();
        for (FocalResult fr : allFocal) {
            double xpehhStd = std > 0 ? (fr.xpehhUnstd - mean) / std : 0.0;

            fr.node.setProperty(propName, xpehhStd);
            fr.node.setProperty(propNameUnstd, fr.xpehhUnstd);

            XPEHHResult xr = new XPEHHResult();
            xr.variantId = fr.variantId;
            xr.pos = fr.pos;
            xr.af_pop1 = fr.afPop1;
            xr.af_pop2 = fr.afPop2;
            xr.xpehh_unstd = fr.xpehhUnstd;
            xr.xpehh = xpehhStd;
            results.add(xr);
        }

        results.sort(Comparator.comparingLong(r -> r.pos));
        return results.stream();
    }

    /**
     * Query min/max variant positions for a chromosome (optionally filtered by type).
     * Uses the composite index on (chr, pos) for efficiency.
     */
    private long[] queryChromosomeBounds(String chr, String variantType) {
        StringBuilder qb = new StringBuilder(
                "MATCH (v:Variant) WHERE v.chr = $chr");
        Map<String, Object> params = new HashMap<>();
        params.put("chr", chr);
        if (variantType != null) {
            qb.append(" AND v.variant_type = $variantType");
            params.put("variantType", variantType);
        }
        qb.append(" RETURN min(v.pos) AS minPos, max(v.pos) AS maxPos");

        var result = tx.execute(qb.toString(), params);
        try {
            if (!result.hasNext()) return null;
            var row = result.next();
            Object minObj = row.get("minPos"), maxObj = row.get("maxPos");
            if (minObj == null || maxObj == null) return null;
            return new long[]{((Number) minObj).longValue(), ((Number) maxObj).longValue()};
        } finally {
            result.close();
        }
    }

    /**
     * Identify focal variant indices within a core window [winStart, winEnd].
     * Margin variants (outside the core window) are excluded — they are present
     * in the matrix only to provide EHH walk context.
     */
    private static int[] identifyFocalInWindow(HaplotypeMatrix m1, HaplotypeMatrix m2,
                                                double minAf, VariantFilter filter,
                                                long winStart, long winEnd) {
        List<Integer> focal = new ArrayList<>();
        String variantType = filter.variantType;
        for (int i = 0; i < m1.nVariants; i++) {
            // Core window check: skip margin variants
            long pos = m1.positions[i];
            if (pos < winStart || pos > winEnd) continue;

            int an1 = m1.ans[i], an2 = m2.ans[i];
            int ac1 = m1.acs[i], ac2 = m2.acs[i];
            double af1 = m1.afs[i];

            if (an1 < 2 || an2 < 2) continue;
            if (ac1 == 0 && ac2 == 0) continue;
            if (ac1 == an1 && ac2 == an2) continue;

            // Require polymorphic in BOTH populations: scikit-allel/selscan
            // work on the intersection of per-population variant sets, so
            // only shared polymorphic sites are focal.
            boolean poly1 = ac1 > 0 && ac1 < an1;
            boolean poly2 = ac2 > 0 && ac2 < an2;
            if (!poly1 || !poly2) continue;

            // Apply MAF filter (if any) to pop1
            if (af1 < minAf || af1 > (1.0 - minAf)) continue;

            if (filter.isActive()) {
                if (af1 < filter.minAf || af1 > filter.maxAf) continue;
            }

            // Check variant_type from node property (defense-in-depth;
            // matrix is already filtered at load time when variantType is set)
            if (variantType != null && m1.nodes != null) {
                Object vt = m1.nodes[i].getProperty("variant_type", null);
                if (vt == null || !variantType.equals(vt)) continue;
            }

            focal.add(i);
        }
        return focal.stream().mapToInt(Integer::intValue).toArray();
    }

    private static int[] allHaplotypeIndices(int nHaplotypes) {
        int[] indices = new int[nHaplotypes];
        for (int i = 0; i < nHaplotypes; i++) indices[i] = i;
        return indices;
    }
}
