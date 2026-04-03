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
 * Computes the Integrated Haplotype Score (iHS) for detecting recent positive selection.
 *
 * <p>Architecture: chunked load &rarr; parallel compute &rarr; AF-bin standardize &rarr; sequential write.
 * <ol>
 *   <li>Determine chromosome bounds</li>
 *   <li>Iterate in genomic windows (default 5 Mb core + 2 Mb margin each side):
 *     <ol>
 *       <li>Load dense haplotype matrix for window + margin via
 *           {@link HaplotypeMatrix} (single-threaded)</li>
 *       <li>Compute iHH_derived and iHH_ancestral in parallel for focal
 *           variants in the core window only (ForkJoinPool)</li>
 *       <li>Collect unstandardized iHS values; release matrix for GC</li>
 *     </ol>
 *   </li>
 *   <li>Standardize within allele-frequency bins across all chunks</li>
 *   <li>Write iHS properties on Variant nodes (single-threaded)</li>
 * </ol>
 * </p>
 *
 * <p>Chunking reduces peak memory from O(nVariants_chr &times; nSamples) to
 * O(nVariants_chunk &times; nSamples). For chr1 with ~700k SNPs, this
 * drops peak usage from ~1.5 GB to ~100&ndash;200 MB per chunk.</p>
 *
 * <p>Usage:
 * <pre>
 * CALL graphpop.ihs('chr22', 'EUR', {min_af: 0.05, start: 16000000, end: 17000000})
 * YIELD variantId, pos, af, ihs_unstd, ihs
 *
 * // Custom chunk size for memory-constrained environments:
 * CALL graphpop.ihs('chr1', 'EUR', {chunk_size: 2000000, ehh_margin: 1000000})
 *
 * // Unstd-only mode for external standardization:
 * CALL graphpop.ihs('chr1', 'EUR', {start: 1, end: 50000000, unstd_only: true})
 * </pre>
 * </p>
 */
public class IHSProcedure {

    private static final double DEFAULT_MIN_AF = 0.05;
    private static final int DEFAULT_AF_BINS = 20;

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

    public static class IHSResult {
        public String variantId;
        public long pos;
        public double af;
        public double ihs_unstd;
        public double ihs;

        public IHSResult() {}
    }

    /** Lightweight per-variant result holder, collected across chunks. */
    private static class IHSFocalResult {
        final String variantId;
        final long pos;
        final double af;
        final double ihsUnstd;
        final Node node;

        IHSFocalResult(String variantId, long pos, double af,
                       double ihsUnstd, Node node) {
            this.variantId = variantId;
            this.pos = pos;
            this.af = af;
            this.ihsUnstd = ihsUnstd;
            this.node = node;
        }
    }

    @Procedure(name = "graphpop.ihs", mode = Mode.WRITE)
    @Description("Compute the integrated haplotype score (iHS) for detecting ongoing positive selection.")
    @SuppressWarnings("unchecked")
    public Stream<IHSResult> ihs(
            @Name("chr") String chr,
            @Name("pop") String pop,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        // Default to SNPs only — iHS EHH walk path must match reference tools
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

        List<String> sampleList = options != null ? (List<String>) options.get("samples") : null;
        boolean useSubset = sampleList != null && !sampleList.isEmpty();
        String variantType = filter.variantType;
        int afBins = options != null && options.containsKey("n_af_bins")
                ? ((Number) options.get("n_af_bins")).intValue() : DEFAULT_AF_BINS;
        EHHComputer ehhComputer = EHHComputer.fromOptions(options);

        // unstd_only mode: write only ihs_unstd, return summary (no standardization,
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

        // Pre-build sample index for custom samples path (reuse across chunks)
        Map<String, Integer> customIdx = null;
        if (useSubset) {
            customIdx = GenotypeLoader.buildSampleIndex(tx, sampleList);
        }

        // Collect results across all chunks
        List<IHSFocalResult> allFocal = new ArrayList<>();

        for (long winStart = chrStart; winStart <= chrEnd; winStart += chunkSize) {
            long winCoreEnd = Math.min(winStart + chunkSize - 1, chrEnd);
            long loadStart = Math.max(winStart - ehhMargin, 0);
            long loadEnd = winCoreEnd + ehhMargin;

            // Load haplotype matrix for this window (with margin for EHH walk context)
            HaplotypeMatrix matrix;
            if (useSubset) {
                matrix = HaplotypeMatrix.load(tx, chr, customIdx,
                        loadStart, loadEnd, variantType);
            } else {
                matrix = HaplotypeMatrix.load(tx, chr, pop,
                        loadStart, loadEnd, variantType);
            }
            if (matrix == null || matrix.nVariants == 0) continue;

            // Only process focal variants in core window [winStart, winCoreEnd]
            // Margin variants are present for EHH walk context but not focal
            int[] focalIndices = identifyFocalInWindow(matrix, minAf, filter,
                    winStart, winCoreEnd);
            if (focalIndices.length == 0) continue;

            // Parallel iHH computation per variant (no DB access)
            double[] ihsUnstd = new double[focalIndices.length];
            boolean[] valid = new boolean[focalIndices.length];

            IntStream.range(0, focalIndices.length).parallel().forEach(i -> {
                int vidx = focalIndices[i];

                byte[] focalHaps = matrix.haplotypes[vidx];
                int nHaps = matrix.nHaplotypes;

                int nDerived = 0;
                for (int h = 0; h < nHaps; h++) {
                    if (((focalHaps[h >> 3] >> (h & 7)) & 1) == 1) nDerived++;
                }
                int nAncestral = nHaps - nDerived;

                if (nDerived < 2 || nAncestral < 2) return;

                int[] derivedIdx = new int[nDerived];
                int[] ancestralIdx = new int[nAncestral];
                int di = 0, ai = 0;
                for (int h = 0; h < nHaps; h++) {
                    if (((focalHaps[h >> 3] >> (h & 7)) & 1) == 1) {
                        derivedIdx[di++] = h;
                    } else {
                        ancestralIdx[ai++] = h;
                    }
                }

                double ihhDerived = ehhComputer.iHH(matrix, vidx, derivedIdx, true);
                double ihhAncestral = ehhComputer.iHH(matrix, vidx, ancestralIdx, true);

                if (ihhDerived > 0 && ihhAncestral > 0) {
                    ihsUnstd[i] = Math.log(ihhAncestral / ihhDerived);
                    valid[i] = true;
                }
            });

            // Collect valid results (lightweight — matrix can be GC'd after this)
            for (int i = 0; i < focalIndices.length; i++) {
                if (!valid[i]) continue;
                int vidx = focalIndices[i];
                allFocal.add(new IHSFocalResult(
                        matrix.variantIds[vidx], matrix.positions[vidx],
                        matrix.afs[vidx], ihsUnstd[i], matrix.nodes[vidx]));
            }
        }

        if (allFocal.isEmpty()) return Stream.empty();

        // Auto-detect: switch to unstd_only if too many variants for result streaming
        if (!unstdOnly && allFocal.size() > AUTO_UNSTD_THRESHOLD) {
            unstdOnly = true;
        }

        String propNameUnstd = "ihs_unstd_" + pop;

        if (unstdOnly) {
            // Write only unstandardized values — caller handles global standardization.
            // Return a single summary row to avoid result streaming OOM.
            for (IHSFocalResult fr : allFocal) {
                fr.node.setProperty(propNameUnstd, fr.ihsUnstd);
            }

            IHSResult summary = new IHSResult();
            summary.variantId = "SUMMARY";
            summary.pos = allFocal.size();
            summary.af = 0;
            summary.ihs_unstd = 0;
            summary.ihs = 0;
            return Stream.of(summary);
        }

        // AF-bin standardization across all collected results
        Map<Integer, List<Integer>> bins = new HashMap<>();
        for (int i = 0; i < allFocal.size(); i++) {
            int bin = Math.min((int) (allFocal.get(i).af * afBins), afBins - 1);
            bins.computeIfAbsent(bin, k -> new ArrayList<>()).add(i);
        }

        double[] ihsStd = new double[allFocal.size()];
        for (List<Integer> binList : bins.values()) {
            double sum = 0.0;
            for (int idx : binList) sum += allFocal.get(idx).ihsUnstd;
            double mean = sum / binList.size();

            double varSum = 0.0;
            for (int idx : binList) {
                double d = allFocal.get(idx).ihsUnstd - mean;
                varSum += d * d;
            }
            double std = binList.size() > 1 ? Math.sqrt(varSum / (binList.size() - 1)) : 1.0;

            for (int idx : binList) {
                ihsStd[idx] = std > 0 ? (allFocal.get(idx).ihsUnstd - mean) / std : 0.0;
            }
        }

        // Sequential write to graph
        String propName = "ihs_" + pop;
        List<IHSResult> results = new ArrayList<>();
        for (int i = 0; i < allFocal.size(); i++) {
            IHSFocalResult fr = allFocal.get(i);

            fr.node.setProperty(propName, ihsStd[i]);
            fr.node.setProperty(propNameUnstd, fr.ihsUnstd);

            IHSResult ir = new IHSResult();
            ir.variantId = fr.variantId;
            ir.pos = fr.pos;
            ir.af = fr.af;
            ir.ihs_unstd = fr.ihsUnstd;
            ir.ihs = ihsStd[i];
            results.add(ir);
        }

        results.sort(Comparator.comparingLong(r -> r.pos));
        return results.stream();
    }

    /**
     * Query min/max variant positions for a chromosome (optionally filtered by type).
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
    private static int[] identifyFocalInWindow(HaplotypeMatrix matrix, double minAf,
                                                VariantFilter filter,
                                                long winStart, long winEnd) {
        List<Integer> focal = new ArrayList<>();
        String variantType = filter.variantType;
        for (int i = 0; i < matrix.nVariants; i++) {
            // Core window check: skip margin variants
            long pos = matrix.positions[i];
            if (pos < winStart || pos > winEnd) continue;

            int ac = matrix.acs[i];
            int an = matrix.ans[i];
            double af = matrix.afs[i];
            if (an < 2 || ac == 0 || ac == an) continue;
            if (af < minAf || af > (1.0 - minAf)) continue;

            if (filter.isActive()) {
                if (af < filter.minAf || af > filter.maxAf) continue;
            }

            // Check variant_type from node property (defense-in-depth;
            // matrix is already filtered at load time when variantType is set)
            if (variantType != null && matrix.nodes != null) {
                Object vt = matrix.nodes[i].getProperty("variant_type", null);
                if (vt == null || !variantType.equals(vt)) continue;
            }

            focal.add(i);
        }
        return focal.stream().mapToInt(Integer::intValue).toArray();
    }
}
