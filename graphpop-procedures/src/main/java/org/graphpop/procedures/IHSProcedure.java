package org.graphpop.procedures;

import org.neo4j.graphdb.Transaction;
import org.neo4j.procedure.Context;
import org.neo4j.procedure.Mode;
import org.neo4j.procedure.Name;
import org.neo4j.procedure.Procedure;

import java.util.*;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Computes the Integrated Haplotype Score (iHS) for detecting recent positive selection.
 *
 * <p>Architecture: load → parallel compute → standardize → sequential write.
 * <ol>
 *   <li>Load dense haplotype matrix via {@link HaplotypeMatrix} (single-threaded)</li>
 *   <li>Compute iHH_derived and iHH_ancestral in parallel per variant (ForkJoinPool)</li>
 *   <li>Standardize within allele-frequency bins</li>
 *   <li>Write iHS properties on Variant nodes (single-threaded)</li>
 * </ol>
 * </p>
 *
 * <p>Usage:
 * <pre>
 * CALL graphpop.ihs('chr22', 'EUR', {min_af: 0.05, start: 16000000, end: 17000000})
 * YIELD variantId, pos, af, ihs_unstd, ihs
 * </pre>
 * </p>
 */
public class IHSProcedure {

    private static final double DEFAULT_MIN_AF = 0.05;
    private static final int AF_BINS = 20;

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

    @Procedure(name = "graphpop.ihs", mode = Mode.WRITE)
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

        // Phase 1: Load dense haplotype matrix (single-threaded)
        // When variantType is set (default "SNP"), only matching variants are loaded
        // so the EHH walk path only contains SNPs.
        HaplotypeMatrix matrix;
        if (useSubset) {
            Map<String, Integer> sampleIndex = GenotypeLoader.buildSampleIndex(tx, sampleList);
            matrix = HaplotypeMatrix.load(tx, chr, sampleIndex, regionStart, regionEnd, variantType);
        } else {
            matrix = HaplotypeMatrix.load(tx, chr, pop, regionStart, regionEnd, variantType);
        }
        if (matrix == null || matrix.nVariants == 0) return Stream.empty();

        // Identify focal variants: common, polymorphic
        int[] focalIndices = identifyFocal(matrix, minAf, filter);
        if (focalIndices.length == 0) return Stream.empty();

        // Phase 2: Parallel iHH computation (no DB access)
        double[] ihsUnstd = new double[focalIndices.length];
        boolean[] valid = new boolean[focalIndices.length];

        IntStream.range(0, focalIndices.length).parallel().forEach(i -> {
            int vidx = focalIndices[i];

            byte[] focalHaps = matrix.haplotypes[vidx];
            int nHaps = matrix.nHaplotypes;

            int nDerived = 0;
            for (int h = 0; h < nHaps; h++) {
                if (focalHaps[h] == 1) nDerived++;
            }
            int nAncestral = nHaps - nDerived;

            if (nDerived < 2 || nAncestral < 2) return;

            int[] derivedIdx = new int[nDerived];
            int[] ancestralIdx = new int[nAncestral];
            int di = 0, ai = 0;
            for (int h = 0; h < nHaps; h++) {
                if (focalHaps[h] == 1) {
                    derivedIdx[di++] = h;
                } else {
                    ancestralIdx[ai++] = h;
                }
            }

            double ihhDerived = EHHComputer.computeIHH(matrix, vidx, derivedIdx);
            double ihhAncestral = EHHComputer.computeIHH(matrix, vidx, ancestralIdx);

            if (ihhDerived > 0 && ihhAncestral > 0) {
                ihsUnstd[i] = Math.log(ihhAncestral / ihhDerived);
                valid[i] = true;
            }
        });

        // Phase 3: Standardize within allele-frequency bins
        Map<Integer, List<Integer>> bins = new HashMap<>();
        for (int i = 0; i < focalIndices.length; i++) {
            if (!valid[i]) continue;
            int bin = Math.min((int) (matrix.afs[focalIndices[i]] * AF_BINS), AF_BINS - 1);
            bins.computeIfAbsent(bin, k -> new ArrayList<>()).add(i);
        }

        double[] ihsStd = new double[focalIndices.length];
        for (List<Integer> binList : bins.values()) {
            double sum = 0.0;
            for (int i : binList) sum += ihsUnstd[i];
            double mean = sum / binList.size();

            double varSum = 0.0;
            for (int i : binList) {
                double d = ihsUnstd[i] - mean;
                varSum += d * d;
            }
            double std = binList.size() > 1 ? Math.sqrt(varSum / (binList.size() - 1)) : 1.0;

            for (int i : binList) {
                ihsStd[i] = std > 0 ? (ihsUnstd[i] - mean) / std : 0.0;
            }
        }

        // Phase 4: Sequential write to graph
        List<IHSResult> results = new ArrayList<>();
        for (int i = 0; i < focalIndices.length; i++) {
            if (!valid[i]) continue;
            int vidx = focalIndices[i];

            matrix.nodes[vidx].setProperty("ihs_" + pop, ihsStd[i]);
            matrix.nodes[vidx].setProperty("ihs_unstd_" + pop, ihsUnstd[i]);

            IHSResult ir = new IHSResult();
            ir.variantId = matrix.variantIds[vidx];
            ir.pos = matrix.positions[vidx];
            ir.af = matrix.afs[vidx];
            ir.ihs_unstd = ihsUnstd[i];
            ir.ihs = ihsStd[i];
            results.add(ir);
        }

        results.sort(Comparator.comparingLong(r -> r.pos));
        return results.stream();
    }

    private static int[] identifyFocal(HaplotypeMatrix matrix, double minAf, VariantFilter filter) {
        List<Integer> focal = new ArrayList<>();
        String variantType = filter.variantType;
        for (int i = 0; i < matrix.nVariants; i++) {
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
