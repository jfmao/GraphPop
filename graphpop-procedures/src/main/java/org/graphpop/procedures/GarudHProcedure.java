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
 * Computes Garud et al. (2015) haplotype homozygosity statistics in sliding windows.
 *
 * <p>H1/H12/H2_H1 detect hard and soft selective sweeps from haplotype
 * frequency distributions. Also computes haplotype diversity Hd.</p>
 *
 * <p>Architecture: load → parallel compute → return.
 * <ol>
 *   <li>Load dense haplotype matrix via {@link HaplotypeMatrix} (single-threaded)</li>
 *   <li>Compute H statistics per window in parallel (ForkJoinPool)</li>
 *   <li>Return results (no DB writes)</li>
 * </ol>
 * </p>
 *
 * <p>Usage:
 * <pre>
 * CALL graphpop.garud_h('chr22', 'EUR', 100000, 50000, {min_af: 0.05})
 * YIELD chr, start, end, h1, h12, h2_h1, hap_diversity, n_haplotypes, n_variants
 * </pre>
 * </p>
 */
public class GarudHProcedure {

    @Context
    public Transaction tx;

    public static class GarudHResult {
        public String chr;
        public long start;
        public long end;
        public double h1;
        public double h12;
        public double h2_h1;
        public double hap_diversity;
        public long n_haplotypes;
        public long n_variants;

        public GarudHResult() {}
    }

    @Procedure(name = "graphpop.garud_h", mode = Mode.READ)
    @SuppressWarnings("unchecked")
    public Stream<GarudHResult> garudH(
            @Name("chr") String chr,
            @Name("pop") String pop,
            @Name(value = "window", defaultValue = "100000") long windowSize,
            @Name(value = "step", defaultValue = "50000") long stepSize,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        // Parse options
        Long start = options != null && options.containsKey("start")
                ? ((Number) options.get("start")).longValue() : null;
        Long end = options != null && options.containsKey("end")
                ? ((Number) options.get("end")).longValue() : null;
        List<String> samples = options != null ? (List<String>) options.get("samples") : null;
        String variantType = options != null && options.containsKey("variant_type")
                ? (String) options.get("variant_type") : "SNP";
        if ("ALL".equals(variantType)) variantType = null;  // explicit override

        // Load haplotype matrix
        HaplotypeMatrix matrix;
        if (samples != null && !samples.isEmpty()) {
            Map<String, Integer> sampleIndex = GenotypeLoader.buildSampleIndex(tx, samples);
            matrix = HaplotypeMatrix.load(tx, chr, sampleIndex, start, end, variantType);
        } else {
            matrix = HaplotypeMatrix.load(tx, chr, pop, start, end, variantType);
        }

        if (matrix == null || matrix.nVariants == 0) {
            return Stream.empty();
        }

        // Build window boundaries
        long minPos = matrix.positions[0];
        long maxPos = matrix.positions[matrix.nVariants - 1];

        List<long[]> windowBounds = new ArrayList<>();
        for (long wStart = minPos; wStart <= maxPos; wStart += stepSize) {
            windowBounds.add(new long[]{wStart, wStart + windowSize - 1});
        }

        int nWindows = windowBounds.size();
        if (nWindows == 0) return Stream.empty();

        // Parallel window computation
        GarudHResult[] results = new GarudHResult[nWindows];
        final HaplotypeMatrix mat = matrix;

        IntStream.range(0, nWindows).parallel().forEach(w -> {
            long wStart = windowBounds.get(w)[0];
            long wEnd = windowBounds.get(w)[1];

            // Binary search for window boundaries in the position array
            int lo = Arrays.binarySearch(mat.positions, wStart);
            if (lo < 0) lo = -(lo + 1);
            int hi = Arrays.binarySearch(mat.positions, wEnd);
            if (hi < 0) hi = -(hi + 1) - 1;

            if (lo > hi || lo >= mat.nVariants) return;

            int nVarInWindow = hi - lo + 1;
            if (nVarInWindow < 1) return;

            GarudH.Result hr = GarudH.compute(mat.haplotypes, lo, hi, mat.nHaplotypes);
            if (hr == null) return;

            GarudHResult r = new GarudHResult();
            r.chr = chr;
            r.start = wStart;
            r.end = wEnd;
            r.h1 = hr.h1;
            r.h12 = hr.h12;
            r.h2_h1 = hr.h2_h1;
            r.hap_diversity = hr.hap_diversity;
            r.n_haplotypes = hr.n_haplotypes;
            r.n_variants = nVarInWindow;

            results[w] = r;
        });

        // Collect non-null results
        List<GarudHResult> output = new ArrayList<>();
        for (GarudHResult r : results) {
            if (r != null) output.add(r);
        }

        return output.stream();
    }
}
