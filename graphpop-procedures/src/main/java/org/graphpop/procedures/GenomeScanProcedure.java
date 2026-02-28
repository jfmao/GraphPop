package org.graphpop.procedures;

import org.graphpop.VectorOps;
import org.neo4j.graphdb.Label;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.RelationshipType;
import org.neo4j.graphdb.Transaction;
import org.neo4j.procedure.Context;
import org.neo4j.procedure.Mode;
import org.neo4j.procedure.Name;
import org.neo4j.procedure.Procedure;

import java.time.Instant;
import java.time.ZoneOffset;
import java.time.format.DateTimeFormatter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Sliding-window genome scan that materializes GenomicWindow nodes.
 *
 * <p>Architecture: load → parallel compute → sequential write.
 * <ol>
 *   <li>Load all variants into memory (single-threaded, one pass)</li>
 *   <li>Pre-compute window boundaries, then compute per-window statistics
 *       in parallel (ForkJoinPool, binary search for variant ranges)</li>
 *   <li>Write GenomicWindow nodes and ON_CHROMOSOME edges (single-threaded)</li>
 * </ol>
 * </p>
 *
 * <p>Each run is tagged with a run_id for versioning. Comparative mode
 * (providing pop2) adds Fst, Dxy, Da to each window.</p>
 *
 * <p>Usage:
 * <pre>
 * CALL graphpop.genome_scan('chr22', 'AFR', 100000, 50000, {})
 * YIELD window_id, chr, start, end, pi, theta_w, tajima_d, n_variants, run_id
 * </pre>
 *
 * Comparative mode:
 * <pre>
 * CALL graphpop.genome_scan('chr22', 'AFR', 100000, 50000, {pop2: 'EUR'})
 * YIELD window_id, fst, dxy, da
 * </pre>
 * </p>
 */
public class GenomeScanProcedure {

    private static final Label GENOMIC_WINDOW_LABEL = Label.label("GenomicWindow");
    private static final Label CHROMOSOME_LABEL = Label.label("Chromosome");
    private static final RelationshipType ON_CHROMOSOME = RelationshipType.withName("ON_CHROMOSOME");

    private static final DateTimeFormatter RUN_TS =
            DateTimeFormatter.ofPattern("yyyyMMdd'T'HHmm").withZone(ZoneOffset.UTC);

    @Context
    public Transaction tx;

    public static class WindowResult {
        public String window_id;
        public String chr;
        public long start;
        public long end;
        public String population;
        public String run_id;
        public long n_variants;
        public long n_segregating;
        public double pi;
        public double theta_w;
        public double tajima_d;
        public double het_exp;
        public double het_obs;
        public double fis;
        // Comparative (only set if pop2 provided)
        public double fst;
        public double dxy;
        public double da;

        public WindowResult() {}
    }

    @Procedure(name = "graphpop.genome_scan", mode = Mode.WRITE)
    public Stream<WindowResult> genomeScan(
            @Name("chr") String chr,
            @Name("pop") String pop,
            @Name(value = "window", defaultValue = "100000") long windowSize,
            @Name(value = "step", defaultValue = "50000") long stepSize,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        String pop2 = options != null ? (String) options.get("pop2") : null;
        String runIdOverride = options != null ? (String) options.get("run_id") : null;

        // Generate run_id
        String runId;
        if (runIdOverride != null) {
            runId = runIdOverride;
        } else {
            String ts = RUN_TS.format(Instant.now());
            String suffix = pop2 != null ? pop + "_" + pop2 : pop;
            runId = "scan_" + suffix + "_w" + (windowSize / 1000) + "k_s"
                    + (stepSize / 1000) + "k_" + ts;
        }

        // Fetch all variants on this chromosome, sorted by position.
        // The ORDER BY is critical for the sliding window logic.
        var result = tx.execute(
                "MATCH (v:Variant) WHERE v.chr = $chr RETURN v ORDER BY v.pos",
                Map.of("chr", chr)
        );

        // Load variants into memory (position + property indices).
        // For chr22 this is ~1M variants — manageable.
        List<VariantData> variants = new ArrayList<>();
        PopulationContext ctx = null;
        PopulationContext ctx2 = null;

        try {
            while (result.hasNext()) {
                Map<String, Object> row = result.next();
                Node variant = (Node) row.get("v");

                if (ctx == null) {
                    ctx = PopulationContext.resolve(tx, variant, pop);
                    if (pop2 != null) {
                        ctx2 = PopulationContext.resolve(tx, variant, pop2);
                    }
                }

                long pos = ((Number) variant.getProperty("pos")).longValue();
                int[] acArr = ArrayUtil.toIntArray(variant.getProperty("ac"));
                int[] anArr = ArrayUtil.toIntArray(variant.getProperty("an"));
                double[] afArr = ArrayUtil.toDoubleArray(variant.getProperty("af"));
                int[] hetArr = ArrayUtil.toIntArray(variant.getProperty("het_count"));

                variants.add(new VariantData(pos, acArr, anArr, afArr, hetArr));
            }
        } finally {
            result.close();
        }

        if (variants.isEmpty() || ctx == null) {
            return Stream.empty();
        }

        // Determine chromosome span
        long minPos = variants.get(0).pos;
        long maxPos = variants.get(variants.size() - 1).pos;

        // Build sorted position array for binary search
        long[] positions = new long[variants.size()];
        for (int i = 0; i < variants.size(); i++) {
            positions[i] = variants.get(i).pos;
        }

        // Phase 1: Pre-compute window boundaries
        List<long[]> windowBounds = new ArrayList<>();
        for (long wStart = minPos; wStart <= maxPos; wStart += stepSize) {
            windowBounds.add(new long[]{wStart, wStart + windowSize - 1});
        }
        int nWindows = windowBounds.size();
        if (nWindows == 0) return Stream.empty();

        // Phase 2: Parallel computation of per-window statistics (no DB access)
        final PopulationContext ctxF = ctx;
        final PopulationContext ctx2F = ctx2;
        WindowResult[] windowResults = new WindowResult[nWindows];

        IntStream.range(0, nWindows).parallel().forEach(w -> {
            long wStart = windowBounds.get(w)[0];
            long wEnd = windowBounds.get(w)[1];

            // Binary search for first variant >= wStart
            int lo = Arrays.binarySearch(positions, wStart);
            if (lo < 0) lo = -(lo + 1); // insertion point

            // Accumulate stats for variants in [wStart, wEnd]
            double piSum = 0.0;
            double heSum = 0.0;
            double hoSum = 0.0;
            long nVariants = 0;
            long nSegregating = 0;

            // Comparative accumulators
            double fstNum = 0.0, fstDen = 0.0, dxySum = 0.0;
            double piW1Sum = 0.0, piW2Sum = 0.0;

            for (int j = lo; j < variants.size(); j++) {
                VariantData v = variants.get(j);
                if (v.pos > wEnd) break;

                int ac = v.ac[ctxF.index];
                int an = v.an[ctxF.index];
                double af = v.af[ctxF.index];
                int het = v.het[ctxF.index];

                if (an < 2) continue;

                nVariants++;
                piSum += VectorOps.piPerSite(af, an);
                heSum += 2.0 * af * (1.0 - af);
                hoSum += (double) het / (an / 2);

                if (ac > 0 && ac < an) {
                    nSegregating++;
                }

                // Comparative stats
                if (ctx2F != null) {
                    int an2 = v.an[ctx2F.index];
                    double af2 = v.af[ctx2F.index];
                    if (an2 >= 2) {
                        double[] fstC = VectorOps.hudsonFstComponents(af, an, af2, an2);
                        fstNum += fstC[0];
                        fstDen += fstC[1];
                        dxySum += VectorOps.dxyPerSite(af, af2);
                        piW1Sum += VectorOps.piPerSite(af, an);
                        piW2Sum += VectorOps.piPerSite(af2, an2);
                    }
                }
            }

            if (nVariants == 0) return;

            double L = nVariants;

            WindowResult wr = new WindowResult();
            wr.chr = chr;
            wr.start = wStart;
            wr.end = wEnd;
            wr.population = pop;
            wr.run_id = runId;
            wr.n_variants = nVariants;
            wr.n_segregating = nSegregating;
            wr.pi = piSum / L;
            wr.theta_w = ctxF.a_n > 0 ? (nSegregating / ctxF.a_n) / L : 0.0;
            wr.tajima_d = TajimaD.compute(piSum, nSegregating, ctxF);
            wr.het_exp = heSum / L;
            wr.het_obs = hoSum / L;
            wr.fis = wr.het_exp > 0 ? 1.0 - wr.het_obs / wr.het_exp : 0.0;

            if (ctx2F != null && fstDen > 0) {
                wr.fst = fstNum / fstDen;
                wr.dxy = dxySum / L;
                wr.da = wr.dxy - (piW1Sum / L + piW2Sum / L) / 2.0;
            }

            wr.window_id = chr + ":" + wStart + "-" + wEnd + ":" + pop + ":" + runId;
            windowResults[w] = wr;
        });

        // Phase 3: Sequential write of GenomicWindow nodes (DB access)
        Node chrNode = tx.findNode(CHROMOSOME_LABEL, "chromosomeId", chr);
        List<WindowResult> results = new ArrayList<>();

        for (int w = 0; w < nWindows; w++) {
            WindowResult wr = windowResults[w];
            if (wr == null) continue; // empty window

            Node gwNode = tx.createNode(GENOMIC_WINDOW_LABEL);
            gwNode.setProperty("windowId", wr.window_id);
            gwNode.setProperty("chr", chr);
            gwNode.setProperty("start", wr.start);
            gwNode.setProperty("end", wr.end);
            gwNode.setProperty("population", pop);
            gwNode.setProperty("run_id", runId);
            gwNode.setProperty("n_variants", wr.n_variants);
            gwNode.setProperty("n_segregating", wr.n_segregating);
            gwNode.setProperty("pi", wr.pi);
            gwNode.setProperty("theta_w", wr.theta_w);
            gwNode.setProperty("tajima_d", wr.tajima_d);
            gwNode.setProperty("het_exp", wr.het_exp);
            gwNode.setProperty("het_obs", wr.het_obs);
            gwNode.setProperty("fis", wr.fis);

            if (pop2 != null) {
                gwNode.setProperty("pop2", pop2);
                gwNode.setProperty("fst", wr.fst);
                gwNode.setProperty("dxy", wr.dxy);
                gwNode.setProperty("da", wr.da);
            }

            if (chrNode != null) {
                gwNode.createRelationshipTo(chrNode, ON_CHROMOSOME);
            }

            results.add(wr);
        }

        return results.stream();
    }

    /**
     * Lightweight holder for variant data needed by the sliding window.
     * Avoids re-reading Node properties on overlapping windows.
     */
    private static class VariantData {
        final long pos;
        final int[] ac;
        final int[] an;
        final double[] af;
        final int[] het;

        VariantData(long pos, int[] ac, int[] an, double[] af, int[] het) {
            this.pos = pos;
            this.ac = ac;
            this.an = an;
            this.af = af;
            this.het = het;
        }
    }
}
