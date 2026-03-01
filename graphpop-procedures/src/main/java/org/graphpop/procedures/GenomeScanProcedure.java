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
 * <p>Architecture: load → pre-compute → parallel window → sequential write.
 * <ol>
 *   <li>Load all variants into memory (single-threaded, one pass)</li>
 *   <li>Pre-compute per-variant statistics into parallel arrays (M1.4 optimization)</li>
 *   <li>Parallel window summation — pure array ops, cache-friendly, JIT-vectorizable</li>
 *   <li>Write GenomicWindow nodes and ON_CHROMOSOME edges (single-threaded)</li>
 * </ol>
 * </p>
 *
 * <p>M1.4 optimization: each per-variant stat is computed exactly once into flat arrays,
 * eliminating ~50% duplicate work from overlapping windows.</p>
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
        public double fay_wu_h;
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
        VariantFilter filter = VariantFilter.fromOptions(options);

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
        var result = tx.execute(
                VariantQuery.buildChromosome(options),
                Map.of("chr", chr)
        );

        // Load variants into memory
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
                int[] homAltArr = ArrayUtil.toIntArray(variant.getProperty("hom_alt_count"));

                // Ancestral allele for Fay & Wu's H
                Object ancestralObj = variant.getProperty("ancestral_allele", null);
                String ancestralAllele = ancestralObj != null ? (String) ancestralObj : null;

                variants.add(new VariantData(pos, acArr, anArr, afArr, hetArr, homAltArr,
                        ancestralAllele, variant));
            }
        } finally {
            result.close();
        }

        if (variants.isEmpty() || ctx == null) {
            return Stream.empty();
        }

        int nVar = variants.size();

        // Build sorted position array for binary search
        long[] positions = new long[nVar];
        for (int i = 0; i < nVar; i++) {
            positions[i] = variants.get(i).pos;
        }

        // ---- M1.4: Pre-compute per-variant statistics into flat arrays ----
        double[] piArr = new double[nVar];
        double[] heArr = new double[nVar];
        double[] hoArr = new double[nVar];
        boolean[] segregating = new boolean[nVar];
        boolean[] valid = new boolean[nVar];
        double[] thetaHArr = new double[nVar];
        boolean[] polarized = new boolean[nVar];

        // Comparative arrays (only allocated if pop2 provided)
        double[] fstNumArr = null, fstDenArr = null, dxyArr = null;
        double[] piW1Arr = null, piW2Arr = null;
        boolean[] validComp = null;

        if (pop2 != null) {
            fstNumArr = new double[nVar];
            fstDenArr = new double[nVar];
            dxyArr = new double[nVar];
            piW1Arr = new double[nVar];
            piW2Arr = new double[nVar];
            validComp = new boolean[nVar];
        }

        for (int i = 0; i < nVar; i++) {
            VariantData v = variants.get(i);
            int ac = v.ac[ctx.index];
            int an = v.an[ctx.index];
            double af = v.af[ctx.index];
            int het = v.het[ctx.index];
            int homAlt = v.homAlt[ctx.index];

            if (an < 2) continue;

            // Apply filter
            if (filter.isActive() && !filter.passes(ac, an, af, het, homAlt, v.node)) {
                continue;
            }

            valid[i] = true;

            piArr[i] = VectorOps.piPerSite(af, an);
            heArr[i] = 2.0 * af * (1.0 - af);
            hoArr[i] = (double) het / (an / 2);
            segregating[i] = ac > 0 && ac < an;

            // Fay & Wu's H
            if (v.ancestralAllele != null) {
                int derivedCount;
                if ("REF".equals(v.ancestralAllele)) {
                    derivedCount = ac;
                } else if ("ALT".equals(v.ancestralAllele)) {
                    derivedCount = an - ac;
                } else {
                    derivedCount = -1;
                }
                if (derivedCount >= 0) {
                    thetaHArr[i] = FayWuH.thetaHPerSite(derivedCount, an);
                    polarized[i] = true;
                }
            }

            if (ctx2 != null) {
                int an2 = v.an[ctx2.index];
                double af2 = v.af[ctx2.index];
                if (an2 >= 2) {
                    validComp[i] = true;
                    double[] fstC = VectorOps.hudsonFstComponents(af, an, af2, an2);
                    fstNumArr[i] = fstC[0];
                    fstDenArr[i] = fstC[1];
                    dxyArr[i] = VectorOps.dxyPerSite(af, af2);
                    piW1Arr[i] = piArr[i];
                    piW2Arr[i] = VectorOps.piPerSite(af2, an2);
                }
            }
        }

        // Phase 1: Pre-compute window boundaries
        long minPos = positions[0];
        long maxPos = positions[nVar - 1];

        List<long[]> windowBounds = new ArrayList<>();
        for (long wStart = minPos; wStart <= maxPos; wStart += stepSize) {
            windowBounds.add(new long[]{wStart, wStart + windowSize - 1});
        }
        int nWindows = windowBounds.size();
        if (nWindows == 0) return Stream.empty();

        // Phase 2: Parallel window summation — pure array ops
        final PopulationContext ctxF = ctx;
        final PopulationContext ctx2F = ctx2;
        final double[] fstNumF = fstNumArr, fstDenF = fstDenArr, dxyF = dxyArr;
        final double[] piW1F = piW1Arr, piW2F = piW2Arr;
        final boolean[] validCompF = validComp;

        WindowResult[] windowResults = new WindowResult[nWindows];

        IntStream.range(0, nWindows).parallel().forEach(w -> {
            long wStart = windowBounds.get(w)[0];
            long wEnd = windowBounds.get(w)[1];

            int lo = Arrays.binarySearch(positions, wStart);
            if (lo < 0) lo = -(lo + 1);

            double piSum = 0.0, heSum = 0.0, hoSum = 0.0;
            long nVariants = 0, nSeg = 0;
            double fstNum = 0.0, fstDen = 0.0, dxySum = 0.0;
            double piW1Sum = 0.0, piW2Sum = 0.0;
            double thetaHSum = 0.0, piPolarSum = 0.0;
            long nPol = 0;

            for (int j = lo; j < nVar; j++) {
                if (positions[j] > wEnd) break;
                if (!valid[j]) continue;

                nVariants++;
                piSum += piArr[j];
                heSum += heArr[j];
                hoSum += hoArr[j];
                if (segregating[j]) nSeg++;

                if (polarized[j]) {
                    thetaHSum += thetaHArr[j];
                    piPolarSum += piArr[j];
                    nPol++;
                }

                if (ctx2F != null && validCompF[j]) {
                    fstNum += fstNumF[j];
                    fstDen += fstDenF[j];
                    dxySum += dxyF[j];
                    piW1Sum += piW1F[j];
                    piW2Sum += piW2F[j];
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
            wr.n_segregating = nSeg;
            wr.pi = piSum / L;
            wr.theta_w = ctxF.a_n > 0 ? (nSeg / ctxF.a_n) / L : 0.0;
            wr.tajima_d = TajimaD.compute(piSum, nSeg, ctxF);
            wr.het_exp = heSum / L;
            wr.het_obs = hoSum / L;
            wr.fis = wr.het_exp > 0 ? 1.0 - wr.het_obs / wr.het_exp : 0.0;

            if (nPol > 0) {
                wr.fay_wu_h = FayWuH.compute(piPolarSum, thetaHSum) / nPol;
            }

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
            if (wr == null) continue;

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
            gwNode.setProperty("fay_wu_h", wr.fay_wu_h);
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
     */
    private static class VariantData {
        final long pos;
        final int[] ac;
        final int[] an;
        final double[] af;
        final int[] het;
        final int[] homAlt;
        final String ancestralAllele;
        final Node node;

        VariantData(long pos, int[] ac, int[] an, double[] af, int[] het, int[] homAlt,
                    String ancestralAllele, Node node) {
            this.pos = pos;
            this.ac = ac;
            this.an = an;
            this.af = af;
            this.het = het;
            this.homAlt = homAlt;
            this.ancestralAllele = ancestralAllele;
            this.node = node;
        }
    }
}
