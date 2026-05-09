package org.graphpop.procedures.pairwise;

import org.graphpop.procedures.GenotypeLoader;
import org.graphpop.procedures.PackedGenotypeReader;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Result;
import org.neo4j.graphdb.Transaction;
import org.neo4j.procedure.Context;
import org.neo4j.procedure.Description;
import org.neo4j.procedure.Mode;
import org.neo4j.procedure.Name;
import org.neo4j.procedure.Procedure;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

/**
 * Identity-by-state (IBS) pairwise statistic on packed genotypes
 * (M4.2). Sibling of {@link KingRobustProcedure}; matrix-only (no
 * ARG required).
 *
 * <p>{@code IBS_ij ∈ [0, 1]}; identical samples = 1, opposite
 * homozygotes average toward 0. Output uses the shared
 * {@link KinshipResult} POJO with {@code phi} = IBS,
 * {@code ibs0} = opposite-homozygote count,
 * {@code het_het} = identical-genotype count, {@code n_snp} =
 * variants used, {@code method} = {@code "ibs"}.</p>
 *
 * <pre>
 * CALL graphpop.kinship.ibs('chr22', 'EUR', {min_snp: 1000})
 *   YIELD sample_a, sample_b, phi, ibs0, het_het, n_snp, n_aa_min, method
 * </pre>
 */
public class IbsProcedure {

    private static final int DEFAULT_MIN_SNP = 1000;
    private static final double DEFAULT_MIN_IBS = 0.0;

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.kinship.ibs", mode = Mode.READ)
    @Description("Identity-by-state (IBS) pairwise statistic on packed "
            + "genotypes. Matrix-only sibling of kinship.king; no ARG needed.")
    @SuppressWarnings("unchecked")
    public Stream<KinshipResult> ibs(
            @Name("chr") String chr,
            @Name("pop") String pop,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        if (options == null) options = new HashMap<>();

        long start = getLong(options, "start", 1L);
        long end = getLong(options, "end", Long.MAX_VALUE);
        int minSnp = (int) getLong(options, "min_snp", DEFAULT_MIN_SNP);
        double minIbs = getDouble(options, "min_ibs", DEFAULT_MIN_IBS);
        boolean includeSelf = getBoolean(options, "include_self", false);
        List<String> sampleList = (List<String>) options.get("samples");

        Map<String, Integer> sampleIndex = (sampleList != null && !sampleList.isEmpty())
                ? GenotypeLoader.buildSampleIndex(tx, sampleList)
                : GenotypeLoader.buildSampleIndex(tx, pop);
        if (sampleIndex.isEmpty()) return Stream.empty();

        int nSamples = sampleIndex.size();
        String[] sampleIds = new String[nSamples];
        for (Map.Entry<String, Integer> e : sampleIndex.entrySet())
            sampleIds[e.getValue()] = e.getKey();
        int[] packedIndices = GenotypeLoader.buildPackedIndices(tx, sampleIndex);

        // Per-pair accumulators.
        int nPairs = nSamples * (nSamples + 1) / 2;
        IbsComputer.PairAccumulator[] acc = new IbsComputer.PairAccumulator[nPairs];
        for (int p = 0; p < nPairs; p++) acc[p] = new IbsComputer.PairAccumulator();
        int[] gtBuf = new int[nSamples];

        Result r = tx.execute(
                "MATCH (v:Variant) "
              + "WHERE v.chr = $chr AND v.pos >= $start AND v.pos <= $end "
              + "RETURN v ORDER BY v.pos",
                Map.of("chr", chr, "start", start, "end", end));
        try {
            while (r.hasNext()) {
                Node variant = (Node) r.next().get("v");
                Object gtPackedObj = variant.getProperty("gt_packed", null);
                if (!(gtPackedObj instanceof byte[])) continue;
                byte[] gtPacked = (byte[]) gtPackedObj;
                for (int s = 0; s < nSamples; s++) {
                    int pi = packedIndices[s];
                    gtBuf[s] = (pi >= 0)
                            ? PackedGenotypeReader.genotype(gtPacked, pi)
                            : PackedGenotypeReader.GT_MISSING;
                }
                int p = 0;
                for (int b = 0; b < nSamples; b++)
                    for (int a = 0; a <= b; a++)
                        acc[p++].accumulate(gtBuf[a], gtBuf[b]);
            }
        } finally {
            r.close();
        }

        List<KinshipResult> rows = new ArrayList<>();
        int p = 0;
        for (int b = 0; b < nSamples; b++) {
            for (int a = 0; a <= b; a++, p++) {
                if (a == b && !includeSelf) continue;
                IbsComputer.PairAccumulator s = acc[p];
                if (s.nUsed < minSnp) continue;
                double ibs = s.ibs();
                if (Double.isNaN(ibs)) continue;
                if (ibs < minIbs) continue;
                rows.add(new KinshipResult(
                        sampleIds[a], sampleIds[b],
                        ibs, s.ibs0, s.ibs2, s.nUsed, 0L, "ibs"));
            }
        }
        return rows.stream();
    }

    private static long getLong(Map<String, Object> opts, String key, long def) {
        Object v = opts.get(key);
        return (v instanceof Number) ? ((Number) v).longValue() : def;
    }

    private static double getDouble(Map<String, Object> opts, String key, double def) {
        Object v = opts.get(key);
        return (v instanceof Number) ? ((Number) v).doubleValue() : def;
    }

    private static boolean getBoolean(Map<String, Object> opts, String key, boolean def) {
        Object v = opts.get(key);
        return (v instanceof Boolean) ? (Boolean) v : def;
    }
}
