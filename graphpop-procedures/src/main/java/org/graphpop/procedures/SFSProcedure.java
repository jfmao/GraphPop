package org.graphpop.procedures;

import org.neo4j.graphdb.Label;
import org.neo4j.graphdb.Node;
import org.neo4j.graphdb.Transaction;
import org.neo4j.procedure.Context;
import org.neo4j.procedure.Mode;
import org.neo4j.procedure.Name;
import org.neo4j.procedure.Procedure;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

/**
 * Computes the Site Frequency Spectrum (SFS) for a population in a genomic region.
 *
 * <p>FAST PATH: uses only allele-count arrays on Variant nodes.</p>
 *
 * <p>Usage:
 * <pre>
 * CALL graphpop.sfs('chr22', 16000000, 17000000, 'AFR', false)
 * YIELD sfs, n_variants, max_ac
 * </pre>
 * </p>
 */
public class SFSProcedure {

    @Context
    public Transaction tx;

    public static class SFSResult {
        public List<Long> sfs;
        public long n_variants;
        public long max_ac;

        public SFSResult() {}
    }

    @Procedure(name = "graphpop.sfs", mode = Mode.READ)
    public Stream<SFSResult> sfs(
            @Name("chr") String chr,
            @Name("start") long start,
            @Name("end") long end,
            @Name("pop") String pop,
            @Name(value = "folded", defaultValue = "false") boolean folded
    ) {
        int popIndex = -1;
        int maxAn = 0;

        // First pass: determine max AN for this population to size the histogram
        // We'll collect in a dynamically-sized structure
        List<int[]> acAnPairs = new ArrayList<>();

        var result = tx.execute(
                "MATCH (v:Variant) WHERE v.chr = $chr AND v.pos >= $start AND v.pos <= $end RETURN v",
                Map.of("chr", chr, "start", start, "end", end)
        );

        try {
            while (result.hasNext()) {
                Map<String, Object> row = result.next();
                Node variant = (Node) row.get("v");

                if (popIndex < 0) {
                    String[] popIds = ArrayUtil.toStringArray(variant.getProperty("pop_ids"));
                    for (int i = 0; i < popIds.length; i++) {
                        if (popIds[i].equals(pop)) { popIndex = i; break; }
                    }
                    if (popIndex < 0) throw new RuntimeException("Population not found: " + pop);
                }

                int[] acArr = ArrayUtil.toIntArray(variant.getProperty("ac"));
                int[] anArr = ArrayUtil.toIntArray(variant.getProperty("an"));

                int ac = acArr[popIndex];
                int an = anArr[popIndex];

                if (an < 2) continue;
                if (an > maxAn) maxAn = an;

                acAnPairs.add(new int[]{ac, an});
            }
        } finally {
            result.close();
        }

        // Build histogram
        long[] histogram;
        if (folded) {
            int foldedSize = maxAn / 2 + 1;
            histogram = new long[foldedSize];
            for (int[] pair : acAnPairs) {
                int ac = pair[0];
                int an = pair[1];
                int foldedAc = Math.min(ac, an - ac);
                if (foldedAc >= 0 && foldedAc < histogram.length) {
                    histogram[foldedAc]++;
                }
            }
        } else {
            histogram = new long[maxAn + 1];
            for (int[] pair : acAnPairs) {
                int ac = pair[0];
                if (ac >= 0 && ac < histogram.length) {
                    histogram[ac]++;
                }
            }
        }

        // Convert to List<Long> for Neo4j compatibility
        List<Long> sfsList = new ArrayList<>(histogram.length);
        for (long count : histogram) sfsList.add(count);

        SFSResult out = new SFSResult();
        out.sfs = sfsList;
        out.n_variants = acAnPairs.size();
        out.max_ac = maxAn;

        return Stream.of(out);
    }
}
