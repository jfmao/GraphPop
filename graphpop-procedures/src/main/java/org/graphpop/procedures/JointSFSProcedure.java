package org.graphpop.procedures;

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
 * Computes the 2D joint site frequency spectrum between two populations.
 *
 * <p>FAST PATH: uses only allele-count arrays on Variant nodes.</p>
 *
 * <p>Usage:
 * <pre>
 * CALL graphpop.joint_sfs('chr22', 16000000, 17000000, 'AFR', 'EUR', false)
 * YIELD joint_sfs, n_variants, max_ac1, max_ac2
 * </pre>
 * </p>
 */
public class JointSFSProcedure {

    @Context
    public Transaction tx;

    public static class JointSFSResult {
        /** Flattened 2D histogram: row-major [ac1][ac2], dimensions (max_ac1+1) x (max_ac2+1). */
        public List<Long> joint_sfs;
        public long n_variants;
        public long max_ac1;
        public long max_ac2;
        public long dim1;
        public long dim2;

        public JointSFSResult() {}
    }

    @Procedure(name = "graphpop.joint_sfs", mode = Mode.READ)
    public Stream<JointSFSResult> jointSfs(
            @Name("chr") String chr,
            @Name("start") long start,
            @Name("end") long end,
            @Name("pop1") String pop1,
            @Name("pop2") String pop2,
            @Name(value = "folded", defaultValue = "false") boolean folded
    ) {
        int popIndex1 = -1;
        int popIndex2 = -1;
        int maxAn1 = 0;
        int maxAn2 = 0;

        List<int[]> entries = new ArrayList<>(); // [ac1, an1, ac2, an2]

        var result = tx.execute(
                "MATCH (v:Variant) WHERE v.chr = $chr AND v.pos >= $start AND v.pos <= $end RETURN v",
                Map.of("chr", chr, "start", start, "end", end)
        );

        try {
            while (result.hasNext()) {
                Map<String, Object> row = result.next();
                Node variant = (Node) row.get("v");

                if (popIndex1 < 0) {
                    String[] popIds = ArrayUtil.toStringArray(variant.getProperty("pop_ids"));
                    for (int i = 0; i < popIds.length; i++) {
                        if (popIds[i].equals(pop1)) popIndex1 = i;
                        if (popIds[i].equals(pop2)) popIndex2 = i;
                    }
                    if (popIndex1 < 0) throw new RuntimeException("Population not found: " + pop1);
                    if (popIndex2 < 0) throw new RuntimeException("Population not found: " + pop2);
                }

                int[] acArr = ArrayUtil.toIntArray(variant.getProperty("ac"));
                int[] anArr = ArrayUtil.toIntArray(variant.getProperty("an"));

                int ac1 = acArr[popIndex1], an1 = anArr[popIndex1];
                int ac2 = acArr[popIndex2], an2 = anArr[popIndex2];

                if (an1 < 2 || an2 < 2) continue;
                if (an1 > maxAn1) maxAn1 = an1;
                if (an2 > maxAn2) maxAn2 = an2;

                entries.add(new int[]{ac1, an1, ac2, an2});
            }
        } finally {
            result.close();
        }

        int dim1, dim2;
        if (folded) {
            dim1 = maxAn1 / 2 + 1;
            dim2 = maxAn2 / 2 + 1;
        } else {
            dim1 = maxAn1 + 1;
            dim2 = maxAn2 + 1;
        }

        long[] histogram = new long[dim1 * dim2];

        for (int[] e : entries) {
            int idx1 = folded ? Math.min(e[0], e[1] - e[0]) : e[0];
            int idx2 = folded ? Math.min(e[2], e[3] - e[2]) : e[2];
            if (idx1 >= 0 && idx1 < dim1 && idx2 >= 0 && idx2 < dim2) {
                histogram[idx1 * dim2 + idx2]++;
            }
        }

        List<Long> flat = new ArrayList<>(histogram.length);
        for (long v : histogram) flat.add(v);

        JointSFSResult out = new JointSFSResult();
        out.joint_sfs = flat;
        out.n_variants = entries.size();
        out.max_ac1 = maxAn1;
        out.max_ac2 = maxAn2;
        out.dim1 = dim1;
        out.dim2 = dim2;

        return Stream.of(out);
    }
}
