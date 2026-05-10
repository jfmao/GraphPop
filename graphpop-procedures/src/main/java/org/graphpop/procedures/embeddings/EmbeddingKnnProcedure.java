package org.graphpop.procedures.embeddings;

import org.neo4j.graphdb.Result;
import org.neo4j.graphdb.Transaction;
import org.neo4j.procedure.Context;
import org.neo4j.procedure.Description;
import org.neo4j.procedure.Mode;
import org.neo4j.procedure.Name;
import org.neo4j.procedure.Procedure;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

/**
 * Top-k nearest-neighbour search over {@code :Sample.embedding}
 * vectors by cosine similarity (M10).
 *
 * <pre>
 * CALL graphpop.embedding.knn($sampleId, $k)
 *   YIELD query_sample_id, neighbor_sample_id, cosine_similarity, rank
 * </pre>
 *
 * <p>Loads every embedded sample's vector into memory in a single
 * Cypher pass, computes pairwise cosine similarity against the
 * query sample, and emits the top-k by descending similarity.
 * The query sample itself is excluded from the results.</p>
 *
 * <p>Linear scan rather than vector index — keeps the test harness
 * uniform and avoids requiring Neo4j 5.13+ for the current
 * benchmark scale (n_samples ≲ 100k). Vector-index acceleration is
 * a drop-in replacement when a downstream needs it.</p>
 */
public class EmbeddingKnnProcedure {

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.embedding.knn", mode = Mode.READ)
    @Description("Top-k nearest neighbours of a sample by cosine "
            + "similarity over :Sample.embedding vectors.")
    public Stream<EmbeddingKnnResult> knn(
            @Name("sampleId") String sampleId,
            @Name("k") long k
    ) {
        if (k <= 0) return Stream.empty();
        EmbeddingTable table = EmbeddingTable.load(tx);
        if (table.size() < 2) return Stream.empty();

        double[] queryVec = table.vectorOf(sampleId);
        if (queryVec == null) return Stream.empty();
        double queryNorm = norm(queryVec);
        if (queryNorm == 0.0) return Stream.empty();

        List<double[]> scored = new ArrayList<>(table.size());
        for (int i = 0; i < table.size(); i++) {
            String sid = table.sampleId(i);
            if (sid.equals(sampleId)) continue;
            double[] v = table.vector(i);
            double n = norm(v);
            if (n == 0.0) continue;
            double sim = dot(queryVec, v) / (queryNorm * n);
            scored.add(new double[]{i, sim});
        }
        scored.sort((a, b) -> Double.compare(b[1], a[1]));

        int kInt = (int) Math.min(k, scored.size());
        List<EmbeddingKnnResult> rows = new ArrayList<>(kInt);
        for (int rank = 0; rank < kInt; rank++) {
            double[] row = scored.get(rank);
            int idx = (int) row[0];
            rows.add(new EmbeddingKnnResult(
                    sampleId, table.sampleId(idx), row[1], rank + 1L));
        }
        return rows.stream();
    }

    static double norm(double[] v) {
        double s = 0.0;
        for (double x : v) s += x * x;
        return Math.sqrt(s);
    }

    static double dot(double[] a, double[] b) {
        int n = Math.min(a.length, b.length);
        double s = 0.0;
        for (int i = 0; i < n; i++) s += a[i] * b[i];
        return s;
    }

    /**
     * Packed in-memory table of per-sample embedding vectors loaded
     * via Cypher. Shared with {@code EmbeddingClusterProcedure}.
     */
    static final class EmbeddingTable {
        private final List<String> ids = new ArrayList<>();
        private final List<double[]> vecs = new ArrayList<>();

        static EmbeddingTable load(Transaction tx) {
            EmbeddingTable t = new EmbeddingTable();
            Result r = tx.execute(
                    "MATCH (s:Sample) WHERE s.embedding IS NOT NULL "
                  + "RETURN s.sampleId AS sid, s.embedding AS emb");
            try {
                while (r.hasNext()) {
                    Map<String, Object> row = r.next();
                    double[] v = coerceToDoubleArray(row.get("emb"));
                    if (v == null) continue;
                    t.ids.add((String) row.get("sid"));
                    t.vecs.add(v);
                }
            } finally {
                r.close();
            }
            return t;
        }

        /**
         * Neo4j stores list properties as primitive arrays internally
         * and may return them as {@code double[]}, {@code float[]},
         * {@code long[]}, {@code Object[]}, or {@link List}. Coerce
         * any of these to a {@code double[]}; return {@code null}
         * for unsupported / null inputs.
         */
        static double[] coerceToDoubleArray(Object v) {
            if (v == null) return null;
            if (v instanceof double[]) return (double[]) v;
            if (v instanceof float[]) {
                float[] f = (float[]) v;
                double[] d = new double[f.length];
                for (int i = 0; i < f.length; i++) d[i] = f[i];
                return d;
            }
            if (v instanceof long[]) {
                long[] l = (long[]) v;
                double[] d = new double[l.length];
                for (int i = 0; i < l.length; i++) d[i] = l[i];
                return d;
            }
            if (v instanceof int[]) {
                int[] a = (int[]) v;
                double[] d = new double[a.length];
                for (int i = 0; i < a.length; i++) d[i] = a[i];
                return d;
            }
            if (v instanceof Object[]) {
                Object[] o = (Object[]) v;
                double[] d = new double[o.length];
                for (int i = 0; i < o.length; i++) {
                    d[i] = ((Number) o[i]).doubleValue();
                }
                return d;
            }
            if (v instanceof List) {
                List<?> list = (List<?>) v;
                double[] d = new double[list.size()];
                for (int i = 0; i < d.length; i++) {
                    d[i] = ((Number) list.get(i)).doubleValue();
                }
                return d;
            }
            return null;
        }

        int size() { return ids.size(); }
        String sampleId(int i) { return ids.get(i); }
        double[] vector(int i) { return vecs.get(i); }

        double[] vectorOf(String sampleId) {
            for (int i = 0; i < ids.size(); i++) {
                if (ids.get(i).equals(sampleId)) return vecs.get(i);
            }
            return null;
        }

        List<String> sampleIds() { return Collections.unmodifiableList(ids); }
        List<double[]> vectors() { return Collections.unmodifiableList(vecs); }
    }
}
