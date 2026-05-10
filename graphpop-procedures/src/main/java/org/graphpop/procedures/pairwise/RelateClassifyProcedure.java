package org.graphpop.procedures.pairwise;

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
 * Per-pair relationship classification (M5 Part A).
 *
 * <p>Aggregates the existing pairwise signal (KING-robust kinship +
 * IBS0 fraction, or IBD-derived phi^IBD) and assigns a Manichaikul
 * et al. 2010 relationship label. Writes
 * {@code (:Sample)-[:RELATIVE {relationship, degree, phi, ibs0_frac,
 * source, created_at}]->(:Sample)} edges idempotently per
 * {@code source}.</p>
 *
 * <pre>
 * CALL graphpop.relate.classify('king', {min_phi: 0.0442})
 *   YIELD sample_a, sample_b, relationship, degree, phi, ibs0_frac, source
 * </pre>
 *
 * <p>Sources:</p>
 * <ul>
 *   <li>{@code "king"} — pulls from persisted
 *       {@code (:Sample)-[r:KINSHIP {method:'king-robust'}]->(:Sample)}
 *       (when present); otherwise the user is expected to have run
 *       {@code graphpop.kinship.king} ahead of time and persisted the
 *       output via the user's pipeline.</li>
 *   <li>any IBD source string — aggregates {@code :IBD_SEGMENT} edges
 *       with that {@code source} via the same logic as
 *       {@link IbdKinshipProcedure}.</li>
 * </ul>
 */
public class RelateClassifyProcedure {

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.relate.classify", mode = Mode.WRITE)
    @Description("Per-pair relationship classification from kinship + IBS0 "
            + "(Manichaikul 2010). Writes :RELATIVE edges idempotently per "
            + "source. Sources: 'king' or any IBD source name.")
    @SuppressWarnings("unchecked")
    public Stream<RelativeResult> classify(
            @Name("source") String source,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        if (options == null) options = new HashMap<>();
        double minPhi = getDouble(options, "min_phi", 0.0442);
        boolean persist = getBoolean(options, "persist", true);

        RelativeClassifier.Cutoffs cutoffs = new RelativeClassifier.Cutoffs(
                getDouble(options, "identical", 0.354),
                getDouble(options, "first_degree", 0.177),
                getDouble(options, "second_degree", 0.0884),
                getDouble(options, "third_degree", 0.0442),
                getDouble(options, "ibs0_threshold", 0.0050));

        // Pull pair signal from the requested source.
        List<PairSignal> pairs = "king".equalsIgnoreCase(source)
                ? loadKing(tx)
                : loadFromIbd(tx, source);

        // Idempotent: clear existing :RELATIVE for this source.
        if (persist) {
            tx.execute(
                "MATCH ()-[r:RELATIVE {source: $source}]->() DELETE r",
                Map.of("source", source));
        }

        long now = System.currentTimeMillis();
        List<Map<String, Object>> writeBatch = new ArrayList<>();
        List<RelativeResult> rows = new ArrayList<>(pairs.size());

        for (PairSignal p : pairs) {
            if (p.phi < minPhi) continue;
            RelativeClassifier.Verdict v = RelativeClassifier.classify(
                    p.phi, p.ibs0Frac, cutoffs);
            rows.add(new RelativeResult(p.a, p.b,
                    v.relationship, v.degree, p.phi, p.ibs0Frac, source));

            if (persist) {
                Map<String, Object> r = new HashMap<>();
                r.put("a", p.a);
                r.put("b", p.b);
                r.put("relationship", v.relationship);
                r.put("degree", (long) v.degree);
                r.put("phi", p.phi);
                r.put("ibs0_frac", p.ibs0Frac);
                writeBatch.add(r);
                if (writeBatch.size() >= 10_000) {
                    flushBatch(tx, source, now, writeBatch);
                    writeBatch.clear();
                }
            }
        }
        if (persist && !writeBatch.isEmpty()) {
            flushBatch(tx, source, now, writeBatch);
        }
        return rows.stream();
    }

    private static void flushBatch(Transaction tx, String source, long now,
                                    List<Map<String, Object>> batch) {
        tx.execute(
                "UNWIND $rows AS r "
              + "MATCH (a:Sample {sampleId: r.a}), (b:Sample {sampleId: r.b}) "
              + "CREATE (a)-[:RELATIVE {"
              + "  relationship: r.relationship, degree: r.degree, "
              + "  phi: r.phi, ibs0_frac: r.ibs0_frac, "
              + "  source: $source, "
              + "  created_at: datetime({epochMillis: $now})"
              + "}]->(b)",
                Map.of("rows", batch, "source", source, "now", now));
    }

    /** Pairwise signal harvested from a source. */
    private static final class PairSignal {
        final String a;
        final String b;
        final double phi;
        final double ibs0Frac;
        PairSignal(String a, String b, double phi, double ibs0Frac) {
            this.a = a; this.b = b; this.phi = phi; this.ibs0Frac = ibs0Frac;
        }
    }

    private static List<PairSignal> loadKing(Transaction tx) {
        // Pull from any KINSHIP relationship tagged method='king-robust'.
        // The user is expected to have run kinship.king and persisted via
        // their own pipeline; we don't inline-rerun here in v1.
        List<PairSignal> out = new ArrayList<>();
        Result r = tx.execute(
                "MATCH (a:Sample)-[k:KINSHIP {method: 'king-robust'}]->(b:Sample) "
              + "RETURN a.sampleId AS a, b.sampleId AS b, "
              + "k.phi AS phi, k.ibs0_frac AS ibs0_frac");
        try {
            while (r.hasNext()) {
                Map<String, Object> row = r.next();
                String a = (String) row.get("a");
                String b = (String) row.get("b");
                double phi = ((Number) row.get("phi")).doubleValue();
                Object i = row.get("ibs0_frac");
                double ibs0 = (i == null) ? Double.NaN : ((Number) i).doubleValue();
                out.add(new PairSignal(a, b, phi, ibs0));
            }
        } finally {
            r.close();
        }
        return out;
    }

    private static List<PairSignal> loadFromIbd(Transaction tx, String source) {
        // Aggregate IBD segments per ordered pair; phi = total length / (2*genome).
        // Genome length inferred from per-chromosome span (matches IbdKinshipProcedure).
        double totalGenome = 0.0;
        Result spans = tx.execute(
                "MATCH ()-[r:IBD_SEGMENT {source: $source}]->() "
              + "WITH r.chr AS chr, max(r.end) AS hi, min(r.start) AS lo "
              + "RETURN sum(hi - lo) AS total_bp",
                Map.of("source", source));
        try {
            if (spans.hasNext()) {
                Object v = spans.next().get("total_bp");
                if (v instanceof Number) totalGenome = ((Number) v).doubleValue();
            }
        } finally {
            spans.close();
        }
        if (totalGenome <= 0) return List.of();
        double denom = 2.0 * totalGenome;

        List<PairSignal> out = new ArrayList<>();
        Result r = tx.execute(
                "MATCH (a:Sample)-[r:IBD_SEGMENT {source: $source}]->(b:Sample) "
              + "WITH a.sampleId AS sa, b.sampleId AS sb, "
              + "sum(r.length_bp) AS total_bp "
              + "RETURN sa, sb, total_bp",
                Map.of("source", source));
        try {
            while (r.hasNext()) {
                Map<String, Object> row = r.next();
                String sa = (String) row.get("sa");
                String sb = (String) row.get("sb");
                double totalBp = ((Number) row.get("total_bp")).doubleValue();
                double phi = totalBp / denom;
                out.add(new PairSignal(sa, sb, phi, Double.NaN));
            }
        } finally {
            r.close();
        }
        return out;
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
