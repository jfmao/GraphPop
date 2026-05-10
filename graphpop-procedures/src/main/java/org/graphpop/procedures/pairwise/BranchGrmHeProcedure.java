package org.graphpop.procedures.pairwise;

import org.neo4j.graphdb.Result;
import org.neo4j.graphdb.Transaction;
import org.neo4j.procedure.Context;
import org.neo4j.procedure.Description;
import org.neo4j.procedure.Mode;
import org.neo4j.procedure.Name;
import org.neo4j.procedure.Procedure;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Stream;

/**
 * Haseman–Elston regression heritability via Algorithm V mat-vec
 * products + Hutchinson trace estimator (Pazokitoroudi et al. 2020
 * RHE-mc, single-component formulation):
 *
 * <pre>
 *   num     = y_c' · G · y_c                       (one mat-vec)
 *   tr(G²)  ≈ (1/M) Σ_k ||G · u_k||², u_k Rademacher (M mat-vecs)
 *   ĥ²      = num / tr(G²)
 *   SE(ĥ²)  ≈ √(2 / tr(G²))                        (Wald)
 * </pre>
 *
 * <p>{@code G} is the centred branch GRM, so {@code tr(G) = 0} and
 * the second-order correction vanishes. Conditional predicates from
 * M4.1 step 4 compose unchanged.</p>
 */
public class BranchGrmHeProcedure {

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.kinship.branch_grm_he", mode = Mode.READ)
    @Description("Haseman-Elston heritability via Algorithm V + "
            + "Hutchinson trace. Single-component RHE-mc formulation. "
            + "Conditional predicates compose.")
    @SuppressWarnings("unchecked")
    public Stream<HeResult> branchGrmHe(
            @Name("run_id") String runId,
            @Name("phenotype") List<Double> phenotype,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        if (options == null) options = new HashMap<>();
        long start = getLong(options, "start", 0L);
        long end = getLong(options, "end", Long.MAX_VALUE);
        long seed = getLong(options, "seed", 42L);
        int nHutch = (int) getLong(options, "n_hutchinson", 50L);
        if (nHutch < 1) {
            throw new RuntimeException("n_hutchinson must be >= 1; got " + nHutch);
        }

        ARG arg = ARGTraversal.load(tx, runId, start, end);
        if (arg.nNodes == 0 || arg.nSamples() == 0) return Stream.empty();
        int n = arg.nSamples();
        if (phenotype.size() != n) {
            throw new RuntimeException(
                "phenotype length " + phenotype.size() + " != n_samples " + n);
        }

        BranchWeightFn weight =
                BranchGrmConditioning.fromOptions(tx, runId, options);

        // Centre y.
        double[] y = new double[n];
        double sum = 0;
        for (int i = 0; i < n; i++) { y[i] = phenotype.get(i); sum += y[i]; }
        double mean = sum / n;
        for (int i = 0; i < n; i++) y[i] -= mean;

        // num = y' · G · y.
        double num;
        double yNormSq = 0;
        for (double v : y) yNormSq += v * v;
        if (yNormSq < 1e-30) {
            // Constant phenotype: no signal.
            return Stream.of(new HeResult(0.0, Double.NaN, 0.0, 0.0,
                    n, nHutch, "he_grm"));
        }
        double[] Gy = BranchGrmMatVec.apply(arg, start, end, weight, y);
        num = 0;
        for (int i = 0; i < n; i++) num += y[i] * Gy[i];

        // tr(G²) ≈ (1/M) Σ ||G·u||² with Rademacher u.
        Random rng = new Random(seed);
        double trGsq = 0;
        double[] u = new double[n];
        for (int k = 0; k < nHutch; k++) {
            for (int i = 0; i < n; i++) {
                u[i] = rng.nextBoolean() ? 1.0 : -1.0;
            }
            double[] Gu = BranchGrmMatVec.apply(arg, start, end, weight, u);
            double s = 0;
            for (int i = 0; i < n; i++) s += Gu[i] * Gu[i];
            trGsq += s;
        }
        trGsq /= nHutch;

        if (trGsq <= 0) {
            return Stream.of(new HeResult(0.0, Double.NaN, num, trGsq,
                    n, nHutch, "he_grm"));
        }

        double h2 = num / trGsq;
        double se = Math.sqrt(2.0 / trGsq);
        return Stream.of(new HeResult(
                h2, se, num, trGsq, n, nHutch, "he_grm"));
    }

    private static long getLong(Map<String, Object> opts, String key, long def) {
        Object v = opts.get(key);
        return (v instanceof Number) ? ((Number) v).longValue() : def;
    }
}
