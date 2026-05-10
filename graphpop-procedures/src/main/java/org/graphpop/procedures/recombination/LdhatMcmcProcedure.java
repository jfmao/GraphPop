package org.graphpop.procedures.recombination;

import org.graphpop.procedures.selection.SelectionUtils;
import org.neo4j.graphdb.Result;
import org.neo4j.graphdb.Transaction;
import org.neo4j.procedure.Context;
import org.neo4j.procedure.Description;
import org.neo4j.procedure.Mode;
import org.neo4j.procedure.Name;
import org.neo4j.procedure.Procedure;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Stream;

/**
 * Per-window LDhat-style MCMC posterior on ρ (M13.A).
 *
 * <p>For each sliding window, runs a Metropolis-Hastings sampler
 * on {@code log(ρ)}:</p>
 *
 * <ul>
 *   <li><b>Prior</b>: {@code log10(ρ) ~ Uniform(prior_log_lo,
 *       prior_log_hi)} (default {@code [-12, -2]}).</li>
 *   <li><b>Likelihood</b>: Gaussian approximation on per-pair r²
 *       residuals around Hudson 1985's {@code E[r²|n, ρ·d]}:
 *       <pre>log L ∝ −(1/2σ²) Σ_pairs (r² − E[r²|ρ·d_pair])²</pre></li>
 *   <li><b>Proposal</b>: {@code log10(ρ_new) = log10(ρ_old)
 *       + Normal(0, prop_sd)} (default {@code prop_sd = 0.5}).</li>
 * </ul>
 *
 * <p>Output per window: posterior mean ρ plus 2.5 % / 97.5 %
 * quantiles, sample count, and acceptance rate. Deterministic
 * given {@code seed}.</p>
 *
 * <pre>
 * CALL graphpop.recombination.ldhat_mcmc($sampleIds,
 *         {window_size: 10000, step: 5000, n_iter: 5000,
 *          burn_in: 1000, prop_sd: 0.5, sigma: 0.1, seed: 42})
 *   YIELD start, end, n_variant_pairs, rho_posterior_mean,
 *         rho_lower_2_5, rho_upper_97_5, n_iter, n_accepted,
 *         n_samples, runId
 * </pre>
 *
 * <p>v1 deviates from the original LDhat (McVean 2002) reversible-
 * jump MCMC in two ways: (1) per-window Bayesian estimate rather
 * than a single piecewise-constant ρ-map with variable breakpoints,
 * and (2) Hudson Gaussian-approximation likelihood rather than the
 * LDhat coalescent likelihood lookup. The two simplifications give
 * a tractable v1; the full reversible-jump variant is a v2
 * follow-up.</p>
 */
public class LdhatMcmcProcedure {

    @Context
    public Transaction tx;

    @Procedure(name = "graphpop.recombination.ldhat_mcmc", mode = Mode.READ)
    @Description("Per-window LDhat-style MCMC posterior on ρ. "
            + "Metropolis-Hastings on log10(ρ) with Hudson 1985 "
            + "Gaussian-likelihood approximation.")
    public Stream<LdhatMcmcResult> ldhatMcmc(
            @Name("sampleIds") List<String> sampleIds,
            @Name(value = "options", defaultValue = "{}") Map<String, Object> options
    ) {
        if (options == null) options = new HashMap<>();
        if (sampleIds == null || sampleIds.isEmpty()) return Stream.empty();
        long windowSize = getLong(options, "window_size", 10_000L);
        long step = getLong(options, "step", windowSize);
        double minMaf = getDouble(options, "min_maf", 0.05);
        long maxPairDistance = getLong(options, "max_pair_distance", 5_000L);
        int nIter = (int) getLong(options, "n_iter", 5_000L);
        int burnIn = (int) getLong(options, "burn_in", 1_000L);
        double propSd = getDouble(options, "prop_sd", 0.5);
        double sigma = getDouble(options, "sigma", 0.1);
        double priorLogLo = getDouble(options, "prior_log_lo", -12.0);
        double priorLogHi = getDouble(options, "prior_log_hi", -2.0);
        long seed = getLong(options, "seed", 42L);

        if (windowSize <= 0 || step <= 0 || maxPairDistance <= 0
                || nIter <= 0 || burnIn < 0) {
            throw new IllegalArgumentException(
                    "window_size, step, max_pair_distance, n_iter must be "
                  + "positive; burn_in must be ≥ 0");
        }
        if (priorLogHi <= priorLogLo) {
            throw new IllegalArgumentException(
                    "prior_log_hi must be > prior_log_lo");
        }

        long seqEnd = maxVariantPosition(tx);
        if (seqEnd <= 0) return Stream.empty();
        long[][] windows = SelectionUtils.slidingWindows(
                seqEnd + 1, windowSize, step);
        if (windows.length == 0) return Stream.empty();

        String runId = peekRunId(tx);

        List<LdhatMcmcResult> rows = new ArrayList<>(windows.length);
        Random rng = new Random(seed);
        for (long[] win : windows) {
            long lo = win[0];
            long hi = win[1];
            List<LdPairLoader.Pair> pairs = LdPairLoader.loadPairs(
                    tx, sampleIds, lo, hi, minMaf, maxPairDistance);
            if (pairs.isEmpty()) {
                rows.add(new LdhatMcmcResult(
                        lo, hi, 0L, Double.NaN, Double.NaN, Double.NaN,
                        nIter, 0L, sampleIds.size(), runId));
                continue;
            }
            McmcSummary s = runMcmcForWindow(pairs, sigma, propSd,
                    priorLogLo, priorLogHi, nIter, burnIn, rng);
            rows.add(new LdhatMcmcResult(
                    lo, hi, pairs.size(),
                    s.posteriorMean, s.lower2_5, s.upper97_5,
                    nIter, s.nAccepted, sampleIds.size(), runId));
        }
        return rows.stream();
    }

    static final class McmcSummary {
        double posteriorMean;
        double lower2_5;
        double upper97_5;
        long nAccepted;
    }

    static McmcSummary runMcmcForWindow(List<LdPairLoader.Pair> pairs,
                                          double sigma, double propSd,
                                          double priorLogLo, double priorLogHi,
                                          int nIter, int burnIn, Random rng) {
        // Initialise log10(ρ) at the centre of the prior.
        double logRho = 0.5 * (priorLogLo + priorLogHi);
        double currLogL = logLikelihood(pairs, Math.pow(10.0, logRho), sigma);
        double[] samples = new double[nIter];
        long nAccepted = 0L;
        double invSigma2 = 1.0 / (sigma * sigma);
        for (int it = 0; it < burnIn + nIter; it++) {
            double prop = logRho + rng.nextGaussian() * propSd;
            if (prop < priorLogLo || prop > priorLogHi) {
                if (it >= burnIn) samples[it - burnIn] = logRho;
                continue;
            }
            double propLogL = logLikelihood(
                    pairs, Math.pow(10.0, prop), sigma);
            double logRatio = propLogL - currLogL;
            if (logRatio >= 0.0 || rng.nextDouble() < Math.exp(logRatio)) {
                logRho = prop;
                currLogL = propLogL;
                if (it >= burnIn) nAccepted++;
            }
            if (it >= burnIn) samples[it - burnIn] = logRho;
        }
        double[] rhoSamples = new double[samples.length];
        for (int i = 0; i < samples.length; i++) {
            rhoSamples[i] = Math.pow(10.0, samples[i]);
        }
        Arrays.sort(rhoSamples);
        McmcSummary s = new McmcSummary();
        s.posteriorMean = meanOf(rhoSamples);
        s.lower2_5 = quantile(rhoSamples, 0.025);
        s.upper97_5 = quantile(rhoSamples, 0.975);
        s.nAccepted = nAccepted;
        return s;
    }

    private static double logLikelihood(List<LdPairLoader.Pair> pairs,
                                         double rho, double sigma) {
        double logL = 0.0;
        double invTwoSigma2 = 1.0 / (2.0 * sigma * sigma);
        for (LdPairLoader.Pair p : pairs) {
            double mu = HudsonRecombination.expectedR2(rho * p.distance);
            double diff = p.r2 - mu;
            logL -= diff * diff * invTwoSigma2;
        }
        return logL;
    }

    private static double meanOf(double[] xs) {
        double s = 0.0;
        for (double x : xs) s += x;
        return xs.length > 0 ? s / xs.length : Double.NaN;
    }

    private static double quantile(double[] sortedXs, double q) {
        if (sortedXs.length == 0) return Double.NaN;
        double idx = q * (sortedXs.length - 1);
        int lo = (int) Math.floor(idx);
        int hi = (int) Math.ceil(idx);
        if (lo == hi) return sortedXs[lo];
        double frac = idx - lo;
        return sortedXs[lo] * (1.0 - frac) + sortedXs[hi] * frac;
    }

    private static long maxVariantPosition(Transaction tx) {
        Result r = tx.execute(
                "MATCH (v:Variant) RETURN max(v.position) AS p");
        try {
            if (r.hasNext()) {
                Object p = r.next().get("p");
                if (p instanceof Number) return ((Number) p).longValue();
            }
        } finally {
            r.close();
        }
        return 0L;
    }

    private static String peekRunId(Transaction tx) {
        Result r = tx.execute(
                "MATCH ()-[m:MUTATED_ON]->() RETURN m.runId AS rid LIMIT 1");
        try {
            if (r.hasNext()) {
                Object rid = r.next().get("rid");
                if (rid != null) return (String) rid;
            }
        } finally {
            r.close();
        }
        return "";
    }

    private static long getLong(Map<String, Object> opts, String key, long def) {
        Object v = opts.get(key);
        return (v instanceof Number) ? ((Number) v).longValue() : def;
    }

    private static double getDouble(Map<String, Object> opts, String key,
                                    double def) {
        Object v = opts.get(key);
        return (v instanceof Number) ? ((Number) v).doubleValue() : def;
    }
}
