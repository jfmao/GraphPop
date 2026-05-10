package org.graphpop.procedures.recombination;

/**
 * Discrete-state HMM over per-window log(ρ) with Gaussian
 * emissions. Pyrho-style smoothing (Spence &amp; Song 2019) at v1.
 *
 * <p>States are bins of {@code log10(ρ)} between
 * {@code stateLogLo} and {@code stateLogHi}. Transition matrix
 * is row-stochastic banded: state {@code i} transitions to itself
 * with prob {@code 1 − switchRate}, to its immediate neighbours
 * with prob {@code switchRate / 2} each (with reflecting
 * boundaries). Emission is Gaussian on {@code log10(ρ̂_window)}
 * given true state with std {@code emissionSd}.</p>
 *
 * <p>Returns the per-window posterior mean of {@code log10(ρ)}
 * (forward-backward), plus the Viterbi state.</p>
 */
public final class HmmSmoothComputer {

    private HmmSmoothComputer() {}

    public static final class Result {
        public final double[] smoothedLog10Rho;  // length = n_windows
        public final int[] viterbiStates;         // length = n_windows

        Result(double[] smoothed, int[] viterbi) {
            this.smoothedLog10Rho = smoothed;
            this.viterbiStates = viterbi;
        }
    }

    /**
     * Run forward-backward + Viterbi over observed log10(ρ̂) values.
     *
     * @param observedLog10Rho per-window observed log10(ρ); NaN for
     *                          missing-data windows (treated as
     *                          uninformative).
     * @param nStates           number of discrete states.
     * @param stateLogLo        lower bound of log10(ρ) state grid.
     * @param stateLogHi        upper bound of log10(ρ) state grid.
     * @param emissionSd        Gaussian emission std on log10 scale.
     * @param switchRate        per-window probability of switching
     *                           to an adjacent state (∈ [0, 1]).
     */
    public static Result smooth(double[] observedLog10Rho, int nStates,
                                  double stateLogLo, double stateLogHi,
                                  double emissionSd, double switchRate) {
        if (nStates < 2) {
            throw new IllegalArgumentException("nStates must be ≥ 2");
        }
        if (stateLogHi <= stateLogLo) {
            throw new IllegalArgumentException(
                    "state_log_hi must be > state_log_lo");
        }
        if (emissionSd <= 0.0) {
            throw new IllegalArgumentException("emission_sd must be > 0");
        }
        if (switchRate < 0.0 || switchRate > 1.0) {
            throw new IllegalArgumentException(
                    "switch_rate must be in [0, 1]");
        }

        int T = observedLog10Rho.length;
        double[] stateLog = new double[nStates];
        for (int s = 0; s < nStates; s++) {
            stateLog[s] = stateLogLo
                    + (stateLogHi - stateLogLo) * s / (nStates - 1);
        }

        // Pre-compute emission log-probabilities (Gaussian) per (t, s).
        double[][] logEmit = new double[T][nStates];
        double inv2Sigma2 = 1.0 / (2.0 * emissionSd * emissionSd);
        double logNorm = -Math.log(emissionSd * Math.sqrt(2.0 * Math.PI));
        for (int t = 0; t < T; t++) {
            double obs = observedLog10Rho[t];
            if (Double.isNaN(obs)) {
                // Uninformative: uniform across states.
                double logUniform = -Math.log(nStates);
                for (int s = 0; s < nStates; s++) logEmit[t][s] = logUniform;
            } else {
                for (int s = 0; s < nStates; s++) {
                    double diff = obs - stateLog[s];
                    logEmit[t][s] = logNorm - diff * diff * inv2Sigma2;
                }
            }
        }

        // Build banded transition log-probability matrix.
        double[][] logTrans = new double[nStates][nStates];
        double stayLog = Math.log(Math.max(1.0 - switchRate, 1e-300));
        double moveLog = Math.log(Math.max(switchRate / 2.0, 1e-300));
        for (int i = 0; i < nStates; i++) {
            for (int j = 0; j < nStates; j++) {
                logTrans[i][j] = Double.NEGATIVE_INFINITY;
            }
            logTrans[i][i] = stayLog;
            if (i > 0) logTrans[i][i - 1] = moveLog;
            if (i < nStates - 1) logTrans[i][i + 1] = moveLog;
            // Reflecting boundary: at i=0 and i=nStates-1, half of the
            // mass that "would have" gone outside reflects back to
            // self. Equivalent: increase the stay-prob at the edges.
            if (i == 0 || i == nStates - 1) {
                logTrans[i][i] = Math.log(1.0 - switchRate / 2.0);
            }
        }

        // Forward.
        double[][] logF = new double[T][nStates];
        double logUniformPrior = -Math.log(nStates);
        for (int s = 0; s < nStates; s++) {
            logF[0][s] = logUniformPrior + logEmit[0][s];
        }
        for (int t = 1; t < T; t++) {
            for (int s = 0; s < nStates; s++) {
                double m = Double.NEGATIVE_INFINITY;
                double[] xs = new double[nStates];
                for (int p = 0; p < nStates; p++) {
                    xs[p] = logF[t - 1][p] + logTrans[p][s];
                    if (xs[p] > m) m = xs[p];
                }
                double sumExp = 0.0;
                for (int p = 0; p < nStates; p++) {
                    if (xs[p] > Double.NEGATIVE_INFINITY) {
                        sumExp += Math.exp(xs[p] - m);
                    }
                }
                logF[t][s] = m + Math.log(sumExp) + logEmit[t][s];
            }
        }

        // Backward.
        double[][] logB = new double[T][nStates];
        // logB[T-1] = 0 (log 1).
        for (int t = T - 2; t >= 0; t--) {
            for (int s = 0; s < nStates; s++) {
                double m = Double.NEGATIVE_INFINITY;
                double[] xs = new double[nStates];
                for (int n = 0; n < nStates; n++) {
                    xs[n] = logTrans[s][n] + logEmit[t + 1][n]
                            + logB[t + 1][n];
                    if (xs[n] > m) m = xs[n];
                }
                double sumExp = 0.0;
                for (int n = 0; n < nStates; n++) {
                    if (xs[n] > Double.NEGATIVE_INFINITY) {
                        sumExp += Math.exp(xs[n] - m);
                    }
                }
                logB[t][s] = m + Math.log(sumExp);
            }
        }

        // Posterior + smoothed log10(ρ) = E[state | data].
        double[] smoothed = new double[T];
        int[] viterbi = new int[T];
        for (int t = 0; t < T; t++) {
            double m = Double.NEGATIVE_INFINITY;
            double[] logPost = new double[nStates];
            for (int s = 0; s < nStates; s++) {
                logPost[s] = logF[t][s] + logB[t][s];
                if (logPost[s] > m) m = logPost[s];
            }
            double sumExp = 0.0;
            for (int s = 0; s < nStates; s++) {
                sumExp += Math.exp(logPost[s] - m);
            }
            double logNormPost = m + Math.log(sumExp);
            double meanLog = 0.0;
            int best = 0;
            double bestLp = Double.NEGATIVE_INFINITY;
            for (int s = 0; s < nStates; s++) {
                double p = Math.exp(logPost[s] - logNormPost);
                meanLog += p * stateLog[s];
                if (logPost[s] > bestLp) {
                    bestLp = logPost[s];
                    best = s;
                }
            }
            smoothed[t] = meanLog;
            viterbi[t] = best;
        }
        return new Result(smoothed, viterbi);
    }
}
