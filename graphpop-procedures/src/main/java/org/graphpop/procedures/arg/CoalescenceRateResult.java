package org.graphpop.procedures.arg;

/**
 * Per-bin row from {@code graphpop.arg.coalescence_rate}.
 */
public class CoalescenceRateResult {

    public double time_lo;
    public double time_hi;
    public double n_coalescent_events;
    public double lineage_pair_time;
    public double rate;
    public String runId;

    public CoalescenceRateResult() {}

    public CoalescenceRateResult(double tLo, double tHi,
                                  double events, double pairTime,
                                  double rate, String runId) {
        this.time_lo = tLo;
        this.time_hi = tHi;
        this.n_coalescent_events = events;
        this.lineage_pair_time = pairTime;
        this.rate = rate;
        this.runId = runId;
    }
}
