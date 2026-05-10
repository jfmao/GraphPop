package org.graphpop.procedures.demography;

/**
 * Per-bin row from {@code graphpop.demography.ne_trajectory}.
 */
public class NeTrajectoryResult {

    public double time_lo;
    public double time_hi;
    public double n_coalescent_events;
    public double lineage_pair_time;
    public double rate;
    public double ne;
    public double ne_se;
    public String flag;
    public String runId;

    public NeTrajectoryResult() {}

    public NeTrajectoryResult(double tLo, double tHi,
                               double events, double pairTime,
                               double rate, double ne, double neSe,
                               String flag, String runId) {
        this.time_lo = tLo;
        this.time_hi = tHi;
        this.n_coalescent_events = events;
        this.lineage_pair_time = pairTime;
        this.rate = rate;
        this.ne = ne;
        this.ne_se = neSe;
        this.flag = flag;
        this.runId = runId;
    }
}
