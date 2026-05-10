package org.graphpop.procedures.arg;

/**
 * Result row for {@code graphpop.arg.tmrca}.
 *
 * <p>For single-position queries: one row per pair with
 * {@code position} set and {@code mean_tmrca} = NaN.
 * For window-mean queries: {@code position} = -1 and
 * {@code mean_tmrca} is the span-weighted mean across marginal
 * trees overlapping the window.</p>
 */
public class TmrcaResult {

    public String sample_a;
    public String sample_b;
    public long position;
    public double tmrca;
    public double mean_tmrca;
    public long mrca_node_id;
    public String runId;

    public TmrcaResult() {}

    public TmrcaResult(String a, String b, long position,
                        double tmrca, double meanTmrca,
                        long mrcaNodeId, String runId) {
        this.sample_a = a;
        this.sample_b = b;
        this.position = position;
        this.tmrca = tmrca;
        this.mean_tmrca = meanTmrca;
        this.mrca_node_id = mrcaNodeId;
        this.runId = runId;
    }
}
