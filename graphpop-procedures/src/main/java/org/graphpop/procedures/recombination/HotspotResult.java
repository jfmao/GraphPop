package org.graphpop.procedures.recombination;

/**
 * Per-window row from {@code graphpop.recombination.hotspots}.
 */
public class HotspotResult {

    public long start;
    public long end;
    public double rho_per_bp;
    public double z_score;
    public double p_value;
    public double adj_p_value;
    public boolean is_hotspot;
    public long n_variant_pairs;
    public String method;
    public String runId;

    public HotspotResult() {}

    public HotspotResult(long start, long end, double rhoPerBp,
                         double zScore, double pValue, double adjPValue,
                         boolean isHotspot, long nVariantPairs,
                         String method, String runId) {
        this.start = start;
        this.end = end;
        this.rho_per_bp = rhoPerBp;
        this.z_score = zScore;
        this.p_value = pValue;
        this.adj_p_value = adjPValue;
        this.is_hotspot = isHotspot;
        this.n_variant_pairs = nVariantPairs;
        this.method = method;
        this.runId = runId;
    }
}
