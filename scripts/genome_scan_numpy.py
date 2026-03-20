#!/usr/bin/env python3
"""
Pure-Python/numpy sliding-window genome scan.

Replaces graphpop.genome_scan (Neo4j WRITE procedure) for the gscan subcommand.
Reads from pre-exported .npz chromosome cache instead of traversing the graph,
enabling parallel execution across chromosomes with no DB lock contention.

Implements exactly the same statistics as GenomeScanProcedure.java (M1.4):
  - pi, theta_w, Tajima's D  (Tajima 1989)
  - Fay & Wu's H, normalised H  (Zeng et al. 2006)
  - het_exp, het_obs, Fis
  - Fst Hudson (numerator/denominator)
  - Fst Weir & Cockerham (a/b/c components)
  - Dxy, Da
  - PBS  (3-population, optional)

Output: list of window dicts matching the existing JSON checkpoint schema,
compatible with all downstream report/summarisation code.
"""

import json
import os

import numpy as np

# ── Default scan parameters ────────────────────────────────────────────────────

WINDOW_SIZE = 100_000
WINDOW_STEP =  50_000
CACHE_DIR   = "data/variants_cache"

# ── Cache I/O ──────────────────────────────────────────────────────────────────

_pop_meta_cache = {}   # module-level: loaded once per process


def load_pop_meta(cache_dir=CACHE_DIR):
    if cache_dir not in _pop_meta_cache:
        with open(os.path.join(cache_dir, "pop_meta.json")) as f:
            _pop_meta_cache[cache_dir] = json.load(f)
    return _pop_meta_cache[cache_dir]


def _pop_index(pop_ids, pop):
    try:
        return list(pop_ids).index(pop)
    except ValueError:
        raise ValueError(f"Population '{pop}' not found in pop_ids: {pop_ids}")


# ── Harmonic numbers ───────────────────────────────────────────────────────────

def _harmonic(n):
    """H(n) = sum(1/i, i=1..n). Exact, matching TajimaD.java."""
    return sum(1.0 / i for i in range(1, n + 1))


def _harmonic2(n):
    """H2(n) = sum(1/i^2, i=1..n)."""
    return sum(1.0 / (i * i) for i in range(1, n + 1))


# ── Core scan function ─────────────────────────────────────────────────────────

def scan(chrom, pop1, pop2=None, pop3=None, consequence=None,
         window_size=WINDOW_SIZE, step_size=WINDOW_STEP,
         cache_dir=CACHE_DIR):
    """
    Sliding-window genome scan on pre-exported numpy cache.

    Parameters
    ----------
    chrom       : str   e.g. 'chr1'
    pop1        : str   focal population
    pop2        : str   comparative population (Fst, Dxy)
    pop3        : str   outgroup population (PBS)
    consequence : str   'missense_variant', 'synonymous_variant', or None
    window_size : int   window size in bp (default 100 000)
    step_size   : int   step size in bp   (default  50 000)
    cache_dir   : str   directory containing .npz files and pop_meta.json

    Returns
    -------
    list of dicts — window results matching the GenomeScanProcedure.WindowResult schema
    """
    meta    = load_pop_meta(cache_dir)
    pop_ids = meta["pop_ids"]

    idx1 = _pop_index(pop_ids, pop1)
    idx2 = _pop_index(pop_ids, pop2) if pop2 else None
    idx3 = _pop_index(pop_ids, pop3) if pop3 else None

    # ── Load cache arrays one at a time (minimise peak RAM) ───────────────────
    npz_path = os.path.join(cache_dir, f"{chrom}.npz")
    if not os.path.exists(npz_path):
        raise FileNotFoundError(
            f"Cache file not found: {npz_path}. "
            "Run scripts/export_variant_cache.py first."
        )

    with np.load(npz_path, allow_pickle=False) as data:
        pos       = data["pos"].copy()                          # int64 (M,)
        anc_code  = data["ancestral_allele"].copy()             # int8  (M,)

        # Apply consequence mask first — reduces working set for filtered scans
        if consequence == "missense_variant":
            csq_mask = data["has_missense"].copy()
        elif consequence == "synonymous_variant":
            csq_mask = data["has_synonymous"].copy()
        else:
            csq_mask = np.ones(len(pos), dtype=bool)

        # Extract only the population columns we need, then free the big arrays
        ac1  = data["ac"][:, idx1].astype(np.float64);   # copy releases (M,26)
        an1  = data["an"][:, idx1].astype(np.float64)
        af1  = data["af"][:, idx1].astype(np.float64)    # float32→float64
        het1 = data["het_count"][:, idx1].astype(np.float64)

        ac2 = an2 = af2 = het2 = None
        if idx2 is not None:
            ac2  = data["ac"][:, idx2].astype(np.float64)
            an2  = data["an"][:, idx2].astype(np.float64)
            af2  = data["af"][:, idx2].astype(np.float64)
            het2 = data["het_count"][:, idx2].astype(np.float64)

        ac3 = an3 = af3 = het3 = None
        if idx3 is not None:
            ac3  = data["ac"][:, idx3].astype(np.float64)
            an3  = data["an"][:, idx3].astype(np.float64)
            af3  = data["af"][:, idx3].astype(np.float64)
            het3 = data["het_count"][:, idx3].astype(np.float64)
    # npz file closed and large (M,26) arrays freed here

    M = len(pos)

    # ── Variant validity mask ─────────────────────────────────────────────────
    valid = csq_mask & (an1 >= 2)
    del csq_mask

    valid_idx = np.where(valid)[0]
    if len(valid_idx) == 0:
        return []

    # Ploidy-aware n: use AN from the first valid variant (matches Java behaviour)
    first_an = int(an1[valid_idx[0]])
    n = first_an                                 # haploid AN = 2 × diploid n
    a_n  = _harmonic(n - 1)   if n > 1 else 0.0
    a_n2 = _harmonic2(n - 1)  if n > 1 else 0.0

    # ── Phase 1: Pre-compute per-variant statistics (vectorised, M1.4) ────────

    # Nucleotide diversity π per site:  2·p·(1−p)·n/(n−1)
    with np.errstate(invalid="ignore", divide="ignore"):
        pi_arr = np.where(valid & (an1 > 1),
                          2.0 * af1 * (1.0 - af1) * an1 / (an1 - 1.0), 0.0)

    # Expected / observed heterozygosity
    he_arr = np.where(valid, 2.0 * af1 * (1.0 - af1), 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        ndip1  = np.where(valid, an1 / 2.0, 1.0)
        ho_arr = np.where(valid & (ndip1 > 0), het1 / ndip1, 0.0)

    # Segregating sites
    seg_arr = valid & (ac1 > 0) & (ac1 < an1)

    # Fay & Wu's H per-site theta_H:  2·i²/(n·(n−1)),  i = derived allele count
    derived = np.where(anc_code == 1, ac1,              # REF is ancestral → derived = ac
               np.where(anc_code == 2, an1 - ac1, -1.0))  # ALT is ancestral → derived = an−ac
    pol_valid = valid & (derived >= 0) & (derived > 0) & (derived < an1)
    with np.errstate(invalid="ignore", divide="ignore"):
        theta_h_arr  = np.where(pol_valid,
                                2.0 * derived**2 / (an1 * (an1 - 1.0)), 0.0)
        pi_polar_arr = np.where(pol_valid, pi_arr, 0.0)

    # ── Comparative statistics (pop2) ─────────────────────────────────────────
    fst_num_arr = fst_den_arr = dxy_arr = pi_w2_arr = None
    wc_a_arr = wc_b_arr = wc_c_arr = None
    # PBS: pop1-vs-pop3 and pop2-vs-pop3 W&C components
    wc_a13 = wc_b13 = wc_c13 = None
    wc_a23 = wc_b23 = wc_c23 = None
    valid_comp = None

    if idx2 is not None:
        valid_comp = valid & (an2 >= 2)

        # Hudson Fst components:
        #   hw1 = 2p1(1−p1)n1/(n1−1),  hw2 = 2p2(1−p2)n2/(n2−1)
        #   hb  = (p1−p2)² − hw1/(2n1) − hw2/(2n2)
        #   num = hb,  den = hb + (hw1+hw2)/2
        with np.errstate(invalid="ignore", divide="ignore"):
            hw1 = 2.0 * af1 * (1.0 - af1) * an1 / (an1 - 1.0)
            hw2 = 2.0 * af2 * (1.0 - af2) * an2 / (an2 - 1.0)
            hb  = (af1 - af2)**2 - hw1 / (2.0 * an1) - hw2 / (2.0 * an2)
            fst_num_arr = np.where(valid_comp, hb, 0.0)
            fst_den_arr = np.where(valid_comp, hb + (hw1 + hw2) / 2.0, 0.0)

        # Dxy per site: p1(1−p2) + p2(1−p1)
        dxy_arr  = np.where(valid_comp, af1 * (1.0 - af2) + af2 * (1.0 - af1), 0.0)
        pi_w2_arr = np.where(valid_comp & (an2 > 1),
                             2.0 * af2 * (1.0 - af2) * an2 / (an2 - 1.0), 0.0)

        # Weir & Cockerham Fst components (r=2, W&C 1984 Eqs. 2–4)
        wc_a_arr, wc_b_arr, wc_c_arr = _wc_components(
            ac1, an1, het1, ac2, an2, het2, valid_comp)

        if idx3 is not None:
            valid_3 = valid & (an2 >= 2) & (an3 >= 2)
            wc_a13, wc_b13, wc_c13 = _wc_components(
                ac1, an1, het1, ac3, an3, het3, valid_3)
            wc_a23, wc_b23, wc_c23 = _wc_components(
                ac2, an2, het2, ac3, an3, het3, valid_3)

    # ── Phase 2: Prefix sums for O(1) window queries ──────────────────────────

    def _cumsum(arr):
        return np.concatenate(([0.0], np.cumsum(arr)))

    cum_pi       = _cumsum(pi_arr)
    cum_he       = _cumsum(he_arr)
    cum_ho       = _cumsum(ho_arr)
    cum_seg      = _cumsum(seg_arr.astype(np.float64))
    cum_theta_h  = _cumsum(theta_h_arr)
    cum_pi_polar = _cumsum(pi_polar_arr)
    cum_valid    = _cumsum(valid.astype(np.float64))
    cum_pol      = _cumsum(pol_valid.astype(np.float64))

    cum_fst_num = cum_fst_den = cum_dxy = cum_pi_w1 = cum_pi_w2 = None
    cum_wca = cum_wcb = cum_wcc = None
    cum_wca13 = cum_wcb13 = cum_wcc13 = None
    cum_wca23 = cum_wcb23 = cum_wcc23 = None
    cum_vcomp = None

    if idx2 is not None:
        cum_fst_num = _cumsum(fst_num_arr)
        cum_fst_den = _cumsum(fst_den_arr)
        cum_dxy     = _cumsum(dxy_arr)
        cum_pi_w1   = _cumsum(np.where(valid_comp, pi_arr, 0.0))
        cum_pi_w2   = _cumsum(pi_w2_arr)
        cum_wca     = _cumsum(wc_a_arr)
        cum_wcb     = _cumsum(wc_b_arr)
        cum_wcc     = _cumsum(wc_c_arr)
        cum_vcomp   = _cumsum(valid_comp.astype(np.float64))
        if idx3 is not None:
            cum_wca13 = _cumsum(wc_a13)
            cum_wcb13 = _cumsum(wc_b13)
            cum_wcc13 = _cumsum(wc_c13)
            cum_wca23 = _cumsum(wc_a23)
            cum_wcb23 = _cumsum(wc_b23)
            cum_wcc23 = _cumsum(wc_c23)

    # ── Phase 3: Window scan ───────────────────────────────────────────────────

    min_pos = int(pos[valid_idx[0]])
    max_pos = int(pos[valid_idx[-1]])

    w_starts = np.arange(min_pos, max_pos + 1, step_size, dtype=np.int64)
    w_ends   = w_starts + window_size - 1

    # Boundary indices into pos array via searchsorted
    lo_idx = np.searchsorted(pos, w_starts, side="left")
    hi_idx = np.searchsorted(pos, w_ends,   side="right")

    run_id = f"scan_{pop1}_np"

    results = []
    for w_start, w_end, lo, hi in zip(
            w_starts.tolist(), w_ends.tolist(),
            lo_idx.tolist(),   hi_idx.tolist()):

        n_variants = int(cum_valid[hi] - cum_valid[lo])
        if n_variants == 0:
            continue

        L = float(n_variants)
        pi_sum   = float(cum_pi[hi]  - cum_pi[lo])
        he_sum   = float(cum_he[hi]  - cum_he[lo])
        ho_sum   = float(cum_ho[hi]  - cum_ho[lo])
        n_seg    = int(  cum_seg[hi] - cum_seg[lo])
        n_pol    = int(  cum_pol[hi] - cum_pol[lo])

        # Tajima's D  (TajimaD.java formula)
        if n_seg > 0 and a_n > 0 and n >= 4:
            theta_w = n_seg / a_n
            d = pi_sum - theta_w
            b1 = (n + 1) / (3.0 * (n - 1))
            b2 = 2.0 * (n**2 + n + 3) / (9.0 * n * (n - 1))
            c1 = b1 - 1.0 / a_n
            c2 = b2 - (n + 2) / (a_n * n) + a_n2 / (a_n**2)
            e1 = c1 / a_n
            e2 = c2 / (a_n**2 + a_n2)
            var_d = e1 * n_seg + e2 * n_seg * (n_seg - 1)
            tajima_d      = d / var_d**0.5 if var_d > 0 else 0.0
            theta_w_site  = theta_w / L
        else:
            tajima_d     = 0.0
            theta_w_site = 0.0

        # Fay & Wu's H  (FayWuH.java formulas)
        fwh = fwh_norm = 0.0
        if n_pol > 0:
            th_sum  = float(cum_theta_h[hi]  - cum_theta_h[lo])
            pip_sum = float(cum_pi_polar[hi] - cum_pi_polar[lo])
            fwh = (pip_sum - th_sum) / n_pol
            if n >= 4 and a_n > 0:
                theta_s = n_pol / a_n
                term1   = theta_s * (n - 2) / (6.0 * (n - 1))
                term2   = theta_s**2 * (
                    18 * n**2 * (3 * n + 2) * a_n2
                    - (88 * n**3 + 9 * n**2 - 13 * n + 6)
                ) / (9.0 * n * (n - 1)**2)
                var_h   = term1 + term2
                fwh_norm = fwh / var_h**0.5 if var_h > 0 else 0.0

        het_exp = he_sum / L
        het_obs = ho_sum / L
        fis     = (1.0 - het_obs / het_exp) if het_exp > 0 else 0.0

        wr = {
            "window_id":      f"{chrom}:{w_start}-{w_end}:{pop1}:{run_id}",
            "chr":            chrom,
            "start":          w_start,
            "end":            w_end,
            "population":     pop1,
            "run_id":         run_id,
            "n_variants":     n_variants,
            "n_segregating":  n_seg,
            "pi":             pi_sum / L,
            "theta_w":        theta_w_site,
            "tajima_d":       tajima_d,
            "fay_wu_h":       fwh      if n_pol > 0 else 0.0,
            "fay_wu_h_norm":  fwh_norm if n_pol > 0 else 0.0,
            "n_polarized":    n_pol,
            "het_exp":        het_exp,
            "het_obs":        het_obs,
            "fis":            fis,
        }

        if idx2 is not None:
            n_comp = int(cum_vcomp[hi] - cum_vcomp[lo])
            if n_comp > 0:
                fst_num = float(cum_fst_num[hi] - cum_fst_num[lo])
                fst_den = float(cum_fst_den[hi] - cum_fst_den[lo])
                dxy_sum = float(cum_dxy[hi]     - cum_dxy[lo])
                piw1    = float(cum_pi_w1[hi]   - cum_pi_w1[lo])
                piw2    = float(cum_pi_w2[hi]   - cum_pi_w2[lo])
                wca     = float(cum_wca[hi]     - cum_wca[lo])
                wcb     = float(cum_wcb[hi]     - cum_wcb[lo])
                wcc     = float(cum_wcc[hi]     - cum_wcc[lo])

                wr["fst"]    = fst_num / fst_den if fst_den > 0 else 0.0
                wc_denom     = wca + wcb + wcc
                wr["fst_wc"] = wca / wc_denom if wc_denom > 0 else 0.0
                wr["dxy"]    = dxy_sum / L
                wr["da"]     = wr["dxy"] - (piw1 / L + piw2 / L) / 2.0
                wr["pop2"]   = pop2

                if idx3 is not None:
                    wca13 = float(cum_wca13[hi] - cum_wca13[lo])
                    wcb13 = float(cum_wcb13[hi] - cum_wcb13[lo])
                    wcc13 = float(cum_wcc13[hi] - cum_wcc13[lo])
                    wca23 = float(cum_wca23[hi] - cum_wca23[lo])
                    wcb23 = float(cum_wcb23[hi] - cum_wcb23[lo])
                    wcc23 = float(cum_wcc23[hi] - cum_wcc23[lo])

                    d12   = wca + wcb + wcc
                    d13   = wca13 + wcb13 + wcc13
                    d23   = wca23 + wcb23 + wcc23
                    fst12 = wca   / d12 if d12 > 0 else 0.0
                    fst13 = wca13 / d13 if d13 > 0 else 0.0
                    fst23 = wca23 / d23 if d23 > 0 else 0.0

                    t12   = -np.log(max(1e-10, 1.0 - fst12))
                    t13   = -np.log(max(1e-10, 1.0 - fst13))
                    t23   = -np.log(max(1e-10, 1.0 - fst23))
                    wr["pbs"]  = (t12 + t13 - t23) / 2.0
                    wr["pop3"] = pop3

        results.append(wr)

    return results


# ── W&C helper ─────────────────────────────────────────────────────────────────

def _wc_components(ac1, an1, het1, ac2, an2, het2, mask):
    """
    Weir & Cockerham (1984) Fst components a, b, c for r=2 populations.
    Vectorised over M variants; masked to `mask`.
    Returns three float64 arrays of shape (M,).
    """
    with np.errstate(invalid="ignore", divide="ignore"):
        n1 = an1 / 2.0         # diploid sample counts
        n2 = an2 / 2.0
        n_total = n1 + n2
        n_bar   = n_total / 2.0
        n_c     = n_total - (n1**2 + n2**2) / n_total  # correction factor

        p1 = np.where(an1 > 0, ac1 / an1, 0.0)
        p2 = np.where(an2 > 0, ac2 / an2, 0.0)
        p_bar = np.where(n_total > 0, (n1 * p1 + n2 * p2) / n_total, 0.0)

        # s² = Σ nᵢ·(pᵢ−p̄)² / ((r−1)·n̄),  r=2
        s2 = np.where(n_bar > 0,
                      (n1 * (p1 - p_bar)**2 + n2 * (p2 - p_bar)**2) / n_bar,
                      0.0)

        h_bar = np.where(n1 * n2 > 0,
                         (het1 / n1 + het2 / n2) / 2.0, 0.0)

        # W&C Eq. 2 (a component)
        a = np.where(
            mask & (n_bar > 1) & (n_c > 0),
            n_bar / n_c * (
                s2 - 1.0 / (n_bar - 1.0) * (
                    p_bar * (1.0 - p_bar) - s2 / 2.0 - h_bar / 4.0
                )
            ),
            0.0,
        )
        # W&C Eq. 3 (b component)
        b = np.where(
            mask & (n_bar > 1),
            n_bar / (n_bar - 1.0) * (
                p_bar * (1.0 - p_bar)
                - s2 / 2.0
                - h_bar * (2.0 * n_bar - 1.0) / (4.0 * n_bar)
            ),
            0.0,
        )
        # W&C Eq. 4 (c component)
        c = np.where(mask, h_bar / 2.0, 0.0)

    return a, b, c


# ── Smoke-test helper ──────────────────────────────────────────────────────────

def smoke_test(chrom="chr22", pop1="YRI", pop2="CEU",
               cache_dir=CACHE_DIR):
    """
    Quick sanity check: run one scan and print summary statistics.
    Compare manually against the Neo4j result for the same call.
    """
    import time
    t0 = time.time()
    windows = scan(chrom, pop1, pop2, cache_dir=cache_dir)
    elapsed = time.time() - t0

    if not windows:
        print(f"No windows returned for {chrom} {pop1} vs {pop2}")
        return

    pis      = [w["pi"]      for w in windows]
    fst_wcs  = [w.get("fst_wc", 0) for w in windows]
    tds      = [w["tajima_d"] for w in windows]

    print(f"Scan: {chrom}  {pop1} vs {pop2}  "
          f"({len(windows)} windows, {elapsed:.2f}s)")
    print(f"  pi     mean={sum(pis)/len(pis):.6f}  "
          f"max={max(pis):.6f}")
    print(f"  fst_wc mean={sum(fst_wcs)/len(fst_wcs):.6f}  "
          f"max={max(fst_wcs):.6f}")
    print(f"  tajima_d mean={sum(tds)/len(tds):.4f}")


if __name__ == "__main__":
    import sys
    chrom  = sys.argv[1] if len(sys.argv) > 1 else "chr22"
    pop1   = sys.argv[2] if len(sys.argv) > 2 else "YRI"
    pop2   = sys.argv[3] if len(sys.argv) > 3 else "CEU"
    smoke_test(chrom, pop1, pop2)
