#!/usr/bin/env python3
"""Quick validation of ld_graphpop_numpy vs PLINK2 for the 'large' region."""
import sys, os, time, subprocess, tempfile, json
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np

# ---- configuration ----
CHR   = "chr22"
START = 16_000_000
END   = 17_000_000
POP   = "CEU"
ROOT  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HAP_CACHE_DIR = os.path.join(ROOT, "data/hap_cache")
VCF = os.path.join(ROOT,
    "data/raw/1000g/CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz")

# Minimal load_hap_cache (mirrors benchmark-vs-tools.py)
def load_hap_cache(pop, start, end):
    meta = json.load(open(os.path.join(HAP_CACHE_DIR, "sample_meta.json")))
    pop_ids   = meta["pop_ids"]
    n_samples = meta["n_samples"]
    sample_idx = np.array([i for i, p in enumerate(pop_ids) if p == pop], dtype=np.int32)
    d          = np.load(os.path.join(HAP_CACHE_DIR, f"{CHR}.npz"), mmap_mode="r")
    pos        = np.array(d["pos"])
    hap_packed = np.array(d["hap"])
    mask       = (pos >= start) & (pos <= end)
    pos        = pos[mask]
    hap_packed = hap_packed[mask]
    hap_full   = np.unpackbits(hap_packed, axis=1, bitorder="big")[:, :2 * n_samples]
    hap_idx    = np.empty(2 * len(sample_idx), dtype=np.int32)
    hap_idx[0::2] = 2 * sample_idx
    hap_idx[1::2] = 2 * sample_idx + 1
    return pos, hap_full[:, hap_idx].astype(np.int8)


def ld_numpy(max_dist_bp, r2_threshold, min_maf):
    """Numpy pairwise r² from hap_cache (mirrors ld_graphpop_numpy)."""
    t0 = time.perf_counter()
    pos, hap = load_hap_cache(POP, START, END)

    hap_f = hap.astype(np.float32)
    n_hap = float(hap_f.shape[1])
    freq  = hap_f.mean(axis=1)
    maf   = np.minimum(freq, 1.0 - freq)
    mask  = maf >= min_maf
    hap_f = hap_f[mask]; pos2 = pos[mask]; freq2 = freq[mask]
    var   = (freq2 * (1.0 - freq2)).clip(min=1e-9).astype(np.float32)
    M     = len(pos2)
    print(f"  Numpy: {M} variants after MAF>={min_maf}", flush=True)

    FULL_THRESH = 5000
    n_pairs = 0
    if M <= FULL_THRESH:
        cov    = hap_f @ hap_f.T / n_hap - np.outer(freq2, freq2)
        r2     = (cov ** 2) / np.outer(var, var)
        i_idx, j_idx = np.triu_indices(M, k=1)
        in_win = (pos2[j_idx] - pos2[i_idx]) <= max_dist_bp
        n_pairs = int((r2[i_idx[in_win], j_idx[in_win]] >= r2_threshold).sum())
    else:
        BLOCK = 512
        j_max = np.searchsorted(pos2, pos2 + max_dist_bp, side="right")
        for i0 in range(0, M, BLOCK):
            i1    = min(i0 + BLOCK, M)
            j_end = int(j_max[i1 - 1])
            hi    = hap_f[i0:i1]
            hj    = hap_f[i0:j_end]
            cov   = hi @ hj.T / n_hap - np.outer(freq2[i0:i1], freq2[i0:j_end])
            denom = np.outer(var[i0:i1], var[i0:j_end])
            r2    = cov ** 2 / denom
            gi    = np.arange(i0, i1)[:, None]
            gj    = np.arange(i0, j_end)[None, :]
            valid = (gj > gi) & (
                (pos2[i0:j_end][None, :] - pos2[i0:i1, None]) <= max_dist_bp)
            n_pairs += int((r2[valid] >= r2_threshold).sum())

    t = time.perf_counter() - t0
    return {"n_pairs": n_pairs, "n_variants": M, "time_s": t}


def ld_plink2(max_dist_kb, r2_threshold, min_maf, use_pgen=True):
    """PLINK2 --r2-unphased for comparison."""
    panel = os.path.join(ROOT, "data/raw/1000g/integrated_call_samples_v3.20130502.ALL.panel")
    pgen  = os.path.join(ROOT, "data/raw/1000g/chr22_ceu")

    with tempfile.TemporaryDirectory() as tmpdir:
        # Write keep file
        keep_file = os.path.join(tmpdir, "keep.txt")
        with open(keep_file, "w") as f:
            for line in open(panel):
                parts = line.split()
                if len(parts) >= 2 and parts[1] == POP:
                    f.write(f"0\t{parts[0]}\n")

        outpfx = os.path.join(tmpdir, "plink_ld")
        if use_pgen and os.path.exists(pgen + ".pgen"):
            input_args = ["--pfile", pgen, "vzs"]
        else:
            input_args = ["--vcf", VCF]

        cmd = ["plink2", *input_args,
               "--chr", CHR,
               "--from-bp", str(START), "--to-bp", str(END),
               "--keep", keep_file,
               "--r2-unphased",
               "--ld-window-kb", str(max_dist_kb),
               "--ld-window-r2", str(r2_threshold),
               "--maf", str(min_maf),
               "--allow-extra-chr",
               "--out", outpfx]

        t0 = time.perf_counter()
        proc = subprocess.run(cmd, capture_output=True)
        t = time.perf_counter() - t0

        if proc.returncode != 0:
            return {"error": proc.stderr.decode()[-200:], "time_s": t}

        # Count pairs
        n_pairs = 0
        for ext in [".vcor", ".ld"]:
            ld_file = outpfx + ext
            if os.path.exists(ld_file):
                with open(ld_file) as f:
                    n_pairs = sum(1 for _ in f) - 1  # skip header
                break

        return {"n_pairs": n_pairs, "time_s": t}


# ============================================================
print(f"LD benchmark: {CHR}:{START:,}-{END:,}  pop={POP}", flush=True)
print("=" * 60, flush=True)

# D1: 100kb window, r²≥0.1, MAF≥0.05
print("\nD1: 100kb window, r²≥0.1, MAF≥0.05", flush=True)
print("  PLINK2 (.pgen)...", end="", flush=True)
p1 = ld_plink2(max_dist_kb=100, r2_threshold=0.1, min_maf=0.05, use_pgen=True)
print(f" n_pairs={p1.get('n_pairs', 'ERR'):,}  t={p1['time_s']:.3f}s", flush=True)

print("  numpy...", end="", flush=True)
n1 = ld_numpy(max_dist_bp=100_000, r2_threshold=0.1, min_maf=0.05)
print(f" n_pairs={n1['n_pairs']:,}  t={n1['time_s']:.3f}s", flush=True)

# Compare
if p1.get("n_pairs") and n1["n_pairs"]:
    delta = abs(n1["n_pairs"] - p1["n_pairs"]) / max(p1["n_pairs"], 1) * 100
    speedup = p1["time_s"] / n1["time_s"] if n1["time_s"] > 0 else float("inf")
    print(f"  -> Pair count diff: {delta:.1f}%  "
          f"Speedup vs PLINK2: {speedup:.1f}x", flush=True)

# D2: 500kb window, r²≥0, MAF≥0.05
print("\nD2: 500kb window, r²≥0 (all pairs), MAF≥0.05", flush=True)
print("  PLINK2 (.pgen)...", end="", flush=True)
p2 = ld_plink2(max_dist_kb=500, r2_threshold=0.0, min_maf=0.05, use_pgen=True)
print(f" n_pairs={p2.get('n_pairs', 'ERR'):,}  t={p2['time_s']:.3f}s", flush=True)

print("  numpy...", end="", flush=True)
n2 = ld_numpy(max_dist_bp=500_000, r2_threshold=0.0, min_maf=0.05)
print(f" n_pairs={n2['n_pairs']:,}  t={n2['time_s']:.3f}s", flush=True)

if p2.get("n_pairs") and n2["n_pairs"]:
    delta2 = abs(n2["n_pairs"] - p2["n_pairs"]) / max(p2["n_pairs"], 1) * 100
    speedup2 = p2["time_s"] / n2["time_s"] if n2["time_s"] > 0 else float("inf")
    print(f"  -> Pair count diff: {delta2:.1f}%  "
          f"Speedup vs PLINK2: {speedup2:.1f}x", flush=True)

print("\nDone.", flush=True)
