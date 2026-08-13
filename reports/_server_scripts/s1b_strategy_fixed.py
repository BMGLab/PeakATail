#!/usr/bin/env python
"""Count-controlled re-analysis of the peak-calling strategy comparison.

Supersedes an earlier ``s1_strategy.py`` (removed), which parsed the wrong column
out of ``bedtools closest`` output: the tool emits 18 columns for this atlas and
the distance is the last one, but the reader asked pandas for ``usecols=[-1]``,
which is not a valid column selector, so every distance came back NaN. The
distance is now extracted with awk on ``$NF`` rather than a hardcoded index.

It is also restructured to be roughly 20x cheaper.  The previous version ran
``bedtools closest`` against the 1.2 GB PolyASite v3 BED once per rarefaction
replicate -- 100 full passes over the atlas.  Subsampling the *calls* and then
measuring distance is statistically identical to measuring distance once and
then subsampling the resulting distance vector, so this version runs
``bedtools closest`` exactly once per strategy and rarefies in memory.

Why rarefy at all: precision against an 18.4M-entry reference saturates near
1.0, so F1 is driven by recall, and recall rises with how many peaks a strategy
emits.  Ranking strategies therefore requires holding the peak count fixed.

Outputs (under ``analysis_extra/``):
  strategy_metrics_full.csv   -- every cutoff + peak-width + atlas-per-call
  strategy_rarefaction.csv    -- all strategies rarefied to a common peak count
  pas_distance_to_atlas.csv   -- binned distance to nearest atlas PAS
  strategy_overlap.csv        -- pairwise agreement between strategies
"""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

SWEEP = Path("/mnt/ssd2/Laugney_Aligned/peakatail_experiments/RERUN_2026-08_fixed")
GRID = SWEEP / "runs" / "grid"
ATLAS = Path("/mnt/ssd2/Laugney_Aligned/refs/polyasite_3.0_GRCh38_ensembl_sorted.bed")
OUT = SWEEP / "analysis_extra"
TMP = Path("/mnt/ssd2/Laugney_Aligned/rt")
OUT.mkdir(parents=True, exist_ok=True)
TMP.mkdir(parents=True, exist_ok=True)

STRATEGY = {
    "lg_annotate": "lambda_gradient",
    "lg_ip_filter": "lambda_gradient (ip=filter)",
    "lg_ip_off": "lambda_gradient (ip=off)",
    "lp_annotate": "lambda_poisson",
    "si_annotate": "sierra_iterative",
}
RUNS = list(STRATEGY)
N_REPLICATES = 200  # cheap now: pure in-memory resampling
RNG_SEED = 0
BINS = [0, 10, 25, 50, 100, 200, 500, 1000, 2000, 5000, 10000, np.inf]


def read_pasbed(run: str) -> pd.DataFrame:
    df = pd.read_csv(
        GRID / run / "pasbed.bed", sep="\t", header=None,
        names=["chrom", "start", "end", "pas_id", "score", "strand"],
        dtype={"chrom": str},
    )
    df["width"] = df["end"] - df["start"]
    return df


def sort_bed(df: pd.DataFrame, path: Path) -> Path:
    df.sort_values(["chrom", "start"]).to_csv(path, sep="\t", header=False, index=False)
    return path


def closest_distances(bed_path: Path, tag: str) -> np.ndarray:
    """Distance from every called PAS to its nearest atlas PAS.

    ``bedtools closest -d`` appends the distance as the FINAL column; the atlas
    contributes a variable number of columns, so the distance is extracted by
    awk on ``$NF`` rather than by a hardcoded index.
    """
    out = TMP / f"{tag}.dist"
    cmd = (
        f"bedtools closest -a {bed_path} -b {ATLAS} -d -t first | "
        f"awk -F'\\t' '{{print $NF}}' > {out}"
    )
    subprocess.run(["bash", "-c", cmd], check=True)
    d = np.loadtxt(out, dtype=np.int64)
    out.unlink(missing_ok=True)
    return d[d >= 0]


# ---------------------------------------------------------------------------
# 1. full metrics table
# ---------------------------------------------------------------------------
rows, beds = [], {}
for run in RUNS:
    bench = json.loads((GRID / run / "benchmark_vs_polyasite_v3.json").read_text())
    bed = read_pasbed(run)
    beds[run] = bed
    for cutoff, m in bench["cutoffs"].items():
        rows.append({
            "run": run, "strategy": STRATEGY[run], "cutoff_bp": int(cutoff),
            "n_predicted": bench["n_predicted"], "n_reference": bench["n_reference"],
            "precision": m["precision"], "recall": m["recall"], "f1": m["f1"],
            "matched_predicted": m["matched_predicted"],
            "matched_reference": m["matched_reference"],
            "atlas_pas_per_call": m["matched_reference"] / bench["n_predicted"],
            "mean_peak_width": float(bed["width"].mean()),
            "median_peak_width": float(bed["width"].median()),
        })
pd.DataFrame(rows).sort_values(["cutoff_bp", "f1"], ascending=[True, False]).to_csv(
    OUT / "strategy_metrics_full.csv", index=False)
print("wrote strategy_metrics_full.csv", flush=True)

# ---------------------------------------------------------------------------
# 2. one distance pass per strategy -> distance histogram + rarefaction
# ---------------------------------------------------------------------------
distances: dict[str, np.ndarray] = {}
dist_rows = []
for run, bed in beds.items():
    # cache the distance vector: it is the only expensive step, and rerunning the
    # downstream summaries must not require another pass over the 1.2 GB atlas
    cache = OUT / f"_dist_{run}.npy"
    if cache.exists():
        d = np.load(cache)
    else:
        p = sort_bed(bed[["chrom", "start", "end", "pas_id", "score", "strand"]], TMP / f"{run}_full.bed")
        d = closest_distances(p, f"{run}_full")
        p.unlink(missing_ok=True)
        np.save(cache, d)
    distances[run] = d
    counts, _ = np.histogram(d, bins=BINS)
    for lo, hi, c in zip(BINS[:-1], BINS[1:], counts):
        dist_rows.append({
            "run": run, "strategy": STRATEGY[run], "bin_lo": lo, "bin_hi": hi,
            "count": int(c), "frac": float(c) / len(d), "n_total": len(d),
            "median_distance": float(np.median(d)),
        })
    print(f"  {run}: n={len(d)} median={np.median(d):.0f} "
          f"frac<=50bp={np.mean(d <= 50):.4f}", flush=True)

pd.DataFrame(dist_rows).to_csv(OUT / "pas_distance_to_atlas.csv", index=False)
print("wrote pas_distance_to_atlas.csv", flush=True)

# rarefy to the smallest *distance vector*, not the smallest BED: a handful of
# calls sit on contigs absent from the atlas and yield no distance at all
n_min = min(len(d) for d in distances.values())
print("rarefying to n =", n_min, flush=True)

rng = np.random.default_rng(RNG_SEED)
rare_rows = []
for run, d in distances.items():
    for rep in range(N_REPLICATES):
        sub = rng.choice(d, size=n_min, replace=False)
        rare_rows.append({
            "run": run, "strategy": STRATEGY[run], "replicate": rep, "n_sampled": n_min,
            "frac_within_50": float(np.mean(sub <= 50)),
            "frac_within_100": float(np.mean(sub <= 100)),
            "frac_within_500": float(np.mean(sub <= 500)),
            "median_distance": float(np.median(sub)),
            "mean_distance": float(np.mean(sub)),
        })
rare = pd.DataFrame(rare_rows)
rare.to_csv(OUT / "strategy_rarefaction_raw.csv", index=False)
rare.groupby(["run", "strategy", "n_sampled"]).agg(
    frac_within_50_mean=("frac_within_50", "mean"),
    frac_within_50_sd=("frac_within_50", "std"),
    frac_within_100_mean=("frac_within_100", "mean"),
    frac_within_500_mean=("frac_within_500", "mean"),
    median_distance_mean=("median_distance", "mean"),
    median_distance_sd=("median_distance", "std"),
).reset_index().sort_values("frac_within_50_mean", ascending=False).to_csv(
    OUT / "strategy_rarefaction.csv", index=False)
print("wrote strategy_rarefaction.csv", flush=True)

# ---------------------------------------------------------------------------
# 3. pairwise strategy overlap (called PAS within 100 bp of each other)
# ---------------------------------------------------------------------------
ov_rows = []
for i, a in enumerate(RUNS):
    for b in RUNS[i + 1:]:
        pa = sort_bed(beds[a][["chrom", "start", "end", "pas_id", "score", "strand"]], TMP / f"ov_{a}.bed")
        pb = sort_bed(beds[b][["chrom", "start", "end", "pas_id", "score", "strand"]], TMP / f"ov_{b}.bed")
        def n_hit(x: Path, y: Path) -> int:
            r = subprocess.run(["bash", "-c", f"bedtools window -w 100 -a {x} -b {y} -u | wc -l"],
                               capture_output=True, text=True, check=True)
            return int(r.stdout.strip())
        na, nb = n_hit(pa, pb), n_hit(pb, pa)
        ov_rows.append({
            "run_a": a, "run_b": b, "n_a": len(beds[a]), "n_b": len(beds[b]),
            "n_a_matched_in_b": na, "n_b_matched_in_a": nb,
            "frac_a_in_b": na / len(beds[a]), "frac_b_in_a": nb / len(beds[b]),
            "jaccard_approx": (na + nb) / (len(beds[a]) + len(beds[b])),
        })
        pa.unlink(missing_ok=True); pb.unlink(missing_ok=True)
        print("  overlap", a, b, na, nb, flush=True)

pd.DataFrame(ov_rows).to_csv(OUT / "strategy_overlap.csv", index=False)
print("DONE s1b", flush=True)
