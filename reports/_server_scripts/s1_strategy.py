#!/usr/bin/env python
"""Count-controlled re-analysis of the peak-calling strategy comparison.

The headline strategy ranking in the cumulative analysis is by F1 against
PolyASite v3.  That reference holds 18,432,135 PAS clusters, so at every cutoff
we test, precision saturates near 1.0 (any call near a gene is within 50 bp of
*some* atlas entry) and F1 is driven almost entirely by recall.  Recall in turn
rises monotonically with how many peaks a strategy emits.  The ranking is
therefore at risk of being a peak-count ranking rather than a PAS-accuracy
ranking -- the same failure mode the project already documented for peak width.

This script quantifies that directly:

1. ``strategy_metrics_full.csv``  -- every cutoff, plus peak-width stats and
   two count-normalised statistics: atlas PAS recovered per called PAS, and
   precision (which is reported for completeness, not for ranking).
2. ``strategy_rarefaction.csv``   -- each strategy randomly subsampled down to
   the smallest strategy's n_PAS (20 replicates) and recall recomputed.  This
   is the count-controlled comparison: if the ranking flips or collapses here,
   the raw F1 ranking was a count artefact.
3. ``pas_distance_to_atlas.csv``  -- binned distance from each called PAS to
   its nearest atlas PAS, per strategy (the honest resolution statistic).
4. ``strategy_overlap.csv``       -- pairwise PAS agreement between strategies.

Writes only under ``analysis_extra/``.  No pipeline recomputation.
"""
from __future__ import annotations

import json
import subprocess
import tempfile
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
N_REPLICATES = 20
RNG_SEED = 0


def read_pasbed(run: str) -> pd.DataFrame:
    df = pd.read_csv(
        GRID / run / "pasbed.bed",
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "pas_id", "score", "strand"],
        dtype={"chrom": str},
    )
    df["width"] = df["end"] - df["start"]
    return df


def run_bedtools(args: list[str], stdout_path: Path) -> None:
    with open(stdout_path, "w") as fh:
        subprocess.run(args, stdout=fh, check=True)


def sort_bed(df: pd.DataFrame, path: Path) -> Path:
    df = df.sort_values(["chrom", "start"])
    df.to_csv(path, sep="\t", header=False, index=False)
    return path


# ---------------------------------------------------------------------------
# 1. full metrics table + peak width
# ---------------------------------------------------------------------------
rows = []
beds: dict[str, pd.DataFrame] = {}
for run in RUNS:
    bench = json.loads((GRID / run / "benchmark_vs_polyasite_v3.json").read_text())
    bed = read_pasbed(run)
    beds[run] = bed
    n_pred = bench["n_predicted"]
    n_ref = bench["n_reference"]
    for cutoff, m in bench["cutoffs"].items():
        rows.append(
            {
                "run": run,
                "strategy": STRATEGY[run],
                "cutoff_bp": int(cutoff),
                "n_predicted": n_pred,
                "n_reference": n_ref,
                "precision": m["precision"],
                "recall": m["recall"],
                "f1": m["f1"],
                "matched_predicted": m["matched_predicted"],
                "matched_reference": m["matched_reference"],
                # count-normalised: how much unique atlas signal each call buys
                "atlas_pas_per_call": m["matched_reference"] / n_pred,
                "mean_peak_width": float(bed["width"].mean()),
                "median_peak_width": float(bed["width"].median()),
            }
        )
metrics = pd.DataFrame(rows).sort_values(["cutoff_bp", "f1"], ascending=[True, False])
metrics.to_csv(OUT / "strategy_metrics_full.csv", index=False)
print("wrote strategy_metrics_full.csv", len(metrics))

# ---------------------------------------------------------------------------
# 2. rarefaction: subsample every strategy to the smallest n_PAS
# ---------------------------------------------------------------------------
n_min = min(len(b) for b in beds.values())
print("rarefying all strategies to n =", n_min)

atlas_sorted = TMP / "atlas_sorted.bed"
if not atlas_sorted.exists():
    # PolyASite v3 bed is already coordinate sorted; reuse directly.
    atlas_sorted = ATLAS

rare_rows = []
rng = np.random.default_rng(RNG_SEED)
for run, bed in beds.items():
    for rep in range(N_REPLICATES):
        idx = rng.choice(len(bed), size=n_min, replace=False)
        sub = bed.iloc[np.sort(idx)]
        with tempfile.NamedTemporaryFile(suffix=".bed", dir=TMP, delete=False) as tf:
            sub_path = Path(tf.name)
        sort_bed(sub[["chrom", "start", "end", "pas_id", "score", "strand"]], sub_path)
        out_path = sub_path.with_suffix(".closest")
        # -d gives distance to nearest atlas feature for each called PAS
        run_bedtools(
            ["bedtools", "closest", "-a", str(sub_path), "-b", str(atlas_sorted), "-d", "-t", "first"],
            out_path,
        )
        dist = pd.read_csv(out_path, sep="\t", header=None, usecols=[-1], names=["d"])["d"]
        dist = dist[dist >= 0]
        rare_rows.append(
            {
                "run": run,
                "strategy": STRATEGY[run],
                "replicate": rep,
                "n_sampled": n_min,
                "frac_within_50": float((dist <= 50).mean()),
                "frac_within_100": float((dist <= 100).mean()),
                "frac_within_500": float((dist <= 500).mean()),
                "median_distance": float(dist.median()),
                "mean_distance": float(dist.mean()),
            }
        )
        sub_path.unlink(missing_ok=True)
        out_path.unlink(missing_ok=True)
    print("  rarefied", run)

rare = pd.DataFrame(rare_rows)
rare.to_csv(OUT / "strategy_rarefaction_raw.csv", index=False)
rare_summary = (
    rare.groupby(["run", "strategy", "n_sampled"])
    .agg(
        frac_within_50_mean=("frac_within_50", "mean"),
        frac_within_50_sd=("frac_within_50", "std"),
        frac_within_100_mean=("frac_within_100", "mean"),
        median_distance_mean=("median_distance", "mean"),
        median_distance_sd=("median_distance", "std"),
    )
    .reset_index()
    .sort_values("frac_within_50_mean", ascending=False)
)
rare_summary.to_csv(OUT / "strategy_rarefaction.csv", index=False)
print("wrote strategy_rarefaction.csv")

# ---------------------------------------------------------------------------
# 3. distance-to-atlas distribution on the FULL call set, per strategy
# ---------------------------------------------------------------------------
BINS = [0, 10, 25, 50, 100, 200, 500, 1000, 2000, 5000, 10000, np.inf]
dist_rows = []
for run, bed in beds.items():
    sub_path = TMP / f"{run}_full.bed"
    sort_bed(bed[["chrom", "start", "end", "pas_id", "score", "strand"]], sub_path)
    out_path = TMP / f"{run}_full.closest"
    run_bedtools(
        ["bedtools", "closest", "-a", str(sub_path), "-b", str(atlas_sorted), "-d", "-t", "first"],
        out_path,
    )
    dist = pd.read_csv(out_path, sep="\t", header=None, usecols=[-1], names=["d"])["d"]
    dist = dist[dist >= 0]
    counts, _ = np.histogram(dist, bins=BINS)
    for lo, hi, c in zip(BINS[:-1], BINS[1:], counts):
        dist_rows.append(
            {
                "run": run,
                "strategy": STRATEGY[run],
                "bin_lo": lo,
                "bin_hi": hi,
                "count": int(c),
                "frac": float(c) / len(dist),
                "n_total": len(dist),
                "median_distance": float(dist.median()),
            }
        )
    sub_path.unlink(missing_ok=True)
    out_path.unlink(missing_ok=True)
    print("  distances", run, "median", float(dist.median()))

pd.DataFrame(dist_rows).to_csv(OUT / "pas_distance_to_atlas.csv", index=False)
print("wrote pas_distance_to_atlas.csv")

# ---------------------------------------------------------------------------
# 4. pairwise strategy overlap (called PAS within 100 bp of each other)
# ---------------------------------------------------------------------------
ov_rows = []
for a in RUNS:
    for b in RUNS:
        if a >= b:
            continue
        pa = TMP / f"ov_{a}.bed"
        pb = TMP / f"ov_{b}.bed"
        sort_bed(beds[a][["chrom", "start", "end", "pas_id", "score", "strand"]], pa)
        sort_bed(beds[b][["chrom", "start", "end", "pas_id", "score", "strand"]], pb)
        out_path = TMP / f"ov_{a}_{b}.txt"
        run_bedtools(["bedtools", "window", "-w", "100", "-a", str(pa), "-b", str(pb), "-u"], out_path)
        n_a_in_b = sum(1 for _ in open(out_path))
        run_bedtools(["bedtools", "window", "-w", "100", "-a", str(pb), "-b", str(pa), "-u"], out_path)
        n_b_in_a = sum(1 for _ in open(out_path))
        ov_rows.append(
            {
                "run_a": a,
                "run_b": b,
                "n_a": len(beds[a]),
                "n_b": len(beds[b]),
                "n_a_matched_in_b": n_a_in_b,
                "n_b_matched_in_a": n_b_in_a,
                "frac_a_in_b": n_a_in_b / len(beds[a]),
                "frac_b_in_a": n_b_in_a / len(beds[b]),
                "jaccard_approx": (n_a_in_b + n_b_in_a) / (len(beds[a]) + len(beds[b])),
            }
        )
        for p in (pa, pb, out_path):
            p.unlink(missing_ok=True)
        print("  overlap", a, b)

pd.DataFrame(ov_rows).to_csv(OUT / "strategy_overlap.csv", index=False)
print("wrote strategy_overlap.csv")
print("DONE s1")
