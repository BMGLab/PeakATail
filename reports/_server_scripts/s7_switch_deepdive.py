#!/usr/bin/env python
"""s7 — per-cell-type switch deep dive (both diff strategies, both length metrics).

Fills the gaps the aggregate figures leave: the report shows *how many* PAS
switch but not *which*, not the distribution of effects, and only one of the two
differential strategies per panel. This harvests, per cell type:

  1. ``diff_strategy_by_celltype.tsv`` — fisher AND nb_multi side by side, with
     the overlap of their significant PAS sets.
  2. ``volcano_points.tsv``            — delta_proportion vs q for every cell type
     (downsampled, but keeping every large-effect significant point).
  3. ``top_switch_genes.tsv``          — strongest switches per cell type x
     contrast, with gene symbols.
  4. ``shannon_by_stage.tsv``          — entropy per cell type x stage, streamed
     in awk (the entropy tables total ~15 GB).
  5. ``hallmark_switch_enrichment.tsv``— Hallmark enrichment of recurrent vs
     private switch-gene sets, against a *tested-gene* background.

Read-only against the real sweep. Nothing is written outside --out.

    python s7_switch_deepdive.py --run <B1_cohort_full> --out <dir> \
        --hallmark msigdb_hallmark.gmt
"""

from __future__ import annotations

import argparse
import logging
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("s7")

FDR = 0.05
# A switch has to be both significant and *large* to be worth naming. 43% of PAS
# clear the q threshold alone (pseudoreplication), so the effect floor is what
# actually selects.
MIN_DELTA = 0.30

# The primary differential strategy is a parameter, not a constant: `fisher` is
# being replaced by a per-cell/pseudobulk test, and this harvest must re-run
# against the replacement without a code edit. The layout it assumes is
#   <celltype>/<strategy>/differential/<strategy>_<contrast>.tsv
# with columns pas_id, gene_id, qvalue, delta_proportion, n_cells_cluster{1,2},
# log2fc. `--diff-strategy` switches it; `--secondary-strategy` is the omnibus
# arm compared against it (nb_multi today).
DEFAULT_DIFF = "fisher"
DEFAULT_SECONDARY = "nb_multi"


def gene_symbols(run: Path) -> pd.DataFrame:
    """gene_id -> gene_name from annotatedpas.bed (cols 5,6)."""
    bed = run / "annotatedpas.bed"
    df = pd.read_csv(
        bed, sep="\t", header=None, usecols=[4, 5], names=["gene_id", "gene_name"],
        dtype=str,
    )
    return df.drop_duplicates("gene_id")


def celltype_dirs(run: Path) -> list[Path]:
    d = run / "B3_switch" / "diff"
    return sorted(p for p in d.iterdir() if p.is_dir()) if d.exists() else []


# ---------------------------------------------------------------------------
# 1 + 2 + 3 — differential strategies, volcano points, top genes
# ---------------------------------------------------------------------------


def harvest_diff(run: Path, out: Path, sym: pd.DataFrame, rng_seed: int = 0,
                 diff: str = DEFAULT_DIFF, secondary: str = DEFAULT_SECONDARY,
                 fdr: float = FDR, min_delta: float = MIN_DELTA) -> None:
    rng = np.random.default_rng(rng_seed)
    strat_rows, volcano, top_rows = [], [], []

    for ct in celltype_dirs(run):
        name = ct.name
        nb_path = next((ct / secondary / "differential").glob("*.tsv"), None) \
                if (ct / secondary / "differential").exists() else None
        nb_sig: set[int] = set()
        nb_tests = nb_sig_n = 0
        if nb_path is not None and nb_path.exists():
            nb = pd.read_csv(nb_path, sep="\t", usecols=["pas_id", "qvalue"])
            nb_tests = len(nb)
            nb_sig = set(nb.loc[nb["qvalue"] <= fdr, "pas_id"].astype(int))
            nb_sig_n = len(nb_sig)

        for tsv in sorted((ct / diff / "differential").glob(f"{diff}_*.tsv")):
            contrast = tsv.stem.replace(f"{diff}_", "")
            try:
                f = pd.read_csv(
                    tsv, sep="\t",
                    usecols=["pas_id", "gene_id", "qvalue", "delta_proportion",
                             "n_cells_cluster1", "n_cells_cluster2", "log2fc"],
                )
            except (OSError, ValueError) as exc:
                log.warning("skip %s: %s", tsv, exc)
                continue
            if f.empty:
                continue
            f["pas_id"] = f["pas_id"].astype(int)
            sig = f["qvalue"] <= fdr
            big = sig & (f["delta_proportion"].abs() >= min_delta)
            f_sig = set(f.loc[sig, "pas_id"])

            inter = len(f_sig & nb_sig)
            union = len(f_sig | nb_sig)
            strat_rows.append({
                "celltype": name, "contrast": contrast,
                f"n_tests_{diff}": int(len(f)),
                f"n_sig_{diff}": int(sig.sum()),
                f"frac_sig_{diff}": round(float(sig.mean()), 4),
                f"n_sig_{diff}_large": int(big.sum()),
                f"frac_sig_{diff}_large": round(float(big.mean()), 4),
                f"n_tests_{secondary}": int(nb_tests),
                f"n_sig_{secondary}": int(nb_sig_n),
                f"frac_sig_{secondary}": round(nb_sig_n / nb_tests, 4) if nb_tests else np.nan,
                "n_sig_both": int(inter),
                "jaccard_sig": round(inter / union, 4) if union else np.nan,
                f"frac_{diff}_sig_also_{secondary}": round(inter / len(f_sig), 4) if f_sig else np.nan,
                "n_cells": int(
                    (f["n_cells_cluster1"].max() or 0) + (f["n_cells_cluster2"].max() or 0)
                ),
            })

            # volcano: keep every large-effect significant point, subsample the rest
            keep = f.index[big.to_numpy()]
            rest = f.index[~big.to_numpy()]
            if len(rest) > 2500:
                rest = rng.choice(rest, 2500, replace=False)
            v = f.loc[np.concatenate([keep.to_numpy(), np.asarray(rest)])].copy()
            v["celltype"], v["contrast"] = name, contrast
            v["neglog10q"] = -np.log10(np.clip(v["qvalue"], 1e-300, None))
            v["is_large_sig"] = v.index.isin(keep)
            volcano.append(
                v[["celltype", "contrast", "delta_proportion", "neglog10q",
                   "is_large_sig", "gene_id"]]
            )

            t = f.loc[big].nlargest(15, "delta_proportion", keep="all")
            b = f.loc[big].nsmallest(15, "delta_proportion", keep="all")
            for _, r in pd.concat([t, b]).iterrows():
                top_rows.append({
                    "celltype": name, "contrast": contrast,
                    "gene_id": r["gene_id"], "pas_id": int(r["pas_id"]),
                    "delta_proportion": round(float(r["delta_proportion"]), 4),
                    "log2fc": round(float(r["log2fc"]), 3),
                    "qvalue": float(r["qvalue"]),
                    "direction": "proximal-up" if r["delta_proportion"] > 0 else "distal-up",
                })
        log.info("diff harvested %s", name)

    pd.DataFrame(strat_rows).to_csv(out / "diff_strategy_by_celltype.tsv", sep="\t", index=False)
    if volcano:
        pd.concat(volcano, ignore_index=True).to_csv(
            out / "volcano_points.tsv", sep="\t", index=False
        )
    top = pd.DataFrame(top_rows)
    if not top.empty:
        top = top.merge(sym, on="gene_id", how="left")
        top["gene_name"] = top["gene_name"].fillna(top["gene_id"])
        top.to_csv(out / "top_switch_genes.tsv", sep="\t", index=False)
    log.info("wrote diff tables")


# ---------------------------------------------------------------------------
# 4 — shannon entropy per stage (streamed; the tables total ~15 GB)
# ---------------------------------------------------------------------------

AWK = r"""
BEGIN { FS = OFS = "\t" }
NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
{
    s = $c["cluster"]; e = $c["entropy"] + 0; ne = $c["normalized_entropy"] + 0
    t = $c["total_reads_gene"] + 0; n = $c["n_pas"] + 0
    rows[s]++; se[s] += e; sne[s] += ne; sn[s] += n
    if (t > 0) { inf[s]++; se_i[s] += e; sne_i[s] += ne }
    if (!((s SUBSEP $c["cell"]) in sc)) { sc[s SUBSEP $c["cell"]]; cells[s]++ }
}
END {
    for (s in rows)
        print s, rows[s], cells[s], se[s], sne[s], sn[s], inf[s] + 0, se_i[s] + 0, sne_i[s] + 0
}
"""


def harvest_shannon(run: Path, out: Path) -> None:
    rows = []
    base = run / "B3_switch" / "length"
    for ct in sorted(p for p in base.iterdir() if p.is_dir()) if base.exists() else []:
        f = ct / "shannon" / "entropy_shannon.tsv"
        if not f.exists():
            continue
        try:
            res = subprocess.run(["awk", AWK, str(f)], capture_output=True,
                                 text=True, check=True)
        except (subprocess.CalledProcessError, OSError) as exc:
            log.warning("awk failed on %s: %s", f, exc)
            continue
        for line in res.stdout.splitlines():
            p = line.split("\t")
            if len(p) != 9:
                continue
            stage, n, cells, se, sne, sn, inf, se_i, sne_i = p
            n, inf = float(n), float(inf)
            if n <= 0:
                continue
            rows.append({
                "celltype": ct.name, "stage": stage,
                "n_gene_cell_pairs": int(n), "n_cells": int(cells or 0),
                "mean_entropy": round(float(se) / n, 5),
                "mean_normalized_entropy": round(float(sne) / n, 5),
                "mean_n_pas": round(float(sn) / n, 3),
                "frac_uninformative": round(1 - inf / n, 5),
                "mean_entropy_informative": round(float(se_i) / inf, 5) if inf else np.nan,
                "mean_norm_entropy_informative": round(float(sne_i) / inf, 5) if inf else np.nan,
            })
        log.info("shannon aggregated %s", ct.name)
    pd.DataFrame(rows).to_csv(out / "shannon_by_stage.tsv", sep="\t", index=False)


# ---------------------------------------------------------------------------
# 5 — Hallmark enrichment, recurrent vs private, tested-gene background
# ---------------------------------------------------------------------------


def harvest_enrichment(run: Path, out: Path, hallmark: Path, sym: pd.DataFrame,
                       diff: str = DEFAULT_DIFF, fdr: float = FDR,
                       min_delta: float = MIN_DELTA) -> None:
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    top_path = out / "top_switch_genes.tsv"
    if not hallmark.exists():
        log.warning("no hallmark gmt at %s", hallmark)
        return

    # recurrence over *large-effect* switches only, and the background is the set
    # of genes actually tested -- using all genes would manufacture enrichment
    # for anything merely well-expressed.
    recs, tested = [], set()
    for ct in celltype_dirs(run):
        for tsv in sorted((ct / diff / "differential").glob(f"{diff}_*.tsv")):
            f = pd.read_csv(tsv, sep="\t", usecols=["gene_id", "qvalue", "delta_proportion"])
            tested.update(f["gene_id"].dropna().unique())
            hit = f[(f["qvalue"] <= fdr) & (f["delta_proportion"].abs() >= min_delta)]
            for g in hit["gene_id"].dropna().unique():
                recs.append({"gene_id": g, "celltype": ct.name})
    if not recs:
        log.warning("no switch genes")
        return
    rec = pd.DataFrame(recs).groupby("gene_id")["celltype"].nunique()
    private = set(rec[rec == 1].index)
    recurrent = set(rec[rec >= 3].index)

    s2n = dict(zip(sym["gene_id"], sym["gene_name"]))
    to_sym = lambda ids: {s2n.get(g) for g in ids if s2n.get(g)}
    bg, priv_s, rec_s = to_sym(tested), to_sym(private), to_sym(recurrent)

    sets = {}
    for line in hallmark.read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) > 2:
            sets[parts[0]] = set(parts[2:]) & bg

    rows = []
    for label, query in (("recurrent", rec_s), ("private", priv_s)):
        for name, members in sets.items():
            if len(members) < 5:
                continue
            a = len(query & members)
            b = len(query) - a
            c = len(members) - a
            d = len(bg) - a - b - c
            if min(a + b, a + c) == 0:
                continue
            odds, p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
            rows.append({
                "set_group": label, "gene_set": name, "n_set_in_bg": len(members),
                "n_query": len(query), "n_overlap": a,
                "odds_ratio": round(float(odds), 3), "pvalue": float(p),
            })
    df = pd.DataFrame(rows)
    if not df.empty:
        for g, idx in df.groupby("set_group").groups.items():
            df.loc[idx, "qvalue"] = multipletests(
                df.loc[idx, "pvalue"], method="fdr_bh"
            )[1]
        df = df.sort_values(["set_group", "pvalue"])
    df.to_csv(out / "hallmark_switch_enrichment.tsv", sep="\t", index=False)
    log.info("enrichment: %d recurrent, %d private genes, bg %d",
             len(rec_s), len(priv_s), len(bg))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--run", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--hallmark", type=Path, default=None)
    ap.add_argument("--diff-strategy", default=DEFAULT_DIFF,
                    help="primary differential strategy dir/prefix (e.g. mwu_percell)")
    ap.add_argument("--secondary-strategy", default=DEFAULT_SECONDARY)
    ap.add_argument("--fdr", type=float, default=FDR)
    ap.add_argument("--min-delta", type=float, default=MIN_DELTA)
    ap.add_argument("--skip-shannon", action="store_true")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    sym = gene_symbols(a.run)
    log.info("%d gene symbols", len(sym))
    harvest_diff(a.run, a.out, sym, diff=a.diff_strategy,
                 secondary=a.secondary_strategy, fdr=a.fdr, min_delta=a.min_delta)
    if not a.skip_shannon:
        harvest_shannon(a.run, a.out)
    if a.hallmark:
        harvest_enrichment(a.run, a.out, a.hallmark, sym, diff=a.diff_strategy,
                           fdr=a.fdr, min_delta=a.min_delta)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
