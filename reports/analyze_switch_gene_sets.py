#!/usr/bin/env python
"""Cancer-gene and pathway characterisation of the recurrent 3'UTR switch genes.

Two questions the biological-findings report raises but never tested against a
proper background:

1. Are the genes whose 3'UTR usage trends across stage enriched for known cancer
   genes (COSMIC Cancer Gene Census / OncoKB)?
2. Do they share pathway themes (MSigDB Hallmark)?

The confound that governs both
-----------------------------
``technical_report.md`` §12.4 shows that recurrence across cell types tracks
**detectability**, not shared regulation: a gene trends in many cell types largely
when it is *detected* in many cell types.  Detectability is driven by expression
level, and cancer-gene catalogues and Hallmark sets are themselves biased toward
well-studied, highly expressed genes.  A naive "recurrent switch genes vs all
genes" enrichment would therefore recover expression bias and report it as cancer
biology.

Two guards are applied throughout:

* **Universe = genes that were actually tested** (those carrying at least one
  stage-trend call), never the whole genome.
* **Detectability is adjusted for**, via logistic regression on
  ``n_celltypes_tested`` and via a detectability-matched resampling control.
  Both the crude and the adjusted effect are reported, so the gap between them is
  visible rather than hidden.

Reference sets are downloaded once and cached under
``reports/cumulative_analysis/refs/`` so the analysis is reproducible offline:

* OncoKB cancer gene list (public API) -- carries the ``sangerCGC`` flag, i.e.
  membership of the COSMIC Cancer Gene Census, plus ONCOGENE/TSG typing.
* MSigDB Hallmark gene sets (50 sets, symbols).

Usage
-----
    python reports/analyze_switch_gene_sets.py [--extra DIR] [--out DIR]
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
EXTRA = HERE / "cumulative_analysis" / "extra"
REFS = HERE / "cumulative_analysis" / "refs"
OUT = EXTRA

ONCOKB_URL = "https://www.oncokb.org/api/v1/utils/cancerGeneList"
HALLMARK_URL = (
    "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/"
    "h.all.v2024.1.Hs.symbols.gmt"
)
RECURRENT_MIN = 3  # cell types; matches technical_report.md §12.4
RNG_SEED = 0
N_MATCHED_DRAWS = 2000


def fetch(url: str, dest: Path) -> Path:
    """Download once, then reuse the cache."""
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=120) as r:
        dest.write_bytes(r.read())
    print(f"  downloaded {dest.name} ({dest.stat().st_size:,} bytes)")
    return dest


def load_cancer_genes() -> pd.DataFrame:
    raw = json.loads(fetch(ONCOKB_URL, REFS / "oncokb_cancer_gene_list.json").read_text())
    df = pd.DataFrame(raw)[["hugoSymbol", "sangerCGC", "oncokbAnnotated", "vogelstein", "geneType"]]
    return df.rename(columns={"hugoSymbol": "symbol"}).drop_duplicates("symbol")


def load_hallmark() -> dict[str, set[str]]:
    text = fetch(HALLMARK_URL, REFS / "msigdb_hallmark.gmt").read_text()
    sets: dict[str, set[str]] = {}
    for line in text.strip().split("\n"):
        f = line.split("\t")
        if len(f) > 2:
            sets[f[0]] = set(f[2:])
    return sets


def fisher_2x2(a: int, b: int, c: int, d: int) -> tuple[float, float]:
    """Odds ratio and two-sided p for [[a,b],[c,d]]; Haldane correction on zeros."""
    table = [[a, b], [c, d]]
    _, p = stats.fisher_exact(table)
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    return (a * d) / (b * c), p


def bh_fdr(p: np.ndarray) -> np.ndarray:
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n, float)
    q[order] = np.minimum.accumulate((p[order] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)


def main(argv: list[str] | None = None) -> int:
    global EXTRA, REFS, OUT
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--extra", type=Path, default=EXTRA)
    ap.add_argument("--refs", type=Path, default=REFS)
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    EXTRA, REFS, OUT = args.extra, args.refs, args.out
    OUT.mkdir(parents=True, exist_ok=True)

    rec_path = EXTRA / "recurrence_decomposition.csv"
    if not rec_path.exists():
        print(f"missing {rec_path}; run reports/_server_scripts/s3_bio.py first")
        return 1

    print("[refs]")
    cancer = load_cancer_genes()
    hallmark = load_hallmark()
    cgc = set(cancer.loc[cancer["sangerCGC"] == True, "symbol"])            # noqa: E712
    onco = set(cancer.loc[cancer["oncokbAnnotated"] == True, "symbol"])     # noqa: E712
    onco_g = set(cancer.loc[cancer["geneType"] == "ONCOGENE", "symbol"])
    tsg = set(cancer.loc[cancer["geneType"] == "TSG", "symbol"])
    print(f"  COSMIC CGC {len(cgc)}, OncoKB-annotated {len(onco)}, "
          f"oncogenes {len(onco_g)}, TSGs {len(tsg)}, Hallmark sets {len(hallmark)}")

    g = pd.read_csv(rec_path)
    g = g[g["symbol"].notna() & (g["symbol"] != "")].copy()
    g["is_recurrent"] = g["n_celltypes_trending"] >= RECURRENT_MIN
    g["in_cgc"] = g["symbol"].isin(cgc)
    g["in_oncokb"] = g["symbol"].isin(onco)
    g["is_oncogene"] = g["symbol"].isin(onco_g)
    g["is_tsg"] = g["symbol"].isin(tsg)
    print(f"\n[universe] {len(g):,} tested genes with a symbol; "
          f"{g['is_recurrent'].sum():,} recurrent (>= {RECURRENT_MIN} cell types)")

    # ---------------------------------------------------------------- overlap
    rows = []
    for label, member in [("COSMIC_CGC", "in_cgc"), ("OncoKB_annotated", "in_oncokb"),
                          ("OncoKB_oncogene", "is_oncogene"), ("OncoKB_TSG", "is_tsg")]:
        rec, priv = g[g["is_recurrent"]], g[~g["is_recurrent"]]
        a, b = int(rec[member].sum()), int((~rec[member]).sum())
        c, d = int(priv[member].sum()), int((~priv[member]).sum())
        or_crude, p_crude = fisher_2x2(a, b, c, d)

        # adjust for detectability: how many cell types even tested the gene
        try:
            import statsmodels.api as sm
            X = sm.add_constant(pd.DataFrame({
                "is_recurrent": g["is_recurrent"].astype(float),
                "n_celltypes_tested": g["n_celltypes_tested"].astype(float),
            }))
            fit = sm.Logit(g[member].astype(float), X).fit(disp=0)
            or_adj = float(np.exp(fit.params["is_recurrent"]))
            p_adj = float(fit.pvalues["is_recurrent"])
        except Exception as exc:  # statsmodels optional
            print(f"  (logistic adjustment unavailable: {type(exc).__name__})")
            or_adj = p_adj = float("nan")

        # Detectability-matched resampling: draw a private-gene set whose
        # n_celltypes_tested distribution matches the recurrent set.
        #
        # This control is only meaningful if the two groups actually overlap in
        # detectability.  Here they barely do -- recurrent genes have a median
        # n_celltypes_tested of 18 against 2 for private genes, and the upper
        # strata hold thousands of recurrent genes against a handful of private
        # ones.  Resampling with replacement from a stratum of 2-3 genes yields a
        # "null" fixed by those few genes, so the estimate is reported together
        # with its effective support and flagged unreliable when support is thin.
        rng = np.random.default_rng(RNG_SEED)
        pool = {k: v[member].to_numpy() for k, v in priv.groupby("n_celltypes_tested")}
        want = rec["n_celltypes_tested"].value_counts()
        weighted_support, weight_total, min_stratum = 0.0, 0, np.inf
        for n_tested, k in want.items():
            arr = pool.get(n_tested)
            size = 0 if arr is None else len(arr)
            weighted_support += size * int(k)
            weight_total += int(k)
            min_stratum = min(min_stratum, size)
        eff_support = weighted_support / max(weight_total, 1)
        draws = []
        for _ in range(N_MATCHED_DRAWS):
            hits = tot = 0
            for n_tested, k in want.items():
                arr = pool.get(n_tested)
                if arr is None or not len(arr):
                    continue
                hits += rng.choice(arr, size=int(k), replace=True).sum()
                tot += int(k)
            if tot:
                draws.append(hits / tot)
        matched = float(np.mean(draws)) if draws else float("nan")
        rate_rec = a / max(a + b, 1)
        emp_p = (float(np.mean(np.array(draws) >= rate_rec)) if draws else float("nan"))
        matched_ok = bool(eff_support >= 20 and min_stratum >= 5)

        rows.append({
            "gene_set": label, "n_recurrent": a + b, "n_recurrent_in_set": a,
            "rate_recurrent": rate_rec,
            "n_private": c + d, "n_private_in_set": c, "rate_private": c / max(c + d, 1),
            "odds_ratio_crude": or_crude, "p_crude": p_crude,
            "odds_ratio_adj_detectability": or_adj, "p_adj_detectability": p_adj,
            "rate_detectability_matched": matched, "p_empirical_matched": emp_p,
            "matched_effective_support": eff_support, "matched_min_stratum": float(min_stratum),
            "matched_control_reliable": matched_ok,
        })
        flag = "" if matched_ok else "  [matched control UNRELIABLE: "\
                                     f"eff. support {eff_support:.1f}, min stratum {min_stratum:.0f}]"
        print(f"  {label:18s} recurrent {a}/{a+b} ({rate_rec:.3%})  "
              f"private {c}/{c+d} ({c/max(c+d,1):.3%})  "
              f"OR {or_crude:.2f} (p={p_crude:.2g}) -> adj {or_adj:.2f} (p={p_adj:.2g}){flag}")

    ov = pd.DataFrame(rows)
    ov.to_csv(OUT / "cancer_gene_overlap.csv", index=False)
    print(f"\nwrote {OUT / 'cancer_gene_overlap.csv'}")

    # -------------------------------------------- direction vs oncogene/TSG
    dec = g[g["dominant_direction"] == "decreasing"]
    inc = g[g["dominant_direction"] == "increasing"]
    a, b = int(dec["is_oncogene"].sum()), int(dec["is_tsg"].sum())
    c, d = int(inc["is_oncogene"].sum()), int(inc["is_tsg"].sum())
    or_dir, p_dir = fisher_2x2(a, b, c, d)
    print(f"[direction] shortening: {a} oncogenes vs {b} TSGs; "
          f"lengthening: {c} vs {d}; OR {or_dir:.2f} p={p_dir:.3g}")
    pd.DataFrame([{
        "test": "oncogene_vs_TSG_by_direction",
        "shortening_oncogene": a, "shortening_TSG": b,
        "lengthening_oncogene": c, "lengthening_TSG": d,
        "odds_ratio": or_dir, "p_value": p_dir,
    }]).to_csv(OUT / "direction_oncogene_tsg.csv", index=False)

    # ------------------------------------------------------ hallmark themes
    universe = set(g["symbol"])
    rec_set = set(g.loc[g["is_recurrent"], "symbol"])
    hrows = []
    for name, members in hallmark.items():
        m = members & universe
        if len(m) < 5:
            continue
        a = len(m & rec_set)
        b = len(rec_set) - a
        c = len(m) - a
        d = len(universe) - len(rec_set) - c
        or_, p = fisher_2x2(a, b, c, d)
        hrows.append({
            "gene_set": name, "n_in_universe": len(m), "n_in_recurrent": a,
            "frac_of_set_recurrent": a / len(m),
            "expected_frac": len(rec_set) / len(universe),
            "odds_ratio": or_, "p_value": p,
        })
    hall = pd.DataFrame(hrows)
    if not hall.empty:
        hall["q_value"] = bh_fdr(hall["p_value"].to_numpy())
        hall = hall.sort_values("p_value")
        hall.to_csv(OUT / "hallmark_enrichment.csv", index=False)
        n_sig = int((hall["q_value"] < 0.05).sum())
        print(f"\n[hallmark] {len(hall)} sets tested, {n_sig} at FDR<0.05")
        print(hall.head(8)[["gene_set", "n_in_universe", "n_in_recurrent",
                            "frac_of_set_recurrent", "odds_ratio", "q_value"]]
              .to_string(index=False))
        print(f"wrote {OUT / 'hallmark_enrichment.csv'}")

    # annotate the recurrent table itself so the report can name genes safely
    g[["gene_id", "symbol", "n_celltypes_trending", "n_celltypes_tested",
       "dominant_direction", "mean_slope", "in_cgc", "in_oncokb",
       "is_oncogene", "is_tsg"]].sort_values(
        "n_celltypes_trending", ascending=False).to_csv(
        OUT / "recurrent_genes_annotated.csv", index=False)
    print(f"wrote {OUT / 'recurrent_genes_annotated.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
