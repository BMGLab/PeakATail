"""Validate that every PDUI / diff / cluster-match strategy produces
output that satisfies its mathematical / biological invariants.

Run after the main pipeline so we have real clustered h5ad data:
    python3 scripts/validate_strategies.py emaout/per_dataset/sampleA/clusters.h5ad
"""

from __future__ import annotations
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from ema.quantification.strategies import get_pdui_strategy, list_pdui_strategies
from ema.switch_test.strategies import get_diff_strategy, list_diff_strategies
from ema.clustering.cross_dataset import get_match_strategy, list_match_strategies
from ema.quantification.marker_selector import select_marker_pas, restrict_count_matrix
from ema.annotate.gtf2isoform_utr import parse_isoform_utrs
from ema.quantification.pas_to_isoform import map_pas_to_isoforms
from ema.utils import ResourceManager


PASS = "\033[32mPASS\033[0m"
FAIL = "\033[31mFAIL\033[0m"
WARN = "\033[33mWARN\033[0m"


class CheckResult:
    def __init__(self):
        self.checks: list[tuple[str, str, str, str]] = []

    def add(self, strategy: str, check: str, status: str, detail: str = ""):
        self.checks.append((strategy, check, status, detail))
        print(f"  [{status}] {strategy} :: {check}  {detail}")

    def summary(self):
        n = len(self.checks)
        passes = sum(1 for _, _, s, _ in self.checks if s == PASS)
        fails = sum(1 for _, _, s, _ in self.checks if s == FAIL)
        warns = sum(1 for _, _, s, _ in self.checks if s == WARN)
        print()
        print("=" * 60)
        print(f"TOTAL: {n} checks — {passes} pass, {fails} fail, {warns} warn")
        if fails:
            print("\nFAILURES:")
            for strat, check, status, detail in self.checks:
                if status == FAIL:
                    print(f"  {strat}: {check} — {detail}")
        return fails == 0


def build_count_dfs(adata, marker_pas):
    X = adata.X.toarray() if sp.issparse(adata.X) else adata.X
    pas_index = [int(v) for v in adata.var_names]
    cell_index = list(adata.obs_names)
    pdui_full = pd.DataFrame(X.T, index=pas_index, columns=cell_index)
    diff_full = pd.DataFrame(X, index=cell_index, columns=pas_index)
    if marker_pas:
        pdui_df = restrict_count_matrix(pdui_full, marker_pas, axis="rows")
        diff_df = restrict_count_matrix(diff_full, marker_pas, axis="cols")
    else:
        pdui_df, diff_df = pdui_full, diff_full
    return pdui_df, diff_df


def check_pdui_classic(df, r):
    if df.empty:
        r.add("pdui:classic", "non-empty", FAIL, "0 rows")
        return
    r.add("pdui:classic", "non-empty", PASS, f"{len(df)} rows")
    pdui_col = next((c for c in df.columns if 'pdui' in c.lower()), None)
    if pdui_col is None:
        r.add("pdui:classic", "has pdui column", FAIL, f"cols={df.columns.tolist()}")
        return
    vals = df[pdui_col].dropna().astype(float)
    in_range = ((vals >= 0) & (vals <= 1)).all()
    r.add("pdui:classic", "PDUI ∈ [0,1]", PASS if in_range else FAIL,
          f"min={vals.min():.3f} max={vals.max():.3f}")


def check_pdui_proportion(df, r):
    if df.empty:
        r.add("pdui:proportion", "non-empty", FAIL, "0 rows")
        return
    r.add("pdui:proportion", "non-empty", PASS, f"{len(df)} rows")
    prop_col = next((c for c in df.columns if 'prop' in c.lower()), None)
    if prop_col is None:
        r.add("pdui:proportion", "has proportion column", FAIL, f"cols={df.columns.tolist()}")
        return
    vals = df[prop_col].dropna().astype(float)
    in_range = ((vals >= 0) & (vals <= 1)).all()
    r.add("pdui:proportion", "proportion ∈ [0,1]", PASS if in_range else FAIL,
          f"min={vals.min():.3f} max={vals.max():.3f}")
    # Sum-to-1 per (gene, cell) — pick a sample
    if 'gene_id' in df.columns and 'cell' in df.columns:
        sample = df.groupby(['gene_id', 'cell'])[prop_col].sum()
        nonzero = sample[sample > 0]
        if len(nonzero) > 0:
            sums_close_to_1 = ((nonzero - 1.0).abs() < 0.01).sum() / len(nonzero)
            r.add("pdui:proportion", "(gene,cell) sums to 1.0",
                  PASS if sums_close_to_1 > 0.95 else WARN,
                  f"{sums_close_to_1*100:.1f}% within ±0.01")


def check_pdui_shannon(df, r):
    if df.empty:
        r.add("pdui:shannon", "non-empty", FAIL, "0 rows")
        return
    r.add("pdui:shannon", "non-empty", PASS, f"{len(df)} rows")
    ent_col = next((c for c in df.columns if 'entropy' in c.lower() and 'norm' not in c.lower()), None)
    norm_col = next((c for c in df.columns if 'normalized_entropy' in c.lower()), None)
    if ent_col:
        vals = df[ent_col].dropna().astype(float)
        in_range = (vals >= 0).all()
        r.add("pdui:shannon", "entropy >= 0", PASS if in_range else FAIL,
              f"min={vals.min():.3f} max={vals.max():.3f}")
    if norm_col:
        vals = df[norm_col].dropna().astype(float)
        in_range = ((vals >= 0) & (vals <= 1.001)).all()
        r.add("pdui:shannon", "normalized_entropy ∈ [0,1]",
              PASS if in_range else FAIL,
              f"min={vals.min():.3f} max={vals.max():.3f}")


def check_diff(df, name, r):
    if df.empty:
        r.add(f"diff:{name}", "non-empty", FAIL, "0 rows")
        return
    r.add(f"diff:{name}", "non-empty", PASS, f"{len(df)} rows")
    if 'pvalue' not in df.columns:
        r.add(f"diff:{name}", "has pvalue", FAIL, f"cols={df.columns.tolist()}")
        return
    p = df['pvalue'].dropna().astype(float)
    p_ok = ((p >= 0) & (p <= 1)).all()
    r.add(f"diff:{name}", "pvalue ∈ [0,1]", PASS if p_ok else FAIL,
          f"min={p.min():.3g} max={p.max():.3g}")
    if 'qvalue' in df.columns:
        q = df['qvalue'].dropna().astype(float)
        q_ok = ((q >= 0) & (q <= 1)).all()
        r.add(f"diff:{name}", "qvalue ∈ [0,1]", PASS if q_ok else FAIL,
              f"min={q.min():.3g} max={q.max():.3g}")
        # qvalue >= pvalue is a property of BH FDR (mostly true for sorted)
        sig = (q < 0.05).sum()
        r.add(f"diff:{name}", "produces results",
              PASS if len(df) > 0 else FAIL,
              f"{sig} sig at q<0.05 ({100*sig/len(df):.1f}%)")


def check_match(df, name, expected_canonical, r):
    if df.empty:
        r.add(f"match:{name}", "non-empty", FAIL, "0 rows")
        return
    r.add(f"match:{name}", "non-empty", PASS, f"{len(df)} rows")
    if 'match_confidence' not in df.columns:
        r.add(f"match:{name}", "has match_confidence", FAIL,
              f"cols={df.columns.tolist()}")
        return
    c = df['match_confidence'].dropna().astype(float)
    c_ok = ((c >= 0) & (c <= 1.001)).all()
    r.add(f"match:{name}", "match_confidence ∈ [0,1]",
          PASS if c_ok else FAIL,
          f"min={c.min():.3f} max={c.max():.3f}")
    if 'canonical_cluster' in df.columns:
        n_canon = df['canonical_cluster'].nunique()
        # For 2 identical-input datasets: should match perfectly to K canonical
        ok = n_canon == expected_canonical
        r.add(f"match:{name}", f"canonical = {expected_canonical} (identical inputs)",
              PASS if ok else WARN,
              f"got {n_canon}")


def main():
    if len(sys.argv) < 2:
        print("Usage: validate_strategies.py <h5ad_path> [<h5ad_b for matching>]")
        sys.exit(1)
    h5ad_a = Path(sys.argv[1])
    h5ad_b = Path(sys.argv[2]) if len(sys.argv) > 2 else h5ad_a

    r = CheckResult()
    rm = ResourceManager()
    n_jobs = rm.get_n_jobs(per_worker_mb=300)
    print(f"Using n_jobs={n_jobs} (free RAM: {rm.free_ram_mb()} MB)")
    print()

    # Load adata + marker selection
    print("--- Loading & selecting markers ---")
    adata = ad.read_h5ad(h5ad_a)
    print(f"adata: {adata.shape[0]} cells × {adata.shape[1]} PAS")
    markers = select_marker_pas(adata, top_n_per_cluster=200)
    print(f"markers: {len(markers)}")
    pdui_df, diff_df = build_count_dfs(adata, markers)
    print(f"restricted: PDUI {pdui_df.shape}, diff {diff_df.shape}")

    # Load isoform map for PDUI
    print("\n--- Building isoform map ---")
    isoform_utrs = parse_isoform_utrs(Path("data/Homo_sapiens.GRCh38.99.gtf"))
    pas_isoform_map = map_pas_to_isoforms(Path("emaout/pasbed.bed"), isoform_utrs)
    print(f"isoform map: {len(pas_isoform_map)} PAS")

    cluster_labels = pd.Series(adata.obs['leiden'].values, index=adata.obs_names)

    # === PDUI strategies ===
    print(f"\n--- PDUI strategies ({list_pdui_strategies()}) ---")
    pdui_checkers = {
        'classic': check_pdui_classic,
        'proportion': check_pdui_proportion,
        'shannon': check_pdui_shannon,
    }
    for method in list_pdui_strategies():
        t0 = time.time()
        try:
            df = get_pdui_strategy(method).compute(
                count_matrix=pdui_df,
                pas_isoform_map=pas_isoform_map,
                aggregation='per_gene',
                isoform_collapse='none',
            )
            print(f"\n{method} ({time.time()-t0:.1f}s):")
            pdui_checkers[method](df, r)
        except Exception as e:
            r.add(f"pdui:{method}", "executes without error", FAIL, str(e)[:80])

    # === Diff strategies ===
    print(f"\n--- Diff strategies ({list_diff_strategies()}) ---")
    for method in list_diff_strategies():
        t0 = time.time()
        try:
            strat = get_diff_strategy(method)
            if strat.supports_multi_condition:
                df = strat.test(diff_df, cluster_labels, n_jobs=n_jobs)
            else:
                # 1 cluster pair only for speed
                df = strat.test(diff_df, cluster_labels,
                                cluster1='0', cluster2='1', n_jobs=n_jobs)
            print(f"\n{method} ({time.time()-t0:.1f}s):")
            check_diff(df, method, r)
        except Exception as e:
            r.add(f"diff:{method}", "executes without error", FAIL, str(e)[:80])

    # === Match strategies ===
    print(f"\n--- Match strategies ({list_match_strategies()}) ---")
    expected_canon = adata.obs['leiden'].nunique()
    for method in list_match_strategies():
        t0 = time.time()
        try:
            df = get_match_strategy(method).match(
                [h5ad_a, h5ad_b], ["sampleA", "sampleB"],
                n_top_markers=50, n_jobs=n_jobs,
            )
            print(f"\n{method} ({time.time()-t0:.1f}s):")
            check_match(df, method, expected_canon, r)
        except Exception as e:
            r.add(f"match:{method}", "executes without error", FAIL, str(e)[:80])

    return 0 if r.summary() else 1


if __name__ == "__main__":
    sys.exit(main())
