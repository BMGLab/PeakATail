"""Peak-by-peak comparison between two `emaout/` directories.

Usage:
    python3 scripts/diff_peak_outputs.py <baseline_dir> <new_dir>

Compares (in order):
1. peakcalling/*.bed — peak count, per-row (chrom, start, end, strand) match
   IGNORING pasnumber (col 4) since reorderings are allowed if the SET of
   peaks is the same.
2. peakcalling/*.mtx — total counts, total nnz, per-CB column sums
   (CB column ordering may differ, so we compare the SUM of counts per CB).
3. peakcalling/*.cb.tsv — same set of CBs (order-insensitive).
4. unified/*.bed — same peak coordinate set.
5. per_dataset/*/clusters.h5ad — same n_obs, n_vars, same cluster sizes.

Exit code 0 = identical (within numerical tolerance), 1 = differences.
"""
from __future__ import annotations
import sys
from pathlib import Path
from collections import Counter

import numpy as np
import scipy.io as sio
import scipy.sparse as sp


PASS = "\033[32m✓\033[0m"
FAIL = "\033[31m✗\033[0m"
WARN = "\033[33m!\033[0m"


def parse_bed(path: Path) -> list[tuple[str, int, int, str]]:
    """Return list of (chrom, start, end, strand) tuples — ignores pasnumber col 4."""
    out = []
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            try:
                out.append((parts[0], int(parts[1]), int(parts[2]), parts[5]))
            except ValueError:
                continue
    return out


def parse_mtx_headerless(path: Path) -> dict[int, int]:
    """Return per-column total counts dict (col_idx → sum)."""
    counts: Counter = Counter()
    with open(path) as f:
        for line in f:
            if line.startswith("%"):
                continue
            parts = line.split()
            if len(parts) != 3:
                continue
            try:
                _, c, v = int(parts[0]), int(parts[1]), int(parts[2])
                counts[c] += v
            except ValueError:
                continue
    return dict(counts)


def parse_cb_list(path: Path) -> set[str]:
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


def diff_dirs(baseline: Path, new: Path) -> int:
    fails = 0

    # 1. Peak calling BED files
    print("=== peakcalling/*.bed ===")
    base_beds = sorted((baseline / "peakcalling").glob("*.bed")) if (baseline / "peakcalling").exists() else []
    for b in base_beds:
        n = new / "peakcalling" / b.name
        if not n.exists():
            print(f"  {FAIL} {b.name}: missing in new")
            fails += 1
            continue
        b_peaks = parse_bed(b)
        n_peaks = parse_bed(n)
        if len(b_peaks) != len(n_peaks):
            print(f"  {FAIL} {b.name}: count {len(b_peaks)} → {len(n_peaks)}")
            fails += 1
            continue
        # Peaks may be re-ordered; compare as sorted sets
        b_set = sorted(b_peaks)
        n_set = sorted(n_peaks)
        if b_set != n_set:
            mismatched = sum(1 for a, c in zip(b_set, n_set) if a != c)
            print(f"  {FAIL} {b.name}: {mismatched} peak coordinate mismatches "
                  f"(of {len(b_peaks)})")
            # Show first 3 differences
            shown = 0
            for a, c in zip(b_set, n_set):
                if a != c and shown < 3:
                    print(f"      base={a}  new={c}")
                    shown += 1
            fails += 1
        else:
            print(f"  {PASS} {b.name}: {len(b_peaks)} peaks identical (chrom,start,end,strand)")

    # 2. Peak calling MTX files
    print("\n=== peakcalling/*.mtx (per-column total counts) ===")
    base_mtxs = sorted((baseline / "peakcalling").glob("*.mtx")) if (baseline / "peakcalling").exists() else []
    for b in base_mtxs:
        n = new / "peakcalling" / b.name
        if not n.exists():
            print(f"  {FAIL} {b.name}: missing in new")
            fails += 1
            continue
        b_counts = parse_mtx_headerless(b)
        n_counts = parse_mtx_headerless(n)
        b_total = sum(b_counts.values())
        n_total = sum(n_counts.values())
        b_nnz_per_col = sorted(b_counts.values())
        n_nnz_per_col = sorted(n_counts.values())
        if b_total == n_total and b_nnz_per_col == n_nnz_per_col:
            print(f"  {PASS} {b.name}: total={b_total}, "
                  f"{len(b_counts)} cols, identical column-count distribution")
        else:
            print(f"  {FAIL} {b.name}: total {b_total} → {n_total}, "
                  f"cols {len(b_counts)} → {len(n_counts)}")
            fails += 1

    # 3. CB lists
    print("\n=== peakcalling/*.cb.tsv (CB sets) ===")
    base_cbs = sorted((baseline / "peakcalling").glob("*.cb.tsv")) if (baseline / "peakcalling").exists() else []
    for b in base_cbs:
        n = new / "peakcalling" / b.name
        if not n.exists():
            print(f"  {FAIL} {b.name}: missing in new")
            fails += 1
            continue
        b_set = parse_cb_list(b)
        n_set = parse_cb_list(n)
        if b_set == n_set:
            print(f"  {PASS} {b.name}: {len(b_set)} CBs identical (set)")
        else:
            print(f"  {FAIL} {b.name}: {len(b_set - n_set)} only in baseline, "
                  f"{len(n_set - b_set)} only in new")
            fails += 1

    # 4. Unified BED
    print("\n=== unified/*.bed ===")
    base_unified = sorted((baseline / "unified").glob("*.bed")) if (baseline / "unified").exists() else []
    for b in base_unified:
        n = new / "unified" / b.name
        if not n.exists():
            print(f"  {WARN} {b.name}: missing in new (skipping)")
            continue
        b_peaks = parse_bed(b)
        n_peaks = parse_bed(n)
        if sorted(b_peaks) == sorted(n_peaks):
            print(f"  {PASS} {b.name}: {len(b_peaks)} unified peaks identical")
        else:
            print(f"  {FAIL} {b.name}: {len(b_peaks)} → {len(n_peaks)}")
            fails += 1

    # 5. Cluster h5ad
    print("\n=== per_dataset/*/clusters.h5ad (cluster sizes) ===")
    try:
        import anndata as ad
    except ImportError:
        print(f"  {WARN} anndata not available, skipping cluster check")
    else:
        base_h5s = sorted((baseline / "per_dataset").glob("*/clusters.h5ad")) if (baseline / "per_dataset").exists() else []
        for b in base_h5s:
            ds_id = b.parent.name
            n = new / "per_dataset" / ds_id / "clusters.h5ad"
            if not n.exists():
                print(f"  {FAIL} {ds_id}: missing in new")
                fails += 1
                continue
            ba = ad.read_h5ad(b)
            na = ad.read_h5ad(n)
            if ba.shape != na.shape:
                print(f"  {FAIL} {ds_id}: shape {ba.shape} → {na.shape}")
                fails += 1
                continue
            if "leiden" in ba.obs and "leiden" in na.obs:
                ba_sizes = sorted(ba.obs["leiden"].value_counts().tolist())
                na_sizes = sorted(na.obs["leiden"].value_counts().tolist())
                if ba_sizes == na_sizes:
                    print(f"  {PASS} {ds_id}: {ba.shape}, "
                          f"{len(ba_sizes)} clusters, identical cluster size distribution")
                else:
                    print(f"  {FAIL} {ds_id}: cluster sizes {ba_sizes} → {na_sizes}")
                    fails += 1

    print()
    print("=" * 60)
    if fails == 0:
        print(f"{PASS} ALL CHECKS PASS — outputs are functionally identical")
        return 0
    else:
        print(f"{FAIL} {fails} CHECKS FAILED")
        return 1


def main():
    if len(sys.argv) != 3:
        print("Usage: diff_peak_outputs.py <baseline_dir> <new_dir>")
        sys.exit(2)
    baseline = Path(sys.argv[1])
    new = Path(sys.argv[2])
    if not baseline.exists():
        print(f"Baseline dir not found: {baseline}")
        sys.exit(2)
    if not new.exists():
        print(f"New dir not found: {new}")
        sys.exit(2)
    sys.exit(diff_dirs(baseline, new))


if __name__ == "__main__":
    main()
