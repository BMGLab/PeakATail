"""Assert every old argparse flag has a Click home in the new CLI.

Run via:
    python3 scripts/audit_cli_coverage.py

Exit 0 = all old flags accounted for. Exit 1 = some flag was dropped
without explicit acknowledgement in the spec's §8 inventory.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

# Mapping per spec §8 — old argparse flag → new Click flag (or "DROPPED" for
# explicitly dead flags).
MAPPING = {
    # ema → ema run
    "--config": "--config",
    "--bamDir": "--bam-dir",
    "--bamFiles": "--bam-files",
    "--sequenceLen": "--seq-len",
    "--CellBarcodeLen": "--cb-len",
    "--BarcodeTag": "--barcode-tag",
    "--gtfDir": "--gtf",
    "--cell_combinations": "--cluster-pairs",
    "--threads": "--threads",
    "--bam-threads": "--bam-threads",
    "--pipeline": "--pipeline",
    "--batch-size": "--batch-size",
    "--tiles": "--tiles",
    "--tile-size": "--tile-size",
    "--tile-overlap": "--tile-overlap",
    "--strategy": "--peak-strategy",  # in run mode
    "--lambda-window": "--lambda-window",
    "--lambda-method": "--lambda-method",
    "--lambda-fold-change": "--lambda-fold-change",
    "--max-pas": "--max-pas",
    "--smoothing-window": "--smoothing-window",
    "--min-prominence": "--min-prominence",
    "--dynamic-threshold": "--dynamic-threshold",
    "--floor-threshold": "--floor-threshold",
    "--internal-priming-filter": "--ip-filter",
    "--genome-fasta": "--genome-fasta",
    "--annotation-filter": "--annot-filter",
    "--ip-a-stretch": "--ip-a-stretch",
    "--min-pas-per-cell": "--min-pas-per-cell",
    "--max-gene-distance": "--max-gene-distance",
    "--utr-multiplier": "--utr-multiplier",
    "--include-extended": "--include-extended",
    "--atlas": "--atlas",
    "--atlas-distance": "--atlas-distance",
    "--clustering-method": "--cluster-method",
    "--resolution": "--resolution",
    "--n-pcs": "--n-pcs",
    "--external-clusters": "--external-clusters",
    "--random-seed": "--random-seed",
    "--cluster-match-method": "--match-method",
    "--n-top-markers": "--n-top-markers",
    "--benchmark": "--benchmark",
    "--validate-db": "--validate-db",
    "--pdui-method": "DROPPED",
    "--pdui-isoform-agg": "DROPPED",
    "--pdui-isoform-collapse": "DROPPED",
    "--diff-method": "DROPPED",
    # ema_switch
    "--h5ad": "--h5ad",
    "--pasbed": "--pasbed",
    "--gtf": "--gtf",
    "--output-dir": "--output",
    "--cluster-pairs": "--cluster-pairs",
    "--cluster-key": "--cluster-key",
    "--marker-top-n": "--marker-top-n",
    "--marker-method": "--marker-method",
    "--fdr-threshold": "--fdr",
    "--max-jobs": "--threads",
    "--per-worker-mb": "--per-worker-mb",
    # ema_merge
    # (--bamFiles and --threads already mapped above)
    # ema_parse_gtf
    "--cache-dir": "--cache-dir",
}


def _scan_old_argparse() -> set[str]:
    """Find every --flag in the LEGACY argparse files."""
    legacy_files = [
        "ema/cli_legacy.py.bak",
        "ema/switch_test/cli.py",
        "ema/merge_bam/cli.py",
        "ema/annotate/cli.py",
    ]
    flags: set[str] = set()
    for path in legacy_files:
        p = Path(path)
        if not p.exists():
            continue
        for m in re.finditer(r'"(\-\-[a-zA-Z][a-zA-Z0-9_-]*)"', p.read_text()):
            flags.add(m.group(1))
    return flags


def _scan_new_click() -> set[str]:
    """Find every --flag in the NEW Click files."""
    flags: set[str] = set()
    for p in Path("ema/cli").glob("*.py"):
        for m in re.finditer(r'"(\-\-[a-zA-Z][a-zA-Z0-9_-]*)"', p.read_text()):
            flags.add(m.group(1))
        for m in re.finditer(r"'(\-\-[a-zA-Z][a-zA-Z0-9_-]*)'", p.read_text()):
            flags.add(m.group(1))
    return flags


def main() -> int:
    old = _scan_old_argparse()
    new = _scan_new_click()
    missing = []
    for o in sorted(old):
        target = MAPPING.get(o)
        if target is None:
            missing.append(f"  no MAPPING entry for {o!r}")
            continue
        if target == "DROPPED":
            continue
        if target not in new:
            missing.append(f"  {o!r} -> {target!r}: not found in new Click CLI")
    if missing:
        print("CLI coverage audit FAILED:")
        for m in missing:
            print(m)
        return 1
    print(f"CLI coverage audit OK: {len(old)} old flags accounted for "
          f"(new CLI exposes {len(new)} total flags).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
