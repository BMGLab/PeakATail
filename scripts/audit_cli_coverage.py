"""Assert every pre-cutover argparse flag has a Click home in the new CLI.

Run via:
    python3 scripts/audit_cli_coverage.py

Exit 0 = all old flags accounted for. Exit 1 = some flag was dropped
without explicit acknowledgement in the spec's §8 inventory.

Implementation note
-------------------
The legacy argparse files were deleted on commit ``fff57d0`` (see commit
message "chore: delete legacy argparse CLIs ..."). To keep this audit
honest after the deletion, the full set of argparse flags is captured in
:file:`tests/fixtures/legacy_flags.txt`. That file is the frozen source
of truth -- DO NOT edit it to silence the audit; instead add the
missing ``--flag`` to the relevant ``ema/cli/*.py`` module.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

# Mapping per spec §8 — old argparse flag → new Click flag (or "DROPPED" for
# explicitly dead flags).
MAPPING: dict[str, str] = {
    # ema → peakatail run
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


_REPO_ROOT = Path(__file__).resolve().parents[1]
LEGACY_SNAPSHOT = _REPO_ROOT / "tests" / "fixtures" / "legacy_flags.txt"


def _scan_old_argparse(snapshot: Path = LEGACY_SNAPSHOT) -> set[str]:
    """Read the frozen pre-cutover flag snapshot.

    The previous implementation searched deleted files on disk and so
    silently returned ``set()``, making the audit a no-op. The snapshot
    file is the authoritative pre-cutover flag list and is committed to
    the repo for exactly this reason.
    """
    if not snapshot.exists():
        raise FileNotFoundError(
            f"Legacy flag snapshot missing: {snapshot}.\n"
            "It is committed to the repo (tests/fixtures/legacy_flags.txt) "
            "and audit_cli_coverage.py cannot validate without it."
        )
    flags: set[str] = set()
    for line in snapshot.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if not line.startswith("--"):
            raise ValueError(
                f"{snapshot}: unexpected line {line!r} (expected --flag, "
                "comment, or blank)"
            )
        flags.add(line)
    if not flags:
        raise ValueError(
            f"{snapshot} contained no flags -- audit cannot be a no-op. "
            "Restore the snapshot from commit fff57d0~1."
        )
    return flags


def _scan_new_click(cli_dir: Path | None = None) -> set[str]:
    """Find every --flag in the NEW Click files."""
    cli_dir = cli_dir or (_REPO_ROOT / "ema" / "cli")
    flags: set[str] = set()
    for p in cli_dir.glob("*.py"):
        text = p.read_text()
        for m in re.finditer(r'"(\-\-[a-zA-Z][a-zA-Z0-9_-]*)"', text):
            flags.add(m.group(1))
        for m in re.finditer(r"'(\-\-[a-zA-Z][a-zA-Z0-9_-]*)'", text):
            flags.add(m.group(1))
    return flags


def audit() -> tuple[int, list[str]]:
    """Run the audit. Return (exit_code, list_of_failures)."""
    old = _scan_old_argparse()
    new = _scan_new_click()
    missing: list[str] = []
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
        return 1, missing
    return 0, []


def main() -> int:
    rc, missing = audit()
    if rc != 0:
        print("CLI coverage audit FAILED:")
        for m in missing:
            print(m)
        return rc
    old = _scan_old_argparse()
    new = _scan_new_click()
    print(
        f"CLI coverage audit OK: {len(old)} legacy flags accounted for "
        f"(new CLI exposes {len(new)} total flags)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
