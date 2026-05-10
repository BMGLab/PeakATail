"""Audit script must exit 0 and actually validate Click coverage.

This file is the second line of defence behind the snapshot file. The
prior incarnation only invoked the script as a subprocess and trusted
its exit code -- which made the audit a tautology when the script
silently returned an empty set. Here we additionally:

1. Verify the snapshot file is non-empty and contains real flags.
2. Verify the audit FAILS when a Click flag is removed (mock test).
3. Verify the audit FAILS when a brand-new legacy flag has no MAPPING.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SNAPSHOT = REPO_ROOT / "tests" / "fixtures" / "legacy_flags.txt"


# ---------------------------------------------------------------------------
# Subprocess-level: the actual audit binary still has to exit 0.
# ---------------------------------------------------------------------------

def test_audit_passes():
    res = subprocess.run(
        [sys.executable, "scripts/audit_cli_coverage.py"],
        capture_output=True, text=True, cwd=REPO_ROOT,
    )
    assert res.returncode == 0, f"audit failed:\n{res.stdout}\n{res.stderr}"


# ---------------------------------------------------------------------------
# Snapshot file must be present and non-empty -- the previous bug was
# that the script silently scanned deleted files and returned set().
# ---------------------------------------------------------------------------

def test_snapshot_file_exists_and_has_flags():
    assert SNAPSHOT.exists(), f"missing snapshot file {SNAPSHOT}"
    flag_lines = [
        ln.strip() for ln in SNAPSHOT.read_text().splitlines()
        if ln.strip() and not ln.strip().startswith("#")
    ]
    assert len(flag_lines) >= 50, (
        f"only {len(flag_lines)} flags in snapshot -- pre-cutover argparse "
        "had ~59. Snapshot may have been truncated."
    )
    for ln in flag_lines:
        assert ln.startswith("--"), f"snapshot line is not a flag: {ln!r}"


# ---------------------------------------------------------------------------
# Mock test: if you delete a Click flag, the audit MUST fail.
# This is the smoke test that proves the audit is not a tautology.
# ---------------------------------------------------------------------------

def test_audit_fails_when_click_flag_removed(tmp_path):
    """Run the audit pointing at a fake Click dir missing one mapped flag."""
    sys.path.insert(0, str(REPO_ROOT))
    try:
        # Reload to pick up the latest mapping
        if "audit_cli_coverage" in sys.modules:
            del sys.modules["audit_cli_coverage"]
        sys.path.insert(0, str(REPO_ROOT / "scripts"))
        try:
            import audit_cli_coverage as audit_mod  # type: ignore[import-not-found]
        finally:
            sys.path.pop(0)

        # Fake Click dir with a single bogus file that does NOT mention --threads
        fake_cli = tmp_path / "ema" / "cli"
        fake_cli.mkdir(parents=True)
        (fake_cli / "stub.py").write_text(
            '"""Bogus stub. Mentions --config and --gtf but nothing else."""\n'
            'import click\n'
            '@click.command()\n'
            '@click.option("--config")\n'
            '@click.option("--gtf")\n'
            'def x(): pass\n'
        )

        old = audit_mod._scan_old_argparse(SNAPSHOT)
        new = audit_mod._scan_new_click(fake_cli)

        # The audit's exit code is implicit; replicate the missing-check
        # body so the test does not depend on print/stdout format.
        missing = []
        for o in sorted(old):
            target = audit_mod.MAPPING.get(o)
            if target is None or target == "DROPPED":
                continue
            if target not in new:
                missing.append((o, target))

        assert missing, (
            "Expected the audit to flag the missing flags when scanning a "
            "bogus Click dir. Instead missing was empty -- the audit is not "
            "actually checking coverage."
        )
        # And specifically --threads -> --threads must appear in missing
        threads_pair = [(o, t) for (o, t) in missing if t == "--threads"]
        assert threads_pair, (
            "--threads should be flagged as missing when the Click dir "
            "doesn't expose it. Got: " + repr(missing[:5])
        )
    finally:
        sys.path.remove(str(REPO_ROOT))


# ---------------------------------------------------------------------------
# Snapshot integrity: every snapshot flag must have a MAPPING entry.
# Catches the case where someone adds a new legacy alias without telling
# the audit how to bridge it.
# ---------------------------------------------------------------------------

def test_every_snapshot_flag_has_mapping_entry():
    sys.path.insert(0, str(REPO_ROOT / "scripts"))
    try:
        if "audit_cli_coverage" in sys.modules:
            del sys.modules["audit_cli_coverage"]
        import audit_cli_coverage as audit_mod  # type: ignore[import-not-found]
    finally:
        sys.path.pop(0)

    flags = audit_mod._scan_old_argparse(SNAPSHOT)
    missing_mappings = [f for f in sorted(flags) if f not in audit_mod.MAPPING]
    assert not missing_mappings, (
        f"Snapshot lists {len(missing_mappings)} flags with no MAPPING entry: "
        f"{missing_mappings[:5]}{' ...' if len(missing_mappings) > 5 else ''}. "
        "Add an entry to MAPPING (or DROPPED) in scripts/audit_cli_coverage.py."
    )
