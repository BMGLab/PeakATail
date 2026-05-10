"""Audit script must exit 0."""
import subprocess
import sys


def test_audit_passes():
    res = subprocess.run(
        [sys.executable, "scripts/audit_cli_coverage.py"],
        capture_output=True, text=True,
    )
    assert res.returncode == 0, f"audit failed:\n{res.stdout}\n{res.stderr}"
