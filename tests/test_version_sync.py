"""The version must agree everywhere it is written down.

`pyproject.toml` is the source of truth; `CITATION.cff`, the bioconda recipe,
the README citation block and two docs pages mirror it. They drift the moment
someone edits one by hand -- which is exactly what happened before
`scripts/bump_version.py` existed. This test is that script's `--check` mode,
so drift fails CI instead of reaching a release (and PyPI, where a version
number can never be reused).
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "scripts" / "bump_version.py"


def test_bump_version_script_exists() -> None:
    assert SCRIPT.is_file(), f"missing {SCRIPT}"


def test_version_is_consistent_across_every_site() -> None:
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--check"],
        capture_output=True, text=True, cwd=REPO,
    )
    assert proc.returncode == 0, (
        "version drift between pyproject.toml and its mirrors:\n"
        f"{proc.stdout}{proc.stderr}"
    )


def test_check_actually_detects_drift(tmp_path: Path) -> None:
    """Guard on the guard: a --check that can never fail is worthless."""
    import re
    import shutil

    work = tmp_path / "repo"
    for rel in ("pyproject.toml", "CITATION.cff", "scripts/bump_version.py"):
        dst = work / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(REPO / rel, dst)

    cff = work / "CITATION.cff"
    cff.write_text(re.sub(r'(?m)^version: "[^"]+"', 'version: "9.9.9"', cff.read_text()))

    proc = subprocess.run(
        [sys.executable, str(work / "scripts" / "bump_version.py"), "--check"],
        capture_output=True, text=True, cwd=work,
    )
    assert proc.returncode == 1, "drifted CITATION.cff was not detected"
    assert "CITATION.cff" in proc.stderr
