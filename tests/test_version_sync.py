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


# ---------------------------------------------------------------------------
# SITES must be COMPLETE, not merely correct.
#
# Every check in this file asks "do the declared sites agree?". None asks
# "is every site declared?" -- so a version hardcoded into a NEW file is
# invisible to the guard, and ships stale. uv.lock was exactly that: it pinned
# `peakatail version = "0.2.0"` while pyproject said 0.3.0, and --check
# reported "consistent across 9 site(s)".
# ---------------------------------------------------------------------------

# Files that legitimately mention the current version without being a mirror
# of it. Each entry needs a reason: the point is that adding one is a decision.
_VERSION_MENTION_ALLOWLIST = {
    "CHANGELOG.md": "release headings and historical entries",
    "docs/RELEASING.md": "SemVer illustrations (`0.3.0` -> `0.4.0`)",
    ".github/workflows/release.yml": "a comment example (`e.g. v0.3.0`)",
    "tests/test_version_sync.py": "this file's own allowlist",
    "docs/cli/switch-diff.md": (
        "historical references ('before 0.3.0 the NB offset was computed "
        "from the restricted matrix'), which date WHEN a behaviour changed "
        "and must NOT move with the version. They collide with this check "
        "only while 0.3.0 is the current version"
    ),
    "tests/test_label_independent_prefilter_i94.py": (
        "same: a docstring dating the NB library-size offset fix"
    ),
    "tests/test_changelog_section.py": (
        "sample version numbers in CHANGELOG fixtures -- these are test data, "
        "not mirrors of the project version, and pinning them to the real one "
        "would make the fixtures change every release"
    ),
}

# Lockfiles pin hundreds of third-party packages, any of which may coincide
# with our version number. uv.lock IS a declared site, but anchored to the
# peakatail block; a whole-file scan would drown in false positives.
_SCAN_SKIP_SUFFIXES = (".lock",)


def test_sites_covers_every_file_that_hardcodes_the_version():
    sys.path.insert(0, str(REPO / "scripts"))
    from bump_version import SITES, current_version

    version = current_version()
    declared = {rel for rel, _ in SITES}

    tracked = subprocess.run(
        ["git", "ls-files"], cwd=REPO,
        capture_output=True, text=True, check=True,
    ).stdout.split()

    undeclared: dict[str, int] = {}
    for rel in tracked:
        if rel in declared or rel in _VERSION_MENTION_ALLOWLIST:
            continue
        if rel.endswith(_SCAN_SKIP_SUFFIXES):
            continue
        path = REPO / rel
        if not path.is_file():
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        hits = text.count(version)
        if hits:
            undeclared[rel] = hits

    assert not undeclared, (
        f"these tracked files hardcode the current version ({version}) but are "
        f"not declared in scripts/bump_version.py's SITES, so a bump would "
        f"leave them stale while --check still reports 'consistent':\n"
        + "\n".join(f"  {rel} ({n}x)" for rel, n in sorted(undeclared.items()))
        + "\n\nAdd a SITES pattern, or add the file to "
          "_VERSION_MENTION_ALLOWLIST with a reason."
    )
