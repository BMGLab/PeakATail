"""`scripts/changelog_section.py` feeds the public GitHub Release body.

The release workflow pipes this script's stdout into
``gh release create --notes-file``, and Zenodo mints a DOI from that Release, so
what it prints ends up in a citable archival record. These tests pin the three
ways it could quietly produce the wrong record: bleeding into the next version's
section, matching a *prefix* of the version, and printing nothing at all when
the CHANGELOG was never collapsed out of ``Unreleased``.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "changelog_section.py"

sys.path.insert(0, str(REPO_ROOT / "scripts"))
from changelog_section import extract_section  # noqa: E402

SAMPLE = """# Changelog

## Unreleased — something in flight

- not part of any release yet

## 0.3.0 (2026-09-10) — Public release

### Fixed
- the thing

## 0.2.0 (2026-05-10) — Product CLI

- the older thing
"""


def test_section_stops_at_the_next_version() -> None:
    body = extract_section(SAMPLE, "0.3.0")
    assert "the thing" in body
    # The previous release's notes must not be appended to this one's.
    assert "the older thing" not in body
    # Nor may an in-flight Unreleased section leak in.
    assert "not part of any release yet" not in body
    # The version heading itself is not repeated in the body (it is already
    # the Release title), but sub-headings like "### Fixed" are kept.
    assert "## 0.3.0" not in body
    assert body.startswith("### Fixed")


def test_leading_v_is_accepted_and_a_prefix_is_not() -> None:
    assert extract_section(SAMPLE, "0.3.0") == extract_section(SAMPLE.replace("## 0.3.0", "## v0.3.0"), "0.3.0")
    # "0.3" must not silently return the 0.3.0 section: a tag typo would
    # otherwise publish the wrong notes under the wrong number.
    with pytest.raises(LookupError):
        extract_section(SAMPLE, "0.3")


def test_missing_version_exits_nonzero(tmp_path: Path) -> None:
    """A CHANGELOG still stuck on `Unreleased` must fail the release, not
    publish empty notes."""
    changelog = tmp_path / "CHANGELOG.md"
    changelog.write_text("# Changelog\n\n## Unreleased\n\n- pending\n", encoding="utf-8")
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "0.3.0", "--changelog", str(changelog)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 1, proc.stdout
    assert "0.3.0" in proc.stderr
    assert proc.stdout.strip() == ""


def test_runs_against_the_real_changelog() -> None:
    """0.2.0 is a released version that is in the file today; if this stops
    working the extractor has drifted from the CHANGELOG's actual heading
    style."""
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "v0.2.0"], capture_output=True, text=True, cwd=REPO_ROOT,
    )
    assert proc.returncode == 0, proc.stderr
    assert "Breaking changes" in proc.stdout
    assert "## 0.1" not in proc.stdout
