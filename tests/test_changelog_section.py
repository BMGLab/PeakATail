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
    # --allow-remaining-unreleased because this reads an ALREADY-RELEASED
    # section from a working tree whose next release has not been collapsed
    # yet; the strict default exists for the release path, not for reading
    # history.
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "v0.2.0", "--allow-remaining-unreleased"],
        capture_output=True, text=True, cwd=REPO_ROOT,
    )
    assert proc.returncode == 0, proc.stderr
    assert "Breaking changes" in proc.stdout
    assert "## 0.1" not in proc.stdout


# ---------------------------------------------------------------------------
# The realistic mistake is a PARTIAL collapse, not a forgotten one.
#
# This CHANGELOG accumulates one `## Unreleased` section per merged branch --
# 23 of them at 0.3.0. A maintainer who collapses the first heading and misses
# the rest gets a green PyPI publish, a green GitHub Release, and a permanent
# Zenodo DOI whose archival notes describe a few percent of the release. The
# original guard only fired when NOTHING had been collapsed, so it did not
# catch this at all.
# ---------------------------------------------------------------------------

def _partially_collapsed(tmp_path: Path) -> Path:
    changelog = tmp_path / "CHANGELOG.md"
    changelog.write_text(
        "# Changelog\n\n"
        "## 0.3.0 (2026-09-10)\n\n- the one section that was collapsed\n\n"
        "## Unreleased — forgotten A\n\n- a change that would vanish\n\n"
        "## Unreleased — forgotten B\n\n- another one\n\n"
        "## 0.2.0 (2026-05-10)\n\n- older\n",
        encoding="utf-8",
    )
    return changelog


def test_partial_collapse_is_refused(tmp_path: Path) -> None:
    changelog = _partially_collapsed(tmp_path)
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "0.3.0", "--changelog", str(changelog)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 1, (
        "a partially-collapsed CHANGELOG was accepted; the release notes (and "
        f"therefore the DOI record) would have omitted 2 sections.\n{proc.stdout}"
    )
    assert "Unreleased" in proc.stderr
    assert proc.stdout.strip() == "", "notes must not be emitted on refusal"


def test_partial_collapse_names_the_offending_lines(tmp_path: Path) -> None:
    """The error has to be actionable — 23 headings is a lot to hunt by eye."""
    changelog = _partially_collapsed(tmp_path)
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "0.3.0", "--changelog", str(changelog)],
        capture_output=True, text=True,
    )
    expected = [
        str(i) for i, line in enumerate(
            changelog.read_text(encoding="utf-8").splitlines(), start=1)
        if line.startswith("## Unreleased")
    ]
    assert len(expected) == 2, f"fixture drifted: {expected}"
    for lineno in expected:
        assert lineno in proc.stderr, (
            f"leftover heading at line {lineno} not named in stderr:\n{proc.stderr}"
        )


def test_deliberate_leftover_can_be_opted_out(tmp_path: Path) -> None:
    """A section genuinely held for the next cycle must still be releasable."""
    changelog = _partially_collapsed(tmp_path)
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "0.3.0", "--changelog", str(changelog),
         "--allow-remaining-unreleased"],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    assert "the one section that was collapsed" in proc.stdout
    assert "a change that would vanish" not in proc.stdout
