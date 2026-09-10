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


# ---------------------------------------------------------------------------
# Duplicate version headings: the way a sed-based collapse actually goes wrong.
#
# With ~two dozen `## Unreleased` sections, the obvious collapse is one sed over
# the headings -- which yields ~two dozen `## <version>` headings and ZERO
# leftover `Unreleased` ones. remaining_unreleased() is therefore silent, and
# returning only the first section published 1.8% of the release as the
# permanent DOI record with every step green.
# ---------------------------------------------------------------------------

def _sed_collapsed(tmp_path: Path) -> Path:
    changelog = tmp_path / "CHANGELOG.md"
    changelog.write_text(
        "# Changelog\n\n"
        "## 0.3.0 (2026-09-10)\n\n- first section\n\n"
        "## 0.3.0 (2026-09-10)\n\n- second section\n\n"
        "## 0.3.0 (2026-09-10)\n\n- third section\n\n"
        "## 0.2.0 (2026-05-10)\n\n- older\n",
        encoding="utf-8",
    )
    return changelog


def test_duplicate_version_headings_are_all_included(tmp_path: Path) -> None:
    changelog = _sed_collapsed(tmp_path)
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "0.3.0", "--changelog", str(changelog)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    for marker in ("first section", "second section", "third section"):
        assert marker in proc.stdout, (
            f"{marker!r} missing from the release notes -- only the first "
            f"`## 0.3.0` section was used. Got:\n{proc.stdout}"
        )
    assert "older" not in proc.stdout, "bled into the previous release"
    assert "joined 3" in proc.stderr, (
        f"the join should be reported on stderr, got: {proc.stderr!r}"
    )


def test_max_chars_is_actually_honoured(tmp_path: Path) -> None:
    """A limit-enforcer that does not enforce is worse than none.

    A negative `keep` made rfind() index from the END of the string, so the
    output silently exceeded the limit while stderr claimed otherwise.
    """
    changelog = _sed_collapsed(tmp_path)
    for limit in (500, 600, 900):
        proc = subprocess.run(
            [sys.executable, str(SCRIPT), "0.3.0", "--changelog", str(changelog),
             "--max-chars", str(limit)],
            capture_output=True, text=True,
        )
        assert proc.returncode == 0, proc.stderr
        assert len(proc.stdout.rstrip("\n")) <= limit, (
            f"--max-chars {limit} produced {len(proc.stdout)} characters"
        )


def test_truncation_never_leaves_an_open_code_fence(tmp_path: Path) -> None:
    """An unclosed fence renders the truncation notice AS CODE -- hiding the
    one line that says the notes are incomplete."""
    changelog = tmp_path / "CHANGELOG.md"
    changelog.write_text(
        "# Changelog\n\n## 0.3.0 (2026-09-10)\n\n"
        + "".join(f"- filler line {i}\n" for i in range(40))
        + "\n```python\n"
        + "".join(f"code_line_{i}()\n" for i in range(40))
        + "```\n\n## 0.2.0\n\n- older\n",
        encoding="utf-8",
    )
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "0.3.0", "--changelog", str(changelog),
         "--max-chars", "900"],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.count("```") % 2 == 0, (
        "truncated body left an odd number of code fences, so the trailing "
        f"notice renders inside a code block:\n{proc.stdout[-300:]}"
    )


def test_max_chars_below_the_notice_length_is_rejected(tmp_path: Path) -> None:
    """The guard that makes the truncation arithmetic safe.

    `keep = max_chars - len(notice)` goes negative for a tiny limit, and
    `str.rfind(sub, 0, negative)` indexes from the END rather than clamping --
    so the output silently EXCEEDED the limit while stderr reported it had been
    truncated to it. `max(0, ...)` bounds the arithmetic; this argparse guard is
    what stops a caller asking for a limit that cannot hold even the notice.
    """
    changelog = _sed_collapsed(tmp_path)
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "0.3.0", "--changelog", str(changelog),
         "--max-chars", "50"],
        capture_output=True, text=True,
    )
    assert proc.returncode != 0, (
        f"--max-chars 50 was accepted and produced {len(proc.stdout)} "
        "characters; a limit smaller than the truncation notice cannot be "
        "honoured and must be refused."
    )
    assert "at least 500" in proc.stderr
