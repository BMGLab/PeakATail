#!/usr/bin/env python3
"""Print one version's section of ``CHANGELOG.md`` on stdout.

The release workflow feeds this to ``gh release create --notes-file``, so the
text this prints becomes the public GitHub Release body -- and, because Zenodo
mints a DOI from the Release, part of a citable archival record. Getting the
wrong section (or an empty one) is therefore not a cosmetic problem, which is
why this is a real script with tests rather than an inline shell one-liner.

Usage::

    python scripts/changelog_section.py X.Y.Z [--changelog PATH]

Exit status is 1 (with the reason on stderr) when no ``## <version>`` heading
exists, so a release whose CHANGELOG was never collapsed from ``Unreleased``
fails loudly instead of publishing empty notes.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


def extract_section(text: str, version: str) -> str:
    """Return the body under ``## <version>``, without the heading itself.

    The heading may carry a trailing date and/or title (``## 0.2.0 (2026-05-10)
    -- Product CLI``); an optional leading ``v`` is accepted. The body runs to
    the next ``## `` heading or to end of file.
    """
    pattern = re.compile(
        r"^##[ \t]+v?" + re.escape(version) + r"(?![0-9A-Za-z.+-]).*$",
        re.MULTILINE,
    )
    match = pattern.search(text)
    if match is None:
        raise LookupError(
            f"CHANGELOG.md has no '## {version}' heading. "
            f"Collapse the 'Unreleased' sections into '## {version}' before tagging."
        )
    rest = text[match.end():]
    nxt = re.search(r"^##[ \t]", rest, re.MULTILINE)
    body = rest[: nxt.start()] if nxt else rest
    return body.strip("\n")


def remaining_unreleased(text: str) -> list[int]:
    """Line numbers of every ``## Unreleased`` heading still in ``text``.

    Collapsing the release notes is a manual step and this CHANGELOG carries
    one ``## Unreleased`` section per merged branch -- two dozen of them in
    the 0.3 cycle. A
    check that only fires when *nothing* was collapsed is therefore the wrong
    check: the realistic mistake is collapsing the first heading and missing
    the other 22, which would publish to PyPI and mint a permanent Zenodo DOI
    whose notes describe a few percent of the release, with every workflow step
    green.
    """
    return [
        i for i, line in enumerate(text.splitlines(), start=1)
        if re.match(r"^##[ \t]+Unreleased\b", line)
    ]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("version", help="version number, with or without a leading 'v'")
    ap.add_argument(
        "--changelog",
        type=Path,
        default=REPO_ROOT / "CHANGELOG.md",
        help="path to CHANGELOG.md (default: the repo's own)",
    )
    ap.add_argument(
        "--max-chars",
        type=int,
        default=0,
        help="truncate the notes to this many characters, appending a pointer "
             "to CHANGELOG.md. GitHub rejects a release body over 125,000 "
             "characters with 'Body is too long' -- and by then the PyPI "
             "upload has already happened. 0 disables truncation.",
    )
    ap.add_argument(
        "--allow-remaining-unreleased",
        action="store_true",
        help="proceed even if `## Unreleased` headings remain. Only for the "
             "rare case where a section was deliberately kept for the NEXT "
             "cycle; the default refusal is what stops a partially-collapsed "
             "CHANGELOG from becoming a permanent DOI record.",
    )
    args = ap.parse_args(argv)

    version = args.version.lstrip("v")
    try:
        body = extract_section(args.changelog.read_text(encoding="utf-8"), version)
    except LookupError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    if not body.strip():
        print(
            f"error: the '## {version}' section of {args.changelog} is empty.",
            file=sys.stderr,
        )
        return 1
    leftovers = remaining_unreleased(args.changelog.read_text(encoding="utf-8"))
    if leftovers and not args.allow_remaining_unreleased:
        print(
            f"error: {args.changelog} still has {len(leftovers)} "
            f"'## Unreleased' heading(s) at line(s) "
            f"{', '.join(str(n) for n in leftovers[:10])}"
            f"{' ...' if len(leftovers) > 10 else ''}.\n"
            f"       The '## {version}' section was found, but those sections "
            f"would be left OUT of the release notes -- and those notes become "
            f"the permanent Zenodo DOI record. Collapse them into "
            f"'## {version}' too, or pass --allow-remaining-unreleased if a "
            f"section is deliberately held for the next cycle.",
            file=sys.stderr,
        )
        return 1
    if args.max_chars and len(body) > args.max_chars:
        notice = (
            "\n\n---\n\n*These notes were truncated to fit GitHub's release-body "
            "limit. The complete changelog for this version is in "
            "[`CHANGELOG.md`](CHANGELOG.md).*"
        )
        keep = args.max_chars - len(notice)
        cut = body.rfind("\n", 0, keep)          # never split mid-line
        body = body[: cut if cut > 0 else keep] + notice
        print(
            f"note: notes truncated to {len(body)} chars "
            f"(limit {args.max_chars}); full text remains in CHANGELOG.md",
            file=sys.stderr,
        )

    print(body)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
