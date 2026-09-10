#!/usr/bin/env python3
"""Print one version's section of ``CHANGELOG.md`` on stdout.

The release workflow feeds this to ``gh release create --notes-file``, so the
text this prints becomes the public GitHub Release body -- and, because Zenodo
mints a DOI from the Release, part of a citable archival record. Getting the
wrong section (or an empty one) is therefore not a cosmetic problem, which is
why this is a real script with tests rather than an inline shell one-liner.

Usage::

    python scripts/changelog_section.py 0.3.0 [--changelog PATH]

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


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("version", help="version number, with or without a leading 'v'")
    ap.add_argument(
        "--changelog",
        type=Path,
        default=REPO_ROOT / "CHANGELOG.md",
        help="path to CHANGELOG.md (default: the repo's own)",
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
    print(body)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
