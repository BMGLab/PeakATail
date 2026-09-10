#!/usr/bin/env python3
"""Bump (or verify) the project version everywhere it is written down.

`pyproject.toml` is the single source of truth. Every other site mirrors it,
and they drift the moment anyone edits one by hand -- which is exactly how
`CITATION.cff` and the docs ended up stale before. Run this instead:

    python scripts/bump_version.py --check      # verify every site agrees
    python scripts/bump_version.py 0.4.0        # rewrite every site

`--check` is what `tests/test_version_sync.py` runs, so a drifted file fails
CI rather than reaching a release.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# (path, [regex]) -- each regex must have exactly 3 groups: prefix, version, suffix.
SITES: list[tuple[str, list[str]]] = [
    ("pyproject.toml", [r'(?m)^(version = ")([^"]+)(")']),
    ("CITATION.cff", [r'(?m)^(version: ")([^"]+)(")']),
    ("recipes/peakatail/meta.yaml", [r'(\{%\s*set version = ")([^"]+)("\s*%\})']),
    ("README.md", [r'(version = \{)([0-9][^}]*)(\})']),
    ("docs/tutorials/01-installation.md", [r'(peakatail, version )([0-9][0-9A-Za-z.+-]*)(\s)']),
    ("docs/concepts/output-files.md", [r'("peakatail_version": ")([^"]+)(")']),
    # Illustrative release commands -- stale numbers here mislead the next releaser.
    ("CONTRIBUTING.md", [
        r'(git tag -a v)([0-9][0-9A-Za-z.+-]*)()',
        r'(git push origin v)([0-9][0-9A-Za-z.+-]*)()',
        r'(PeakATail v)([0-9][0-9A-Za-z.+-]*)()',
    ]),
]

SEMVER = re.compile(r"^\d+\.\d+\.\d+([0-9A-Za-z.+-]*)$")


def current_version() -> str:
    text = (ROOT / "pyproject.toml").read_text()
    m = re.search(r'(?m)^version = "([^"]+)"', text)
    if not m:
        sys.exit("pyproject.toml has no `version = \"...\"` line")
    return m.group(1)


def scan() -> dict[str, list[str]]:
    """{path: [versions found]} across every declared site."""
    found: dict[str, list[str]] = {}
    for rel, patterns in SITES:
        p = ROOT / rel
        if not p.exists():
            continue
        text = p.read_text()
        hits: list[str] = []
        for pat in patterns:
            hits.extend(m.group(2) for m in re.finditer(pat, text))
        found[rel] = hits
    return found


def check() -> int:
    want = current_version()
    bad: list[str] = []
    for rel, hits in scan().items():
        if not hits:
            bad.append(f"  {rel}: no version found (pattern stopped matching?)")
        for h in hits:
            if h != want:
                bad.append(f"  {rel}: {h!r} != pyproject {want!r}")
    if bad:
        print(f"version drift (pyproject.toml says {want!r}):", file=sys.stderr)
        print("\n".join(bad), file=sys.stderr)
        print("\nrun: python scripts/bump_version.py " + want, file=sys.stderr)
        return 1
    n = sum(len(v) for v in scan().values())
    print(f"version {want} consistent across {n} site(s) in {len(scan())} file(s)")
    return 0


def bump(new: str) -> int:
    if not SEMVER.match(new):
        sys.exit(f"{new!r} is not a semver-looking version (expected e.g. 0.4.0)")
    old = current_version()
    changed = 0
    for rel, patterns in SITES:
        p = ROOT / rel
        if not p.exists():
            print(f"  skip (missing): {rel}")
            continue
        text = orig = p.read_text()
        for pat in patterns:
            text = re.sub(pat, lambda m: m.group(1) + new + m.group(3), text)
        if text != orig:
            p.write_text(text)
            changed += 1
            print(f"  updated: {rel}")
        else:
            print(f"  unchanged: {rel}")
    print(f"\n{old} -> {new} across {changed} file(s)")
    print("Next: convert the CHANGELOG `## Unreleased` headings to "
          f"`## {new} - <date>`, then tag v{new} on main.")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("version", nargs="?", help="new version, e.g. 0.4.0")
    ap.add_argument("--check", action="store_true",
                    help="verify every site matches pyproject.toml; exit 1 on drift")
    a = ap.parse_args()
    if a.check:
        return check()
    if not a.version:
        ap.error("give a version to bump to, or --check")
    return bump(a.version)


if __name__ == "__main__":
    raise SystemExit(main())
