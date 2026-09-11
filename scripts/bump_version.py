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
    # Whitespace-tolerant: the BibTeX entry aligns its `=`, and a fixed
    # single-space pattern silently stopped matching when it was reformatted.
    ("README.md", [r'(version\s*=\s*\{)([0-9][^}]*)(\})']),
    ("docs/tutorials/01-installation.md", [r'(peakatail, version )([0-9][0-9A-Za-z.+-]*)(\s)']),
    ("docs/concepts/output-files.md", [r'("peakatail_version": ")([^"]+)(")']),
    # uv.lock records the project's OWN version alongside every dependency's,
    # so the pattern is anchored to the peakatail package block -- a bare
    # version regex would match the 32 unrelated packages that happen to share
    # a number. It had already drifted a full minor release behind pyproject
    # before this entry existed.
    ("uv.lock", [r'(?ms)(^\[\[package\]\]\nname = "peakatail"\nversion = ")([^"]+)(")']),
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


def scan() -> tuple[dict[str, list[str]], list[str]]:
    """``({path: [versions found]}, [declared-but-missing paths])``."""
    found: dict[str, list[str]] = {}
    missing: list[str] = []
    for rel, patterns in SITES:
        p = ROOT / rel
        if not p.exists():
            # A declared mirror that no longer exists means the guard silently
            # stopped guarding it: --check happily reported "consistent across
            # 8 site(s)" after recipes/peakatail/meta.yaml was moved away.
            missing.append(rel)
            continue
        text = p.read_text()
        hits: list[str] = []
        for pat in patterns:
            hits.extend(m.group(2) for m in re.finditer(pat, text))
        found[rel] = hits
    return found, missing


def check() -> int:
    want = current_version()
    found, missing = scan()
    bad: list[str] = [f"  {rel}: declared in SITES but the file is missing"
                      for rel in missing]
    for rel, hits in found.items():
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
    n = sum(len(v) for v in found.values())
    print(f"version {want} consistent across {n} site(s) in {len(found)} file(s)")
    return 0


def bump(new: str) -> int:
    if not SEMVER.match(new):
        sys.exit(f"{new!r} is not a semver-looking version (expected e.g. 0.4.0)")
    old = current_version()
    changed = 0
    for rel, patterns in SITES:
        p = ROOT / rel
        if not p.exists():
            sys.exit(f"declared site is missing: {rel}. Remove it from SITES "
                     "deliberately, or restore the file -- silently skipping "
                     "leaves that mirror unguarded.")
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
    print(f"Next: collapse the CHANGELOG `## Unreleased` headings into ONE "
          f"`## {new} - <date>` section, then tag v{new} on main.")
    print("      (scripts/changelog_section.py joins duplicates if you sed "
          "them all, but one section reads better in the Release body.)")
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
