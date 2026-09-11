"""Every CLI flag the docs name must actually exist.

The docs told users to run `peakatail switch length --pdui-method proportion`
and `peakatail switch geneview --length-tsv ...`. Neither flag exists: the
first was renamed to `--strategy`, the second was removed. Both appeared in
copy-pasteable examples, so following the documentation produced an error --
and nothing noticed, because no test compared the documentation against the
CLI.

This walks the real Click command tree and fails on any backticked `--flag` in
docs/ that is not a real option somewhere.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

DOCS = Path(__file__).resolve().parents[1] / "docs"

# Flags that belong to OTHER tools, or to repo scripts rather than the CLI.
# Each needs a reason: the point is that adding one is a decision.
_FOREIGN_FLAGS = {
    "--check": "scripts/bump_version.py, not the peakatail CLI",
    "--rad": "alevin-fry, named when describing input formats",
    "--sketch": "alevin-fry",
    "--solo": "STARsolo",
    "--soloCBwhitelist": "STARsolo",
    # Documented precisely because it does NOT exist: run.md states "there is
    # no `--compat` flag; the list is a command line you type."
    "--compat": "named only to say it does not exist (docs/cli/run.md)",
}


def _real_flags() -> set[str]:
    """Every option string in the actual Click tree, walked recursively."""
    import click

    # `ema.cli.main` is the Click Group; `ema.cli.cli` is a plain wrapper
    # function, so walking that finds no options at all.
    from ema.cli import main as root

    found: set[str] = set()

    def walk(cmd: click.Command) -> None:
        for param in cmd.params:
            for opt in list(param.opts) + list(param.secondary_opts):
                if opt.startswith("--"):
                    found.add(opt)
        if isinstance(cmd, click.Group):
            for sub in cmd.commands.values():
                walk(sub)

    walk(root)
    # Click adds --help itself; it is not in cmd.params but is very real.
    found.add("--help")
    return found


def _documented_flags() -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    for path in DOCS.rglob("*.md"):
        for flag in re.findall(r"`(--[a-z0-9][a-z0-9-]+)", path.read_text(errors="ignore")):
            out.setdefault(flag, set()).add(str(path.relative_to(DOCS)))
    return out


def test_docs_do_not_name_flags_that_do_not_exist():
    real = _real_flags()
    assert real, "walked the CLI and found no options; the walker is broken"

    ghosts = {
        flag: files
        for flag, files in _documented_flags().items()
        if flag not in real and flag not in _FOREIGN_FLAGS
    }
    assert not ghosts, (
        "docs name flags the CLI does not have, so the documented command "
        "fails when a user runs it:\n"
        + "\n".join(f"  {f}  ({', '.join(sorted(w))})" for f, w in sorted(ghosts.items()))
        + "\n\nRename to the real flag, drop the reference, or add it to "
          "_FOREIGN_FLAGS with a reason if it belongs to another tool."
    )


def test_the_foreign_flag_allowlist_stays_honest():
    """An allowlist entry that IS a real flag hides a rename."""
    real = _real_flags()
    wrong = sorted(f for f in _FOREIGN_FLAGS if f in real)
    assert not wrong, (
        f"{wrong} are listed as foreign but are real peakatail flags; remove "
        "them from _FOREIGN_FLAGS so they are checked properly"
    )
