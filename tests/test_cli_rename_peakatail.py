"""The command is `peakatail`; `ema` must keep working as a deprecated alias.

A rename that silently breaks `ema` breaks every published pipeline that calls it,
so these pin both halves: the new name is the one the tool reports, and the old
name still runs and says so on stderr.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
from click.testing import CliRunner

sys.path.insert(0, str(Path(__file__).parent.parent))

from ema.cli import ema_deprecated_main, main

ROOT = Path(__file__).parent.parent


try:                       # 3.11+
    import tomllib as _toml
except ModuleNotFoundError:  # 3.10, which the project venv still uses
    import tomli as _toml


def _scripts() -> dict[str, str]:
    with open(ROOT / "pyproject.toml", "rb") as fh:
        return _toml.load(fh)["project"]["scripts"]


def test_peakatail_is_an_installed_console_script():
    assert _scripts().get("peakatail") == "ema.cli:main"


def test_ema_alias_is_still_installed():
    """Removing it would break every existing script that calls `ema`."""
    assert "ema" in _scripts(), (
        "the `ema` console script must remain as a deprecated alias; dropping it "
        "is a breaking change for published pipelines"
    )
    assert _scripts()["ema"] == "ema.cli:ema_deprecated_main"


def test_group_reports_the_new_name():
    assert main.name == "peakatail"


def test_help_shows_peakatail_not_ema():
    out = CliRunner().invoke(main, ["--help"]).output
    assert "peakatail" in out
    assert "Usage: ema " not in out


def test_alias_warns_on_stderr_only(capsys, monkeypatch):
    """The notice must not contaminate stdout, which callers parse."""
    monkeypatch.setattr(sys, "argv", ["ema", "--help"])
    with pytest.raises(SystemExit):
        ema_deprecated_main()
    cap = capsys.readouterr()
    assert "deprecated" in cap.err
    assert "peakatail" in cap.err
    assert "deprecated" not in cap.out


def test_alias_actually_dispatches_the_same_cli(capsys, monkeypatch):
    """Behavioural, not structural: `ema --version` must still print the version."""
    monkeypatch.setattr(sys, "argv", ["ema", "--version"])
    with pytest.raises(SystemExit) as exc:
        ema_deprecated_main()
    assert exc.value.code == 0
    cap = capsys.readouterr()
    assert "version" in cap.out.lower() or any(ch.isdigit() for ch in cap.out)
    assert "deprecated" in cap.err
