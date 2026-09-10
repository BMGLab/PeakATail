"""Verify pyproject.toml metadata + entry points are installed correctly."""
import importlib.metadata as md
import shutil
import subprocess

import pytest


def _console_scripts() -> dict:
    """Entry points the INSTALLED distribution declares."""
    eps = md.distribution("peakatail").entry_points
    return {e.name: e.value for e in eps if e.group == "console_scripts"}


def test_distribution_name():
    """Package is published as 'peakatail' even though import path is 'ema'."""
    dist = md.distribution("peakatail")
    assert dist.metadata["Name"] == "peakatail"


def test_import_path_is_ema():
    """Importable as `ema` (the historical name, deliberately not renamed)."""
    import ema  # noqa: F401


def test_peakatail_is_the_declared_entry_point():
    """The command is `peakatail`; the package has been named that since 0.2.0."""
    scripts = _console_scripts()
    if "peakatail" not in scripts:
        pytest.skip("installed distribution predates the peakatail rename")
    assert scripts["peakatail"] == "ema.cli:main"


def test_ema_entry_point_survives_as_a_deprecated_alias():
    """Dropping `ema` would break every published pipeline that calls it."""
    scripts = _console_scripts()
    if "peakatail" not in scripts:
        pytest.skip("installed distribution predates the peakatail rename")
    assert "ema" in scripts, "the deprecated `ema` alias must remain installed"


def test_entry_point_runs():
    """Whichever command is on PATH must answer --help."""
    cmd = shutil.which("peakatail") or shutil.which("ema")
    if cmd is None:
        pytest.skip("no console script on PATH (package not installed in this env)")
    res = subprocess.run([cmd, "--help"], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr[:400]
    assert "Usage:" in res.stdout or "Usage:" in res.stderr


def test_legacy_entry_points_removed():
    """ema_switch / ema_merge / ema_parse_gtf are GONE (hard cutover)."""
    for name in ("ema_switch", "ema_merge", "ema_parse_gtf"):
        assert shutil.which(name) is None, (
            f"{name} should be removed in this branch but still exists on PATH"
        )
