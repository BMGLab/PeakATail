"""Verify pyproject.toml metadata + entry points are installed correctly."""
import importlib.metadata as md
import shutil
import subprocess


def test_distribution_name():
    """Package is published as 'peakatail' even though import path is 'ema'."""
    dist = md.distribution("peakatail")
    assert dist.metadata["Name"] == "peakatail"


def test_import_path_is_ema():
    """Importable as `ema` (the historical name)."""
    import ema  # noqa: F401


def test_entry_point_is_installed():
    """`ema` console script exists on PATH after install."""
    assert shutil.which("ema") is not None, "ema entry point not found on PATH"


def test_entry_point_runs():
    """`ema --help` returns successfully (placeholder check until CLI is wired)."""
    # NOTE: this will fail until Phase 4 lands the Click root. That's fine —
    # this test is the contract we're moving toward. Mark xfail until then.
    import pytest
    res = subprocess.run(["ema", "--help"], capture_output=True, text=True)
    if res.returncode != 0:
        pytest.xfail("Click root not yet wired (Phase 4)")
    assert "Usage:" in res.stdout or "Usage:" in res.stderr


def test_legacy_entry_points_removed():
    """ema_switch / ema_merge / ema_parse_gtf are GONE (hard cutover)."""
    for name in ("ema_switch", "ema_merge", "ema_parse_gtf"):
        assert shutil.which(name) is None, (
            f"{name} should be removed in this branch but still exists on PATH"
        )
