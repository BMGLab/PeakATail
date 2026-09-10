"""Verify pyproject.toml metadata + entry points are installed correctly."""
import importlib.metadata as md
import shutil
import subprocess

import pytest


def _project_version() -> str:
    """The version this source tree declares.

    Read with a regex rather than tomllib so this file still works on the
    interpreter a contributor happens to have (tomllib is 3.11+, and the point
    of the guard below is to behave sanely in an environment that is NOT the
    supported one). Same single-line shape scripts/bump_version.py matches.
    """
    import re
    from pathlib import Path

    pyproject = Path(__file__).resolve().parents[1] / "pyproject.toml"
    match = re.search(r'(?m)^version = "([^"]+)"',
                      pyproject.read_text(encoding="utf-8"))
    assert match, "pyproject.toml has no top-level `version = \"...\"` line"
    return match.group(1)


def _distribution():
    """The installed `peakatail` distribution, or skip -- but never mask a rename.

    These tests assert on installed *metadata*, so they only mean anything
    where the package was installed, which is what CI does (`pip install -e
    '.[test]'`). From a bare source checkout there is no distribution and
    `md.distribution()` raises, which used to surface as three hard failures
    that said nothing about the code.

    Skipping on *any* PackageNotFoundError would be too blunt: "published
    under the wrong name" also raises it, and that is precisely the regression
    test_distribution_name exists to catch -- the skip would swallow it. So
    before skipping, ask which distribution actually provides the `ema` module.
    If one does and it carries THIS tree's version, the package is installed
    under the wrong name and that is a failure, not a skip. A stale install
    from an older version is neither, and is skipped.
    """
    try:
        return md.distribution("peakatail")
    except md.PackageNotFoundError:
        pass

    providers = md.packages_distributions().get("ema", [])
    for name in providers:
        try:
            installed = md.version(name)
        except md.PackageNotFoundError:  # pragma: no cover -- racy uninstall
            continue
        if installed == _project_version():
            pytest.fail(
                f"the `ema` package is installed as distribution {name!r} "
                f"at this tree's version ({installed}), not as 'peakatail'. "
                "The distribution name is what users `pip install` and what "
                "the PyPI project is called; renaming it silently breaks "
                "every install instruction and the bioconda recipe."
            )
    pytest.skip(
        "the 'peakatail' distribution is not installed in this interpreter "
        f"(providers of `ema`: {providers or 'none'}); run "
        "`pip install -e '.[test]'` to exercise the packaging metadata tests"
    )


def _console_scripts() -> dict:
    """Entry points the INSTALLED distribution declares."""
    eps = _distribution().entry_points
    return {e.name: e.value for e in eps if e.group == "console_scripts"}


def test_distribution_name():
    """Package is published as 'peakatail' even though import path is 'ema'."""
    dist = _distribution()
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
