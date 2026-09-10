"""Click root group contract."""
from click.testing import CliRunner

from ema.cli import main


def test_help_runs():
    runner = CliRunner()
    result = runner.invoke(main, ["--help"])
    assert result.exit_code == 0
    assert "Usage:" in result.output
    assert "Commands:" in result.output


def test_version_flag():
    runner = CliRunner()
    result = runner.invoke(main, ["--version"])
    assert result.exit_code == 0
    assert "peakatail" in result.output.lower() or "ema" in result.output.lower()


def test_unknown_command_exits_nonzero():
    runner = CliRunner()
    result = runner.invoke(main, ["nonsense"])
    assert result.exit_code != 0


def test_bare_invocation_invokes_wizard(monkeypatch):
    """Bare `ema` (no subcommand) launches the wizard."""
    called = {}

    def fake_wizard():
        called["yes"] = True
        return 0

    monkeypatch.setattr("ema.cli.wizard.run", fake_wizard)
    runner = CliRunner()
    result = runner.invoke(main, [])
    assert called.get("yes"), "wizard.run() should have been invoked"


from ema.cli.common import parse_log_overrides, resolve_output_dir


def test_parse_log_overrides_empty():
    assert parse_log_overrides(None) == {}
    assert parse_log_overrides("") == {}


def test_parse_log_overrides_single_level():
    assert parse_log_overrides("DEBUG") == {"ema": "DEBUG"}


def test_parse_log_overrides_per_logger():
    assert parse_log_overrides("ema.countmatrix=WARNING") == {"ema.countmatrix": "WARNING"}


def test_parse_log_overrides_combined():
    assert parse_log_overrides("DEBUG,ema.cm=WARNING") == {"ema": "DEBUG", "ema.cm": "WARNING"}


def test_resolve_output_dir_appends_timestamp(tmp_path):
    base = tmp_path / "out"
    result = resolve_output_dir(str(base))
    # Format: out_2026-05-10_181523
    import re
    assert re.match(r".*out_\d{4}-\d{2}-\d{2}_\d{6}$", str(result))
