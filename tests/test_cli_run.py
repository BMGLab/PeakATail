"""Contract tests for `ema run`."""
from click.testing import CliRunner

from ema.cli import main


def test_run_help_lists_all_groups():
    """--help must mention every flag group from spec §8.1."""
    runner = CliRunner()
    result = runner.invoke(main, ["run", "--help"])
    assert result.exit_code == 0
    out = result.output
    # Inputs
    assert "--config" in out
    assert "--gtf" in out
    assert "--atlas" in out
    # Concurrency
    assert "--threads" in out
    assert "--tiles" in out
    assert "--tile-size" in out
    # Peak
    assert "--peak-strategy" in out
    assert "--pas-gap" in out
    # Filters
    assert "--ip-filter" in out
    # Annotation
    assert "--max-gene-distance" in out
    # Clustering
    assert "--cluster-method" in out
    # Matching
    assert "--match-method" in out
    # Logging / common
    assert "--verbose" in out or "-v" in out


def test_run_list_strategies(monkeypatch):
    runner = CliRunner()
    result = runner.invoke(main, ["run", "--list-strategies"])
    assert result.exit_code == 0
    out = result.output
    assert "peak" in out.lower()
    assert "cluster" in out.lower()
    assert "match" in out.lower()


def test_run_invalid_strategy_rejected():
    runner = CliRunner()
    result = runner.invoke(main, ["run", "--peak-strategy", "BOGUS"])
    assert result.exit_code != 0
    assert "BOGUS" in result.output or "Invalid value" in result.output


def test_run_no_input_errors_clearly():
    """No --config and no --bam-dir → clear error, not stack trace."""
    runner = CliRunner()
    result = runner.invoke(main, ["run"])
    assert result.exit_code != 0
    assert "config" in result.output.lower() or "bam" in result.output.lower()
