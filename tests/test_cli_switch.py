"""Contract tests for `ema switch *`."""
from click.testing import CliRunner

from ema.cli import main


def test_switch_help_lists_subcommands():
    runner = CliRunner()
    result = runner.invoke(main, ["switch", "--help"])
    assert result.exit_code == 0
    assert "diff" in result.output
    assert "length" in result.output
    assert "match" in result.output


def test_switch_diff_help_lists_strategy_choices():
    runner = CliRunner()
    result = runner.invoke(main, ["switch", "diff", "--help"])
    assert result.exit_code == 0
    out = result.output
    assert "--strategy" in out
    assert "--h5ad" in out
    assert "--fdr" in out
    # at least one of the strategies named
    assert "fisher" in out or "nb_pairwise" in out


def test_switch_length_help():
    runner = CliRunner()
    result = runner.invoke(main, ["switch", "length", "--help"])
    assert result.exit_code == 0
    assert "--strategy" in result.output
    assert "--isoform-agg" in result.output


def test_switch_match_help():
    runner = CliRunner()
    result = runner.invoke(main, ["switch", "match", "--help"])
    assert result.exit_code == 0
    assert "--strategy" in result.output
    assert "--n-top-markers" in result.output


def test_switch_diff_invalid_strategy_rejected():
    runner = CliRunner()
    result = runner.invoke(main, ["switch", "diff", "--h5ad", "x.h5ad", "--strategy", "BOGUS"])
    assert result.exit_code != 0
