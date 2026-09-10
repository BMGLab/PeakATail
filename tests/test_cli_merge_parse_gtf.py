"""Contract tests for `peakatail merge` and `peakatail parse-gtf`."""
from click.testing import CliRunner

from ema.cli import main


def test_merge_help():
    runner = CliRunner()
    result = runner.invoke(main, ["merge", "--help"])
    assert result.exit_code == 0
    assert "--bam-files" in result.output
    assert "--output" in result.output


def test_merge_no_input_errors():
    runner = CliRunner()
    result = runner.invoke(main, ["merge"])
    assert result.exit_code != 0


def test_parse_gtf_help():
    runner = CliRunner()
    result = runner.invoke(main, ["parse-gtf", "--help"])
    assert result.exit_code == 0
    assert "--gtf" in result.output
    assert "--cache-dir" in result.output


def test_parse_gtf_no_gtf_errors():
    runner = CliRunner()
    result = runner.invoke(main, ["parse-gtf"])
    assert result.exit_code != 0
