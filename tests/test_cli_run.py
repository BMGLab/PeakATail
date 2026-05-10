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


def test_run_threads_flag_resets_resource_manager(tmp_path, monkeypatch):
    """--threads must override the singleton ResourceManager.user_max_threads.

    Regression: the Click --threads flag landed in pipeline kwargs but was
    never bridged into ema.utils.get_resource_manager(); it only worked via
    the legacy argparse shim. Result: silently ignored on the new CLI.
    """
    # Force a clean singleton, supply --threads, capture state before pipeline
    # starts by aborting via a YAML loader patch.
    from ema.utils import reset_resource_manager, get_resource_manager
    reset_resource_manager()

    yaml_path = tmp_path / "minimal.yaml"
    yaml_path.write_text(
        "datasets:\n"
        "  - id: x\n"
        "    merge_strategy: none\n"
        "    bams:\n"
        "      - /nonexistent.bam\n"
    )

    captured: dict = {}

    def _abort_setup_logging(*args, **kwargs):
        # Capture ResourceManager state at the precise moment that
        # the --threads wiring (which happens BEFORE setup_logging) has run.
        captured["rm_user_max"] = get_resource_manager().user_max_threads
        raise SystemExit(123)

    monkeypatch.setattr("ema.logging_config.setup_logging", _abort_setup_logging)
    runner = CliRunner()
    result = runner.invoke(
        main, ["run", "--config", str(yaml_path), "--threads", "7"],
        standalone_mode=False,
    )
    # Either Click caught the SystemExit or the test marker propagated.
    assert captured.get("rm_user_max") == 7, (
        f"--threads not bridged into ResourceManager singleton; "
        f"got {captured.get('rm_user_max')!r}"
    )

    reset_resource_manager()


def test_single_sample_preprocessing_passes_filter_kwargs():
    """main._run_pipeline_body() must forward filter_config.min_cells/min_genes
    into matrixfilter.preprocessing() on the single-sample path.

    Regression: preprocessing()'s default kwargs are evaluated at function-def
    time so they snapshot filter_config from import. Without the explicit
    forward, YAML overrides applied to filter_config were silently ignored.
    """
    import inspect
    import ema.main as main_mod

    body = inspect.getsource(main_mod._run_pipeline_body)
    assert "min_cells=filter_config.min_cells" in body, (
        "single-sample preprocessing() call must explicitly pass min_cells "
        "from filter_config or YAML overrides will be silently ignored"
    )
    assert "min_genes=filter_config.min_genes" in body, (
        "single-sample preprocessing() call must explicitly pass min_genes "
        "from filter_config or YAML overrides will be silently ignored"
    )
