"""Contract tests for `peakatail run`."""
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


def test_run_threads_flag_resets_resource_manager():
    """--threads must override the singleton ResourceManager.user_max_threads.

    Regression: the Click --threads flag landed in pipeline kwargs but was
    never bridged into ema.utils.get_resource_manager(); it only worked via
    the legacy argparse shim. Result: silently ignored on the new CLI.

    We check the wiring at the source-code level (no pipeline invocation)
    so this test cannot pollute module-level config state for subsequent
    tests. The matching e2e behaviour is exercised by the full atlas_full
    BAM run in the audit harness.
    """
    import inspect
    from ema.cli import run as run_mod

    # ``run_mod.run`` is the click.Command object; the actual callable is
    # exposed on ``.callback``.
    src = inspect.getsource(run_mod.run.callback)
    # Must reset_resource_manager + set _RM_INSTANCE = ResourceManager(user_max_threads=...)
    # gated on the user supplying --threads.
    assert "reset_resource_manager()" in src, (
        "ema.cli.run must reset the ResourceManager singleton when --threads is set"
    )
    assert 'kwargs["threads"]' in src or "kwargs.get(\"threads\")" in src, (
        "ema.cli.run must read kwargs['threads']"
    )
    assert "user_max_threads=kwargs[\"threads\"]" in src, (
        "ema.cli.run must construct ResourceManager(user_max_threads=kwargs['threads'])"
    )
    # And the drop set must include 'threads' so it is not forwarded as
    # an unrecognised pipeline kwarg.
    drop_src = inspect.getsource(run_mod._pipeline_kwargs)
    assert "\"threads\"" in drop_src, (
        "_pipeline_kwargs must drop 'threads' (already consumed by run.py itself)"
    )


def test_user_supplied_via_parameter_source(monkeypatch):
    """``_user_supplied_params`` must use Click's ParameterSource.COMMANDLINE.

    Regression: ``_apply_cli_overrides`` previously gated on
    ``v != DEFAULTS[key]`` -- which silently dropped ``--min-read 1500``
    when 1500 also happened to be the default, even though the user
    explicitly typed it.
    """
    import click
    from click.testing import CliRunner

    captured: dict = {}

    @click.command()
    @click.option("--min-read", "min_read", type=int, default=1500)
    @click.pass_context
    def fake(ctx, min_read):
        from ema.cli.run import _user_supplied_params
        captured["user_set"] = _user_supplied_params(ctx)
        captured["min_read"] = min_read

    runner = CliRunner()

    # Case 1: user passes --min-read 1500 (same as default).
    captured.clear()
    res = runner.invoke(fake, ["--min-read", "1500"])
    assert res.exit_code == 0
    assert "min_read" in captured["user_set"], (
        "user-typed value (even if it matches the default) must be reported "
        "as COMMANDLINE-sourced; got: " + repr(captured["user_set"])
    )

    # Case 2: user does NOT pass --min-read (Click fills the default).
    captured.clear()
    res = runner.invoke(fake, [])
    assert res.exit_code == 0
    assert "min_read" not in captured["user_set"], (
        "default-filled value must NOT be reported as COMMANDLINE-sourced; "
        "got: " + repr(captured["user_set"])
    )


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
