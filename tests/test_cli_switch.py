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


def test_switch_length_isoform_agg_choices_match_strategy_vocab():
    """Click --isoform-agg vocabulary MUST match strategy code branches.

    Regression: when the CLI accepted "gene"/"isoform" but strategy code
    branched on "per_gene"/"per_isoform", `--isoform-agg gene` silently
    invoked the per_isoform code path with no warning.
    """
    runner = CliRunner()
    result = runner.invoke(main, ["switch", "length", "--help"])
    assert result.exit_code == 0
    out = result.output
    # Choice rendered in help — both per_* tokens must appear
    assert "per_gene" in out, f"per_gene missing from help:\n{out}"
    assert "per_isoform" in out, f"per_isoform missing from help:\n{out}"
    # Bad legacy values must be rejected loudly by Click
    bad = runner.invoke(main, [
        "switch", "length", "--h5ad", "/dev/null",
        "--isoform-agg", "gene",  # legacy short form is no longer accepted at CLI
    ])
    assert bad.exit_code != 0


def test_switch_length_isoform_collapse_choices():
    """--isoform-collapse must reject `weighted` (was an invalid default)."""
    runner = CliRunner()
    bad = runner.invoke(main, [
        "switch", "length", "--h5ad", "/dev/null",
        "--isoform-collapse", "weighted",
    ])
    assert bad.exit_code != 0
    # And the valid tokens are listed in --help
    h = runner.invoke(main, ["switch", "length", "--help"])
    for token in ("none", "mean", "majority"):
        assert token in h.output, f"expected {token} in help output"


def test_isoform_agg_per_gene_dispatches_to_per_gene_branch(monkeypatch):
    """`--isoform-agg per_gene` must invoke the per_gene branch in strategies.

    Direct regression for BLOCKER 1 — strategy code branches on the literal
    string. We capture the value passed to compute() to assert it.
    """
    captured: dict = {}

    class _StubStrategy:
        name = "stub"

        def compute(self, count_matrix, pas_isoform_map, aggregation,
                    isoform_collapse, **_):
            captured["aggregation"] = aggregation
            captured["isoform_collapse"] = isoform_collapse
            import pandas as pd
            return pd.DataFrame()

    # Patch the strategy registry so the call goes to our stub
    import ema.switch_test.runner as _runner_mod
    monkeypatch.setattr(_runner_mod, "get_pdui_strategy", lambda _name: _StubStrategy())

    # Stub anndata + count matrix path so we don't hit disk
    class _FakeAdata:
        pass

    monkeypatch.setattr(_runner_mod.ad, "read_h5ad", lambda _p: _FakeAdata())
    monkeypatch.setattr(
        _runner_mod, "build_count_dfs",
        lambda _adata, which="both": (__import__("pandas").DataFrame(), None, None, None),
    )

    # Use a path that exists so the loop body runs (we patched read_h5ad anyway)
    _runner_mod.run_length(
        h5ad_paths=["/dev/null"],
        gtf=None,
        output_dir="/tmp/_isoform_test_out",
        cluster_pairs=None,
        cluster_key="leiden",
        strategy="stub",
        isoform_agg="per_gene",
        isoform_collapse="none",
        threads=1,
    )

    assert captured.get("aggregation") == "per_gene", (
        f"strategy.compute() received aggregation={captured.get('aggregation')!r} "
        "but the user requested per_gene — vocabulary drift would silently "
        "invoke the wrong branch"
    )
