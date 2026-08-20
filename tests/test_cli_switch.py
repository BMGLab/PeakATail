"""Contract tests for `ema switch *`."""
import pytest
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


def test_isoform_agg_per_gene_dispatches_to_per_gene_branch(monkeypatch, tmp_path):
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
        def __init__(self):
            import pandas as pd
            self.var = pd.DataFrame({"gene_id": ["G1", "G1"]}, index=["1", "2"])
            self.obs = pd.DataFrame(index=[])

    monkeypatch.setattr(_runner_mod.ad, "read_h5ad", lambda _p: _FakeAdata())
    monkeypatch.setattr(
        _runner_mod, "build_count_dfs",
        lambda _adata, which="both", **_kw: (
            __import__("pandas").DataFrame(), None, None, None),
    )

    # run_length resolves pasbed.bed next to the h5ad and REFUSES to run
    # without it (ranking PAS by input order inverts every minus-strand gene).
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "pasbed.bed").write_text(
        "chr1\t100\t101\t1\t0\t+\n"
        "chr1\t500\t501\t2\t0\t+\n"
    )

    # Use a path whose directory exists so the loop body runs (we patched
    # read_h5ad anyway)
    _runner_mod.run_length(
        h5ad_paths=[str(run_dir / "clusters.h5ad")],
        gtf=None,
        output_dir=str(tmp_path / "out"),
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


def test_per_gene_pas_isoform_map_rank_from_pasbed(monkeypatch, tmp_path):
    """The per_gene ``pas_isoform_map`` synthesis must rank each gene's PAS
    by real (strand-aware) genomic position, not the old hardcoded
    ``rank=1, total=1`` sentinel that made every PAS of every gene
    indistinguishable from a true single-PAS gene.

    GENE_PLUS (+ strand): PAS 10 (start=100), PAS 11 (start=500)
        -> rank 1 (proximal, low coord), rank 2 (distal, high coord).
    GENE_MINUS (- strand): PAS 20 (start=500), PAS 21 (start=100)
        -> rank 1 (proximal, HIGH coord), rank 2 (distal, LOW coord).
    """
    import pandas as pd

    captured: dict = {}

    class _StubStrategy:
        name = "stub"

        def compute(self, count_matrix, pas_isoform_map, aggregation,
                    isoform_collapse, **_):
            captured["pas_isoform_map"] = pas_isoform_map
            return pd.DataFrame()

    import ema.switch_test.runner as _runner_mod
    monkeypatch.setattr(_runner_mod, "get_pdui_strategy", lambda _name: _StubStrategy())

    class _FakeAdata:
        def __init__(self):
            self.var = pd.DataFrame(
                {"gene_id": ["GENE_PLUS", "GENE_PLUS", "GENE_MINUS", "GENE_MINUS"]},
                index=["10", "11", "20", "21"],
            )
            self.obs = pd.DataFrame(index=[])

    monkeypatch.setattr(_runner_mod.ad, "read_h5ad", lambda _p: _FakeAdata())
    monkeypatch.setattr(
        _runner_mod, "build_count_dfs",
        lambda _adata, which="both", **_kw: (pd.DataFrame(), None, None, None),
    )

    # Real pasbed.bed sibling to the (stubbed) h5ad path — this is what
    # `run_length` walks up from to find genomic coordinates.
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    pasbed = run_dir / "pasbed.bed"
    pasbed.write_text(
        "chr1\t100\t101\t10\t0\t+\n"
        "chr1\t500\t501\t11\t0\t+\n"
        "chr1\t500\t501\t20\t0\t-\n"
        "chr1\t100\t101\t21\t0\t-\n"
    )
    h5ad_path = run_dir / "clusters.h5ad"

    _runner_mod.run_length(
        h5ad_paths=[str(h5ad_path)],
        gtf=None,
        output_dir=str(tmp_path / "out"),
        cluster_pairs=None,
        cluster_key="leiden",
        strategy="stub",
        isoform_agg="per_gene",
        isoform_collapse="none",
        threads=1,
    )

    pas_map = captured["pas_isoform_map"]
    assert pas_map[10][0][3] == 1 and pas_map[10][0][4] == 2   # GENE_PLUS proximal
    assert pas_map[11][0][3] == 2 and pas_map[11][0][4] == 2   # GENE_PLUS distal
    assert pas_map[20][0][3] == 1 and pas_map[20][0][4] == 2   # GENE_MINUS proximal (high coord)
    assert pas_map[21][0][3] == 2 and pas_map[21][0][4] == 2   # GENE_MINUS distal (low coord)


def test_switch_length_exposes_pasbed_and_counts_flags():
    """0b/0c: `switch length` must be able to be TOLD where the truth is.

    Without ``--pasbed`` the only way to supply PAS strand was the directory
    walk-up, and when that missed, ranking fell back to input order (a silent
    plus-strand convention).  Without ``--counts-layer`` there was no way to
    point PDUI at raw counts once clustering had overwritten ``.X``.
    """
    runner = CliRunner()
    result = runner.invoke(main, ["switch", "length", "--help"])
    assert result.exit_code == 0
    for flag in ("--pasbed", "--counts-layer", "--allow-non-count-matrix"):
        assert flag in result.output, f"{flag} missing from `switch length --help`"


def test_switch_diff_exposes_counts_flags():
    runner = CliRunner()
    result = runner.invoke(main, ["switch", "diff", "--help"])
    assert result.exit_code == 0
    for flag in ("--counts-layer", "--allow-non-count-matrix"):
        assert flag in result.output, f"{flag} missing from `switch diff --help`"


def test_switch_length_cli_forwards_new_flags_to_run_length(monkeypatch, tmp_path):
    """The flags must actually reach ``run_length`` (not just parse)."""
    import ema.cli.switch_length as _cli_mod

    captured: dict = {}

    def _fake_run_length(**kw):
        captured.update(kw)
        return None, None

    monkeypatch.setattr(
        "ema.switch_test.runner.run_length", _fake_run_length, raising=True
    )
    monkeypatch.setattr(
        _cli_mod, "resolve_subcommand_output_dir",
        lambda *a, **k: tmp_path / "out", raising=True,
    )
    h5 = tmp_path / "clusters.h5ad"
    h5.write_bytes(b"")
    pasbed = tmp_path / "pasbed.bed"
    pasbed.write_text("chr1\t1\t2\t1\t0\t+\n")

    runner = CliRunner()
    result = runner.invoke(main, [
        "switch", "length", "--h5ad", str(h5), "--pasbed", str(pasbed),
        "--counts-layer", "counts", "--allow-non-count-matrix",
        "--no-plots", "--no-log-file",
    ])
    assert result.exit_code == 0, result.output
    assert captured["pasbed"] == str(pasbed)
    assert captured["counts_layer"] == "counts"
    assert captured["allow_non_count_matrix"] is True
