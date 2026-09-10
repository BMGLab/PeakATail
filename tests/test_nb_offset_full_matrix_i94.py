"""End-to-end: a pre-selection must not move NB p-values (issue #94).

The unit tests in ``test_label_independent_prefilter_i94.py`` call the strategy
directly, so they pin the strategy's own behaviour but say nothing about
whether ``run_diff`` actually FORWARDS the unrestricted matrix. Deleting the
``full_count_matrix=diff_df_denom`` argument in ``runner.py`` left every one of
them green while real runs silently regressed -- so this file drives the whole
``run_diff`` path instead.

The strategy is ``nb_pairwise`` and the scope is ``per_gene``, the combination
whose offset the restriction used to change.
"""
from __future__ import annotations

import pytest

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import ema.switch_test.runner as runner_mod  # noqa: E402
from test_marker_top_n_utr_denominator_i94 import (  # noqa: E402
    ISOFORM_UTRS, MARKERS, PASBED_RECORDS, _make_adata,
)


@pytest.fixture(scope="module")
def nb_run(tmp_path_factory):
    """``run(marker_top_n) -> pair DataFrame`` through the real run_diff."""
    tmp_path = tmp_path_factory.mktemp("i94_nb_offset")
    (tmp_path / "pasbed.bed").write_text(
        "".join(f"{c}\t{s}\t{e}\t{p}\t{sc}\t{st}\n"
                for c, s, e, p, sc, st in PASBED_RECORDS)
    )
    gtf_path = tmp_path / "genome.gtf"
    gtf_path.write_text("")

    import ema.annotate.gtf2isoform_utr as gtf2isoform_mod

    real_read_h5ad = runner_mod.ad.read_h5ad
    real_parse = gtf2isoform_mod.parse_isoform_utrs
    real_select = runner_mod.select_marker_pas
    runner_mod.ad.read_h5ad = lambda _p: _make_adata()
    gtf2isoform_mod.parse_isoform_utrs = lambda *a, **kw: ISOFORM_UTRS
    # Pin the selection: the marker RANKING is the label double-dip half of
    # #94 and is not under test. The offset is the only moving part.
    runner_mod.select_marker_pas = lambda *a, **kw: [int(p) for p in MARKERS]

    def run(marker_top_n: int):
        pair_results = runner_mod.run_diff(
            h5ad_paths=[str(tmp_path / "clustered.h5ad")],
            pasbed=None,
            gtf=str(gtf_path),
            output_dir=str(tmp_path / f"nb_top{marker_top_n}"),
            cluster_pairs=None,
            cluster_key="leiden",
            marker_top_n=marker_top_n,
            marker_method="wilcoxon",
            strategy="nb_pairwise",
            fdr=0.05,
            threads=1,
            per_worker_mb=300,
            min_cells_per_group=2,
            isoform_agg="per_gene",
            utr_unmatched="drop",
        )
        assert list(pair_results) == [("A", "B")], list(pair_results)
        return pair_results[("A", "B")]

    try:
        yield run
    finally:
        runner_mod.ad.read_h5ad = real_read_h5ad
        gtf2isoform_mod.parse_isoform_utrs = real_parse
        runner_mod.select_marker_pas = real_select


def test_marker_top_n_does_not_move_nb_pvalues_end_to_end(nb_run):
    """run_diff must hand nb_pairwise the unrestricted matrix for its offset."""
    off = nb_run(0).set_index("pas_id")
    on = nb_run(len(MARKERS)).set_index("pas_id")

    shared = [p for p in on.index if p in off.index]
    assert shared, f"no shared rows: off={list(off.index)} on={list(on.index)}"

    mismatched = {
        p: (off.loc[p, "pvalue"], on.loc[p, "pvalue"])
        for p in shared
        if off.loc[p, "pvalue"] != on.loc[p, "pvalue"]
    }
    assert not mismatched, (
        "--marker-top-n changed nb_pairwise p-values through run_diff: "
        f"{mismatched}. run_diff must pass full_count_matrix=diff_df_denom so "
        "the per-cell library-size offset stays the cell's sequencing depth "
        "(issue #94)."
    )


@pytest.fixture(scope="module")
def nb_multi_run(tmp_path_factory):
    """Same, for the OMNIBUS branch.

    ``nb_pairwise`` is not multi-condition, so the test above never reaches
    ``runner.py``'s omnibus call -- which was a separate, independently
    un-forwarded site. Deleting its ``full_count_matrix`` left the pairwise
    test green.
    """
    tmp_path = tmp_path_factory.mktemp("i94_nb_multi_offset")
    (tmp_path / "pasbed.bed").write_text(
        "".join(f"{c}\t{s}\t{e}\t{p}\t{sc}\t{st}\n"
                for c, s, e, p, sc, st in PASBED_RECORDS)
    )
    gtf_path = tmp_path / "genome.gtf"
    gtf_path.write_text("")

    import ema.annotate.gtf2isoform_utr as gtf2isoform_mod

    real_read_h5ad = runner_mod.ad.read_h5ad
    real_parse = gtf2isoform_mod.parse_isoform_utrs
    real_select = runner_mod.select_marker_pas
    runner_mod.ad.read_h5ad = lambda _p: _make_adata()
    gtf2isoform_mod.parse_isoform_utrs = lambda *a, **kw: ISOFORM_UTRS
    runner_mod.select_marker_pas = lambda *a, **kw: [int(p) for p in MARKERS]

    def run(marker_top_n: int):
        pair_results = runner_mod.run_diff(
            h5ad_paths=[str(tmp_path / "clustered.h5ad")],
            pasbed=None,
            gtf=str(gtf_path),
            output_dir=str(tmp_path / f"nbmulti_top{marker_top_n}"),
            cluster_pairs=None,
            cluster_key="leiden",
            marker_top_n=marker_top_n,
            marker_method="wilcoxon",
            strategy="nb_multi",
            fdr=0.05,
            threads=1,
            per_worker_mb=300,
            min_cells_per_group=2,
            isoform_agg="per_gene",
            utr_unmatched="drop",
        )
        key = ("omnibus", "")
        assert key in pair_results, list(pair_results)
        return pair_results[key]

    try:
        yield run
    finally:
        runner_mod.ad.read_h5ad = real_read_h5ad
        gtf2isoform_mod.parse_isoform_utrs = real_parse
        runner_mod.select_marker_pas = real_select


def test_marker_top_n_does_not_move_nb_multi_omnibus_pvalues(nb_multi_run):
    off = nb_multi_run(0)
    on = nb_multi_run(len(MARKERS))

    shared = [p for p in on.index if p in off.index]
    assert shared, f"no shared rows: off={list(off.index)} on={list(on.index)}"

    mismatched = {
        p: (off.loc[p, "pvalue"], on.loc[p, "pvalue"])
        for p in shared
        if off.loc[p, "pvalue"] != on.loc[p, "pvalue"]
    }
    assert not mismatched, (
        "--marker-top-n changed nb_multi omnibus p-values: "
        f"{mismatched}. The omnibus call in run_diff must also pass "
        "full_count_matrix=diff_df_denom (issue #94)."
    )
