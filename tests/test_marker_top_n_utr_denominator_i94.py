"""Issue #94 regression: ``--marker-top-n`` must not narrow the UTR background.

PR #105 repaired the ``per_gene`` scope: the marker pre-selection decides only
WHICH PAS are reported, while the within-gene Fisher denominator is computed
from the FULL count matrix (``full_count_matrix``).  The UTR-scoped modes did
not get that treatment -- ``_build_diff_isoform_groups`` was handed the
already marker-restricted matrix, so under ``--isoform-agg within_utr`` a
UTR's background silently became "the marker-selected PAS of this UTR", and
under ``between_utr`` each UTR column summed only its selected member PAS.
Opting into a *speed* filter therefore changed the statistic itself.

These tests pin the same contract the per_gene fix pins, one scope down:
a marker-restricted run reproduces the unrestricted run's denominators and
p-values exactly for every row it still reports, and a ``--marker-top-n 0``
run is untouched.

Mechanics mirror ``tests/test_run_diff_isoform.py``: the real ``fisher``
strategy and the real ``bedtools intersect`` inside ``map_pas_to_isoforms``
run; only GTF parsing and the (label-double-dipping, hence irrelevant here)
marker ranking are stubbed, so the marker set is a fixed list.

Fixture layout
---------------
GENE_X, two 3'UTRs:
    TRANSCRIPT_X1 (chr1:1000-1200): PAS 1, 2, 3, 4
    TRANSCRIPT_X2 (chr1:2000-2200): PAS 5, 6
Markers: PAS 1, 2, 5 -- so PAS 3, 4 and 6 carry counts that belong in the
background of a group whose reported rows are all markers.  Counts are
deliberately SPARSE so the default ``count_mode="cells"`` denominator
("cells where ANY member is detected") does not saturate.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

import ema.switch_test.runner as runner_mod

ISOFORM_UTRS = {
    "GENE_X": {
        "TRANSCRIPT_X1": [("chr1", 1000, 1200, "+", 0)],
        "TRANSCRIPT_X2": [("chr1", 2000, 2200, "+", 0)],
    },
}

# (chrom, start, end, pas_id, score, strand)
PASBED_RECORDS = [
    ("chr1", 1050, 1052, 1, 0, "+"),   # TRANSCRIPT_X1
    ("chr1", 1100, 1102, 2, 0, "+"),   # TRANSCRIPT_X1
    ("chr1", 1150, 1152, 3, 0, "+"),   # TRANSCRIPT_X1
    ("chr1", 1180, 1182, 4, 0, "+"),   # TRANSCRIPT_X1
    ("chr1", 2050, 2052, 5, 0, "+"),   # TRANSCRIPT_X2
    ("chr1", 2150, 2152, 6, 0, "+"),   # TRANSCRIPT_X2
]

PAS_IDS = ["1", "2", "3", "4", "5", "6"]
UTR1_PAS = ["1", "2", "3", "4"]
UTR2_PAS = ["5", "6"]
MARKERS = ["1", "2", "5"]
NON_MARKERS = ["3", "4", "6"]

UTR1 = "GENE_X::TRANSCRIPT_X1"
UTR2 = "GENE_X::TRANSCRIPT_X2"

N_CELLS = 120
COUNT_SEED = 4242

_GENE_COLS = ["n_reads_gene_cluster1", "n_reads_gene_cluster2"]


def _counts() -> np.ndarray:
    rng = np.random.RandomState(COUNT_SEED)
    X = rng.poisson(0.25, size=(N_CELLS, len(PAS_IDS))).astype(float)
    # A little real signal at PAS 1 so the p-values are not a wall of 1.0.
    X[: N_CELLS // 2, 0] += rng.poisson(0.6, size=N_CELLS // 2)
    return X


def _make_adata() -> ad.AnnData:
    var = pd.DataFrame({"gene_id": ["GENE_X"] * len(PAS_IDS)}, index=PAS_IDS)
    obs = pd.DataFrame(
        {"leiden": ["A"] * (N_CELLS // 2) + ["B"] * (N_CELLS // 2)},
        index=[f"c{i}" for i in range(N_CELLS)],
    )
    return ad.AnnData(X=_counts(), obs=obs, var=var)


def _read_sums(pas_ids: list[str]) -> tuple[int, int]:
    """``(cluster A reads, cluster B reads)`` summed over ``pas_ids``."""
    X = _counts()
    cols = [PAS_IDS.index(p) for p in pas_ids]
    half = N_CELLS // 2
    return (
        int(X[:half][:, cols].sum()),
        int(X[half:][:, cols].sum()),
    )


@pytest.fixture(scope="module")
def run_env(tmp_path_factory):
    """``run(isoform_agg, marker_top_n) -> pair DataFrame`` on one input."""
    tmp_path = tmp_path_factory.mktemp("i94_utr_denominator")
    (tmp_path / "pasbed.bed").write_text(
        "".join(f"{c}\t{s}\t{e}\t{p}\t{sc}\t{st}\n" for c, s, e, p, sc, st in PASBED_RECORDS)
    )
    gtf_path = tmp_path / "genome.gtf"
    gtf_path.write_text("")  # only needs to exist

    import ema.annotate.gtf2isoform_utr as gtf2isoform_mod

    real_read_h5ad = runner_mod.ad.read_h5ad
    real_parse = gtf2isoform_mod.parse_isoform_utrs
    real_select = runner_mod.select_marker_pas
    runner_mod.ad.read_h5ad = lambda _p: _make_adata()
    gtf2isoform_mod.parse_isoform_utrs = lambda *a, **kw: ISOFORM_UTRS
    # The marker RANKING is the (unfixable-by-code) label double-dip half of
    # #94 and is not what is under test here; pin the selection so the
    # denominator is the only moving part.
    # ``select_marker_pas`` returns INT pas ids (and ``build_count_dfs``
    # keys the count frame the same way), so the stub must too.
    runner_mod.select_marker_pas = lambda *a, **kw: [int(p) for p in MARKERS]

    def run(isoform_agg: str, marker_top_n: int) -> pd.DataFrame:
        pair_results = runner_mod.run_diff(
            h5ad_paths=[str(tmp_path / "clustered.h5ad")],
            pasbed=None,
            gtf=str(gtf_path),
            output_dir=str(tmp_path / f"{isoform_agg}_top{marker_top_n}"),
            cluster_pairs=None,
            cluster_key="leiden",
            marker_top_n=marker_top_n,
            marker_method="wilcoxon",
            strategy="fisher",
            fdr=0.05,
            threads=1,
            per_worker_mb=300,
            min_cells_per_group=5,
            isoform_agg=isoform_agg,
            utr_unmatched="drop",
        )
        assert list(pair_results) == [("A", "B")]
        return pair_results[("A", "B")]

    try:
        yield run
    finally:
        runner_mod.ad.read_h5ad = real_read_h5ad
        gtf2isoform_mod.parse_isoform_utrs = real_parse
        runner_mod.select_marker_pas = real_select


@pytest.fixture(scope="module")
def within_utr_frames(run_env):
    """``(marker_top_n=0, marker_top_n=3)`` within_utr results, keyed by row."""
    off = run_env("within_utr", 0)
    on = run_env("within_utr", len(MARKERS))
    return off, on


@pytest.fixture(scope="module")
def between_utr_frames(run_env):
    off = run_env("between_utr", 0)
    on = run_env("between_utr", len(MARKERS))
    return off, on


def _row(df: pd.DataFrame, row_id: str) -> pd.Series:
    rows = df[df["pas_id"] == row_id]
    assert len(rows) == 1, f"{row_id}: expected 1 row, got {len(rows)}"
    return rows.iloc[0]


# ---------------------------------------------------------------------------
# within_utr
# ---------------------------------------------------------------------------


class TestWithinUtrDenominator:
    def test_marker_restriction_only_gates_reported_rows(self, within_utr_frames):
        off, on = within_utr_frames
        assert set(off["pas_id"]) == set(PAS_IDS)
        assert set(on["pas_id"]) == set(MARKERS)

    def test_background_is_the_whole_utr_not_the_selected_pas(self, within_utr_frames):
        """PAS 1's UTR total must count PAS 3 and 4 too, marker run or not."""
        _, on = within_utr_frames
        expected_c1, expected_c2 = _read_sums(UTR1_PAS)
        row = _row(on, "1")
        assert int(row["n_reads_gene_cluster1"]) == expected_c1
        assert int(row["n_reads_gene_cluster2"]) == expected_c2

    def test_marker_run_reproduces_unrestricted_pvalues(self, within_utr_frames):
        off, on = within_utr_frames
        for pas_id in MARKERS:
            off_row, on_row = _row(off, pas_id), _row(on, pas_id)
            for col in _GENE_COLS:
                assert int(on_row[col]) == int(off_row[col]), (
                    f"{col} for PAS {pas_id} is {on_row[col]} under "
                    f"--marker-top-n but {off_row[col]} without it -- the "
                    f"UTR background must come from ALL member PAS (issue #94)"
                )
            assert on_row["pvalue"] == pytest.approx(off_row["pvalue"], rel=1e-12), (
                f"p-value for PAS {pas_id} moved merely because non-marker PAS "
                f"of its UTR were not selected for testing (issue #94)"
            )
            for col in ("delta_proportion", "log2fc", "odds_ratio"):
                assert on_row[col] == pytest.approx(off_row[col], rel=1e-12)

    def test_fixture_would_expose_a_narrowed_background(self):
        """Guard: the non-marker PAS really do carry counts.

        Without this the assertions above could pass vacuously on a fixture
        where "the selected PAS of this UTR" happens to equal "this UTR".
        """
        full_c1, full_c2 = _read_sums(UTR1_PAS)
        sel_c1, sel_c2 = _read_sums([p for p in UTR1_PAS if p in MARKERS])
        assert sel_c1 < full_c1 and sel_c2 < full_c2


# ---------------------------------------------------------------------------
# between_utr
# ---------------------------------------------------------------------------


class TestBetweenUtrDenominator:
    def test_utr_columns_sum_all_member_pas(self, between_utr_frames):
        """A UTR row is reported when a member PAS was selected, but its
        counts stay the sum of ALL its members."""
        _, on = between_utr_frames
        assert set(on["pas_id"]) == {UTR1, UTR2}
        expected_c1, expected_c2 = _read_sums(PAS_IDS)  # gene = both UTRs
        for row_id in (UTR1, UTR2):
            row = _row(on, row_id)
            assert int(row["n_reads_gene_cluster1"]) == expected_c1
            assert int(row["n_reads_gene_cluster2"]) == expected_c2

    def test_marker_run_reproduces_unrestricted_pvalues(self, between_utr_frames):
        off, on = between_utr_frames
        assert set(off["pas_id"]) == {UTR1, UTR2}
        for row_id in (UTR1, UTR2):
            off_row, on_row = _row(off, row_id), _row(on, row_id)
            for col in _GENE_COLS + ["n_reads_pas_cluster1", "n_reads_pas_cluster2"]:
                assert int(on_row[col]) == int(off_row[col]), (
                    f"{col} for {row_id} is {on_row[col]} under "
                    f"--marker-top-n but {off_row[col]} without it (issue #94)"
                )
            assert on_row["pvalue"] == pytest.approx(off_row["pvalue"], rel=1e-12)


# ---------------------------------------------------------------------------
# --marker-top-n 0 is untouched
# ---------------------------------------------------------------------------


class TestMarkerTopN0Unchanged:
    def test_within_utr_top_n_0_reports_every_pas_at_utr_scope(self, within_utr_frames):
        off, _ = within_utr_frames
        assert set(off["pas_id"]) == set(PAS_IDS)
        for pas_id in PAS_IDS:
            row = _row(off, pas_id)
            members = UTR1_PAS if pas_id in UTR1_PAS else UTR2_PAS
            expected_c1, expected_c2 = _read_sums(members)
            assert row["diff_group_id"] == (UTR1 if pas_id in UTR1_PAS else UTR2)
            assert int(row["n_reads_gene_cluster1"]) == expected_c1
            assert int(row["n_reads_gene_cluster2"]) == expected_c2

    def test_between_utr_top_n_0_reports_both_utrs(self, between_utr_frames):
        off, _ = between_utr_frames
        assert set(off["pas_id"]) == {UTR1, UTR2}
        expected_c1, expected_c2 = _read_sums(PAS_IDS)
        for row_id in (UTR1, UTR2):
            row = _row(off, row_id)
            assert row["diff_group_id"] == "GENE_X"
            assert int(row["n_reads_gene_cluster1"]) == expected_c1
            assert int(row["n_reads_gene_cluster2"]) == expected_c2


# ---------------------------------------------------------------------------
# the obsolete warning
# ---------------------------------------------------------------------------


def test_no_longer_warns_that_the_utr_background_is_narrowed(run_env, caplog):
    """The warning #105 left behind described a defect that is now fixed."""
    with caplog.at_level(logging.WARNING, logger="ema.switch_test.runner"):
        run_env("within_utr", len(MARKERS))
    text = "\n".join(
        r.getMessage() for r in caplog.records if r.levelno >= logging.WARNING
    )
    assert "narrows the within-group denominator" not in text
