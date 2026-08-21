"""Tests for `ema switch diff`'s ``--isoform-agg`` scoping.

Covers the 3-way differential-test scope added on top of the plain
``per_gene`` background:

- ``per_gene`` (default): unchanged/byte-identical -- each PAS is tested
  against the rest of its GENE (legacy behaviour).
- ``within_utr`` (a.k.a. legacy alias ``per_isoform``): each PAS is tested
  against the OTHER PAS SHARING ITS 3'UTR isoform (tandem-UTR APA), not the
  whole gene.  ``--utr-unmatched`` controls PAS with no annotated UTR
  overlap: ``"gene"`` (default) keeps them via a gene-level fallback
  bucket; ``"drop"`` omits them.  A PAS overlapping >=1 UTR is tested once
  per UTR (one output row per UTR).
- ``between_utr``: PAS are collapsed to 3'UTR-level counts and the UNIT of
  the test becomes the 3'UTR itself -- "does 3'UTR preference differ
  between groups" (genes with >= 2 UTRs only).

All three scopes flow through the single shared choke point
``ema.switch_test.runner._run_grouped_diff`` (fed by
``_build_diff_isoform_groups``), so this is exercised end-to-end with the
REAL ``fisher`` strategy + the REAL ``bedtools intersect`` call inside
``map_pas_to_isoforms`` -- only GTF parsing (``parse_isoform_utrs``) is
stubbed, mirroring ``tests/test_run_length_utr_unmatched.py``.

Fixture layout
---------------
GENE_X: two clean, non-overlapping 3'UTRs, two PAS each --
    TRANSCRIPT_X1 (chr1:1000-1200): PAS 1, PAS 2
    TRANSCRIPT_X2 (chr1:2000-2200): PAS 3, PAS 4
    PAS 5: outside both UTRs (orphan; utr_unmatched fallback target)

GENE_Z: two transcripts sharing the SAME first UTR exon, so PAS 6 and PAS 7
both intersect BOTH transcripts' UTRs -- multi-UTR PAS:
    TRANSCRIPT_Z1 (chr3:3000-3500)
    TRANSCRIPT_Z2 (chr3:3000-3500) + (chr3:4000-4500)
    PAS 6, PAS 7: inside chr3:3000-3500 -> map to BOTH T_Z1 and T_Z2.

GENE_W: a single 3'UTR with one PAS -- proves single-UTR genes/UTRs never
get tested (nothing to contrast).
    TRANSCRIPT_W1 (chr4:5000-5200): PAS 8

Two clusters "A" (cellA1, cellA2) and "B" (cellB1, cellB2) with counts
chosen so a UTR-level background sum provably differs from the whole-gene
background sum.
"""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

import ema.switch_test.runner as runner_mod

# ---------------------------------------------------------------------------
# Shared fixture data
# ---------------------------------------------------------------------------

ISOFORM_UTRS = {
    "GENE_X": {
        "TRANSCRIPT_X1": [("chr1", 1000, 1200, "+", 0)],
        "TRANSCRIPT_X2": [("chr1", 2000, 2200, "+", 0)],
    },
    "GENE_Z": {
        "TRANSCRIPT_Z1": [("chr3", 3000, 3500, "+", 0)],
        "TRANSCRIPT_Z2": [("chr3", 3000, 3500, "+", 0), ("chr3", 4000, 4500, "+", 1)],
    },
    "GENE_W": {
        "TRANSCRIPT_W1": [("chr4", 5000, 5200, "+", 0)],
    },
}

# (chrom, start, end, pas_id, score, strand)
PASBED_RECORDS = [
    ("chr1", 1050, 1052, 1, 0, "+"),   # GENE_X / TRANSCRIPT_X1
    ("chr1", 1150, 1152, 2, 0, "+"),   # GENE_X / TRANSCRIPT_X1
    ("chr1", 2050, 2052, 3, 0, "+"),   # GENE_X / TRANSCRIPT_X2
    ("chr1", 2150, 2152, 4, 0, "+"),   # GENE_X / TRANSCRIPT_X2
    ("chr1", 9000, 9002, 5, 0, "+"),   # GENE_X / no UTR overlap (orphan)
    ("chr3", 3200, 3202, 6, 0, "+"),   # GENE_Z / both TRANSCRIPT_Z1 + Z2
    ("chr3", 3250, 3252, 7, 0, "+"),   # GENE_Z / both TRANSCRIPT_Z1 + Z2
    ("chr4", 5100, 5102, 8, 0, "+"),   # GENE_W / TRANSCRIPT_W1 (single UTR)
]

PAS_GENE_ASSIGNMENT = {
    "1": "GENE_X", "2": "GENE_X", "3": "GENE_X", "4": "GENE_X", "5": "GENE_X",
    "6": "GENE_Z", "7": "GENE_Z",
    "8": "GENE_W",
}

CELLS = ["cellA1", "cellA2", "cellB1", "cellB2"]
CLUSTERS = ["A", "A", "B", "B"]

# counts[pas_id] = [cellA1, cellA2, cellB1, cellB2]
COUNTS = {
    "1": [10, 10, 2, 2],   # UTR1 (X1) cluster totals: A=20, B=4
    "2": [5, 5, 5, 5],     # UTR1 (X1) cluster totals: A=10, B=10
    "3": [1, 1, 1, 1],     # UTR2 (X2) cluster totals: A=2, B=2
    "4": [1, 1, 1, 1],     # UTR2 (X2) cluster totals: A=2, B=2
    "5": [3, 3, 3, 3],     # orphan, GENE_X: A=6, B=6
    "6": [4, 4, 1, 1],     # GENE_Z, multi-UTR: A=8, B=2
    "7": [2, 2, 6, 6],     # GENE_Z, multi-UTR: A=4, B=12
    "8": [1, 1, 1, 1],     # GENE_W, single UTR: A=2, B=2
}
# GENE_X UTR1 (PAS1+PAS2) totals: A=30, B=14
# GENE_X UTR2 (PAS3+PAS4) totals: A=4,  B=4
# GENE_X whole-gene (PAS1..5) totals: A=40, B=24


def _write_pasbed(tmp_path: Path) -> None:
    (tmp_path / "pasbed.bed").write_text(
        "".join(f"{c}\t{s}\t{e}\t{p}\t{sc}\t{st}\n" for c, s, e, p, sc, st in PASBED_RECORDS)
    )


def _make_adata() -> ad.AnnData:
    pas_ids = list(PAS_GENE_ASSIGNMENT.keys())
    var = pd.DataFrame({"gene_id": [PAS_GENE_ASSIGNMENT[p] for p in pas_ids]}, index=pas_ids)
    obs = pd.DataFrame({"leiden": CLUSTERS}, index=CELLS)
    X = np.array([[COUNTS[p][ci] for p in pas_ids] for ci in range(len(CELLS))], dtype=float)
    return ad.AnnData(X=X, obs=obs, var=var)


@pytest.fixture
def _run_diff_env(tmp_path, monkeypatch):
    """Wire up run_diff() so it exercises the real bedtools intersect path
    (via map_pas_to_isoforms) but stubs GTF parsing.  Returns a callable
    ``run(isoform_agg, utr_unmatched="gene", with_gtf=True) -> pair_results``.
    """
    _write_pasbed(tmp_path)
    adata = _make_adata()
    gtf_path = tmp_path / "genome.gtf"
    gtf_path.write_text("")  # only needs to exist

    monkeypatch.setattr(runner_mod.ad, "read_h5ad", lambda _p: adata)

    import ema.annotate.gtf2isoform_utr as gtf2isoform_mod
    monkeypatch.setattr(gtf2isoform_mod, "parse_isoform_utrs", lambda *a, **kw: ISOFORM_UTRS)

    def run(isoform_agg: str, utr_unmatched: str = "gene", with_gtf: bool = True,
            out_subdir: str = "out") -> dict:
        return runner_mod.run_diff(
            h5ad_paths=[str(tmp_path / "clustered.h5ad")],
            pasbed=None,
            gtf=str(gtf_path) if with_gtf else None,
            output_dir=str(tmp_path / out_subdir),
            cluster_pairs=None,
            cluster_key="leiden",
            marker_top_n=0,
            marker_method="wilcoxon",
            strategy="fisher",
            fdr=0.05,
            threads=1,
            per_worker_mb=300,
            min_cells_per_group=1,
            isoform_agg=isoform_agg,
            utr_unmatched=utr_unmatched,
        )

    return run


def _pair_df(pair_results: dict) -> pd.DataFrame:
    assert list(pair_results.keys()) == [("A", "B")]
    return pair_results[("A", "B")]


def _rows_for(df: pd.DataFrame, pas_id: str) -> pd.DataFrame:
    return df[df["pas_id"] == pas_id]


# ---------------------------------------------------------------------------
# (a) per_gene: unchanged / byte-identical to today
# ---------------------------------------------------------------------------

class TestPerGeneUnchanged:
    def test_per_gene_explicit_matches_default_omission(self, _run_diff_env):
        """isoform_agg='per_gene' produces the SAME output as leaving the
        argument at its default -- the per_gene code path is untouched."""
        import inspect
        sig = inspect.signature(runner_mod.run_diff)
        assert sig.parameters["isoform_agg"].default == "per_gene"

        df_explicit = _pair_df(_run_diff_env("per_gene", out_subdir="explicit"))
        df_default = _pair_df(_run_diff_env("per_gene", out_subdir="default"))
        pd.testing.assert_frame_equal(
            df_explicit.reset_index(drop=True), df_default.reset_index(drop=True)
        )

    def test_per_gene_has_no_diff_group_id_column(self, _run_diff_env):
        """per_gene never touches _run_grouped_diff, so no 'diff_group_id'
        column is added -- structural proof the path is unchanged."""
        df = _pair_df(_run_diff_env("per_gene"))
        assert "diff_group_id" not in df.columns

    def test_per_gene_background_is_whole_gene(self, _run_diff_env):
        """Under per_gene, PAS 1 (GENE_X) is tested against ALL of GENE_X
        (PAS 1-5), not just its UTR-mates."""
        df = _pair_df(_run_diff_env("per_gene"))
        row = _rows_for(df, "1")
        assert len(row) == 1
        assert int(row["n_reads_gene_cluster1"].iloc[0]) == 40
        assert int(row["n_reads_gene_cluster2"].iloc[0]) == 24


# ---------------------------------------------------------------------------
# (b) within_utr: background is the UTR's other PAS, not the gene's
# ---------------------------------------------------------------------------

class TestWithinUtrScopesToUtr:
    def test_utr1_background_excludes_utr2_and_orphan(self, _run_diff_env):
        df = _pair_df(_run_diff_env("within_utr", utr_unmatched="drop"))
        assert "diff_group_id" in df.columns

        row1 = _rows_for(df, "1")
        assert len(row1) == 1
        assert row1["diff_group_id"].iloc[0] == "GENE_X::TRANSCRIPT_X1"
        # UTR1 total (PAS1+PAS2), NOT the whole-gene total (40/24).
        assert int(row1["n_reads_gene_cluster1"].iloc[0]) == 30
        assert int(row1["n_reads_gene_cluster2"].iloc[0]) == 14

        row3 = _rows_for(df, "3")
        assert len(row3) == 1
        assert row3["diff_group_id"].iloc[0] == "GENE_X::TRANSCRIPT_X2"
        assert int(row3["n_reads_gene_cluster1"].iloc[0]) == 4
        assert int(row3["n_reads_gene_cluster2"].iloc[0]) == 4

    def test_single_pas_utr_gets_no_row(self, _run_diff_env):
        """GENE_W's single UTR has only 1 PAS (PAS 8) -- nothing to
        contrast against, so it must not appear at all."""
        df = _pair_df(_run_diff_env("within_utr", utr_unmatched="drop"))
        assert _rows_for(df, "8").empty


# ---------------------------------------------------------------------------
# (c) utr_unmatched: gene fallback keeps orphans (at gene level); drop omits
# ---------------------------------------------------------------------------

class TestUtrUnmatchedFallback:
    def test_gene_mode_keeps_orphan_tested_at_gene_level(self, _run_diff_env):
        df = _pair_df(_run_diff_env("within_utr", utr_unmatched="gene"))
        row5 = _rows_for(df, "5")
        assert len(row5) == 1
        assert row5["diff_group_id"].iloc[0] == "GENE_X::_gene_"
        # Background = the WHOLE gene (40/24), i.e. genuinely "gene level".
        assert int(row5["n_reads_gene_cluster1"].iloc[0]) == 40
        assert int(row5["n_reads_gene_cluster2"].iloc[0]) == 24

    def test_gene_mode_does_not_duplicate_utr_pas_into_fallback_bucket(self, _run_diff_env):
        """PAS 1-4 (which DO have a UTR) must keep exactly ONE row each
        (from their own UTR group) -- not an extra row from the '_gene_'
        fallback bucket."""
        df = _pair_df(_run_diff_env("within_utr", utr_unmatched="gene"))
        for pas_id in ("1", "2", "3", "4"):
            rows = _rows_for(df, pas_id)
            assert len(rows) == 1, f"PAS {pas_id} got {len(rows)} rows, expected 1"
            assert rows["diff_group_id"].iloc[0].endswith("_gene_") is False

    def test_drop_mode_omits_orphan(self, _run_diff_env):
        df = _pair_df(_run_diff_env("within_utr", utr_unmatched="drop"))
        assert _rows_for(df, "5").empty


# ---------------------------------------------------------------------------
# (d) multi-UTR PAS: tested once per UTR
# ---------------------------------------------------------------------------

class TestMultiUtrPas:
    def test_multi_utr_pas_gets_one_row_per_utr(self, _run_diff_env):
        df = _pair_df(_run_diff_env("within_utr", utr_unmatched="drop"))
        for pas_id in ("6", "7"):
            rows = _rows_for(df, pas_id)
            assert len(rows) == 2, f"PAS {pas_id} got {len(rows)} rows, expected 2"
            assert set(rows["diff_group_id"]) == {
                "GENE_Z::TRANSCRIPT_Z1", "GENE_Z::TRANSCRIPT_Z2",
            }


# ---------------------------------------------------------------------------
# between_utr: differential 3'UTR PREFERENCE (genes with >= 2 UTRs only)
# ---------------------------------------------------------------------------

class TestBetweenUtr:
    def test_two_utr_gene_produces_utr_level_rows(self, _run_diff_env):
        df = _pair_df(_run_diff_env("between_utr", utr_unmatched="drop"))
        # GENE_Z also has 2 UTRs, so its rows are present too -- restrict to
        # GENE_X's own group before checking the UTR-id set.
        gene_x_df = df[df["diff_group_id"] == "GENE_X"]
        assert set(gene_x_df["pas_id"]) == {
            "GENE_X::TRANSCRIPT_X1", "GENE_X::TRANSCRIPT_X2",
        }

        row_utr1 = gene_x_df[gene_x_df["pas_id"] == "GENE_X::TRANSCRIPT_X1"].iloc[0]
        row_utr2 = gene_x_df[gene_x_df["pas_id"] == "GENE_X::TRANSCRIPT_X2"].iloc[0]
        # UTR1's OWN total (30/14) vs the gene's two-UTR background total
        # (30+4=34 / 14+4=18) -- the unit under test is the UTR, and the
        # background is "all UTRs of this gene", not individual PAS.
        assert int(row_utr1["n_reads_pas_cluster1"]) == 30
        assert int(row_utr1["n_reads_pas_cluster2"]) == 14
        assert int(row_utr1["n_reads_gene_cluster1"]) == 34
        assert int(row_utr1["n_reads_gene_cluster2"]) == 18
        assert int(row_utr2["n_reads_pas_cluster1"]) == 4
        assert int(row_utr2["n_reads_pas_cluster2"]) == 4

    def test_single_utr_gene_excluded(self, _run_diff_env):
        """GENE_W has only one 3'UTR -- nothing to contrast, so it must be
        absent from the between_utr output entirely."""
        df = _pair_df(_run_diff_env("between_utr", utr_unmatched="drop"))
        assert "GENE_W" not in set(df["diff_group_id"])


# ---------------------------------------------------------------------------
# Fallback + validation
# ---------------------------------------------------------------------------

class TestFallbackAndValidation:
    def test_missing_gtf_falls_back_to_per_gene(self, _run_diff_env, caplog):
        """within_utr/between_utr without --gtf falls back to per_gene with
        a warning, mirroring run_length's per_isoform fallback."""
        df = _pair_df(_run_diff_env("within_utr", with_gtf=False))
        assert "diff_group_id" not in df.columns

    def test_legacy_per_isoform_alias_behaves_like_within_utr(self, _run_diff_env):
        df_alias = _pair_df(_run_diff_env("per_isoform", utr_unmatched="drop", out_subdir="alias"))
        df_within = _pair_df(_run_diff_env("within_utr", utr_unmatched="drop", out_subdir="within"))
        pd.testing.assert_frame_equal(
            df_alias.reset_index(drop=True), df_within.reset_index(drop=True)
        )

    def test_invalid_isoform_agg_rejected(self, tmp_path, monkeypatch):
        adata = _make_adata()
        monkeypatch.setattr(runner_mod.ad, "read_h5ad", lambda _p: adata)
        with pytest.raises(ValueError, match="isoform_agg"):
            runner_mod.run_diff(
                h5ad_paths=[str(tmp_path / "clustered.h5ad")],
                pasbed=None, gtf=None, output_dir=str(tmp_path / "out"),
                cluster_pairs=None, cluster_key="leiden",
                marker_top_n=0, marker_method="wilcoxon",
                strategy="fisher", fdr=0.05, threads=1, per_worker_mb=300,
                min_cells_per_group=1, isoform_agg="bogus",
            )

    def test_invalid_utr_unmatched_rejected(self, tmp_path, monkeypatch):
        adata = _make_adata()
        monkeypatch.setattr(runner_mod.ad, "read_h5ad", lambda _p: adata)
        with pytest.raises(ValueError, match="utr_unmatched"):
            runner_mod.run_diff(
                h5ad_paths=[str(tmp_path / "clustered.h5ad")],
                pasbed=None, gtf=None, output_dir=str(tmp_path / "out"),
                cluster_pairs=None, cluster_key="leiden",
                marker_top_n=0, marker_method="wilcoxon",
                strategy="fisher", fdr=0.05, threads=1, per_worker_mb=300,
                min_cells_per_group=1, isoform_agg="per_gene",
                utr_unmatched="bogus",
            )
