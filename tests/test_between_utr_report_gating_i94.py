"""`between_utr` must withhold a UTR none of whose PAS were selected (#94).

`_build_diff_isoform_groups` decides which UTRs are REPORTED under
``--marker-top-n``::

    reported_utrs = {utr_id for utr_id, members in utr_pas_members.items()
                     if report_pas is None or any(p in report_pas for p in members)}

That comprehension chooses which rows a differential table contains. It was
untested: the fixture in ``test_marker_top_n_utr_denominator_i94.py`` gives
*every* UTR a marker, so replacing the whole expression with
``set(utr_pas_members)`` -- gating nothing -- passed the entire suite.

This file supplies the missing geometry: a third UTR whose member PAS are all
non-markers. It pins both halves of the contract, which pull in opposite
directions:

  a) the marker-free UTR is NOT reported, and
  b) its reads still count towards the gene's background.
"""
from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import ema.switch_test.runner as runner_mod

ISOFORM_UTRS = {
    "GENE_X": {
        "TRANSCRIPT_X1": [("chr1", 1000, 1200, "+", 0)],
        "TRANSCRIPT_X2": [("chr1", 2000, 2200, "+", 0)],
        # No member of X3 is ever selected -- the case that was never exercised.
        "TRANSCRIPT_X3": [("chr1", 3000, 3200, "+", 0)],
    },
}

PASBED_RECORDS = [
    ("chr1", 1050, 1052, 1, 0, "+"),
    ("chr1", 1100, 1102, 2, 0, "+"),
    ("chr1", 2050, 2052, 3, 0, "+"),
    ("chr1", 2150, 2152, 4, 0, "+"),
    ("chr1", 3050, 3052, 5, 0, "+"),   # marker-free UTR
    ("chr1", 3150, 3152, 6, 0, "+"),   # marker-free UTR
]

PAS_IDS = ["1", "2", "3", "4", "5", "6"]
MARKERS = ["1", "3"]                    # one in X1, one in X2, none in X3
UTR1, UTR2, UTR3 = (f"GENE_X::TRANSCRIPT_X{i}" for i in (1, 2, 3))

N_CELLS = 120
COUNT_SEED = 909


def _counts() -> np.ndarray:
    rng = np.random.RandomState(COUNT_SEED)
    X = rng.poisson(0.3, size=(N_CELLS, len(PAS_IDS))).astype(float)
    X[: N_CELLS // 2, 0] += rng.poisson(0.7, size=N_CELLS // 2)
    # Give the marker-free UTR real mass, so dropping it from the background
    # would be plainly visible in the gene totals.
    X[:, 4] += rng.poisson(1.5, size=N_CELLS)
    X[:, 5] += rng.poisson(1.5, size=N_CELLS)
    return X


def _make_adata() -> ad.AnnData:
    var = pd.DataFrame({"gene_id": ["GENE_X"] * len(PAS_IDS)}, index=PAS_IDS)
    obs = pd.DataFrame(
        {"leiden": ["A"] * (N_CELLS // 2) + ["B"] * (N_CELLS // 2)},
        index=[f"c{i}" for i in range(N_CELLS)],
    )
    return ad.AnnData(X=_counts(), obs=obs, var=var)


def _all_pas_read_sums() -> tuple[int, int]:
    X = _counts()
    half = N_CELLS // 2
    return int(X[:half].sum()), int(X[half:].sum())


@pytest.fixture(scope="module")
def run_between(tmp_path_factory):
    tmp_path = tmp_path_factory.mktemp("i94_between_gating")
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

    def run(marker_top_n: int) -> pd.DataFrame:
        pair_results = runner_mod.run_diff(
            h5ad_paths=[str(tmp_path / "clustered.h5ad")],
            pasbed=None,
            gtf=str(gtf_path),
            output_dir=str(tmp_path / f"between_top{marker_top_n}"),
            cluster_pairs=None,
            cluster_key="leiden",
            marker_top_n=marker_top_n,
            marker_method="wilcoxon",
            strategy="fisher",
            fdr=0.05,
            threads=1,
            per_worker_mb=300,
            min_cells_per_group=5,
            isoform_agg="between_utr",
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


def test_marker_free_utr_is_not_reported(run_between):
    """(a) `--marker-top-n` must actually gate which UTRs come back."""
    off = set(run_between(0)["pas_id"])
    on = set(run_between(len(MARKERS))["pas_id"])

    assert {UTR1, UTR2, UTR3} <= off, (
        f"with selection OFF every UTR should be reported, got {sorted(off)}"
    )
    assert UTR3 not in on, (
        f"{UTR3} has no selected member PAS but was still reported: "
        f"{sorted(on)}. The reported_utrs filter in "
        "_build_diff_isoform_groups is not gating anything."
    )
    assert {UTR1, UTR2} <= on, (
        f"UTRs WITH a selected member must still be reported, got {sorted(on)}"
    )


def test_marker_free_utr_still_counts_towards_the_background(run_between):
    """(b) ...while its reads stay in the gene's denominator (the #94 fix)."""
    on = run_between(len(MARKERS)).set_index("pas_id")
    expected_a, expected_b = _all_pas_read_sums()
    for utr in (UTR1, UTR2):
        got = (int(on.loc[utr, "n_reads_gene_cluster1"]),
               int(on.loc[utr, "n_reads_gene_cluster2"]))
        assert got == (expected_a, expected_b), (
            f"{utr}: gene background is {got}, expected "
            f"{(expected_a, expected_b)} (the sum over ALL PAS including the "
            "unreported UTR). Withholding a UTR from the report must not "
            "remove it from the background."
        )
