"""Regression tests for `switch diff --isoform-agg within_utr` on a SPLICED 3'UTR.

`_build_diff_isoform_groups` builds each 3'UTR's member list by walking the
`(gene_id, transcript_id, ...)` entries `map_pas_to_isoforms` returns for a PAS.
A spliced 3'UTR contributes several `three_prime_utr` records for ONE transcript,
so those entries can name the same `(gene, transcript)` more than once and the
same `pas_id` was appended repeatedly to a single UTR group.

The duplicate columns then reach the fisher strategy through
`_run_grouped_diff`'s `test_matrix[bg_cols]`, where `int(agg1[p])` receives a
Series instead of a scalar and raises
``TypeError: cannot convert the series to <class 'int'>`` -- so every
`--isoform-agg within_utr` run over a GTF with a spliced 3'UTR aborts.

This is the remaining half of PR #109; the `NameError` half of that PR is
already on develop (commit 7b50a5f).  `tests/test_run_diff_isoform.py` covers
the clean, unspliced `within_utr` / `between_utr` / `utr_unmatched` semantics --
only the duplicate-member case is added here.

`map_pas_to_isoforms` now de-duplicates at the source as well (see
`tests/test_pas_to_isoform_spliced_dedup.py`), which fixes the strategies that
read the map raw.  The consumer-side guard exercised here is kept as defence in
depth, so these tests still pin the `switch diff` behaviour.

Fixture
-------
GENE_S / TRANSCRIPT_S1 has a 3'UTR annotated by TWO `three_prime_utr` records
that both cover chr5:6100-6200:

    record 0: chr5:6000-6200
    record 1: chr5:6100-6400

    PAS 10 (chr5:6150) -> intersects BOTH records -> the SAME (GENE_S,
                          TRANSCRIPT_S1) is named twice in its entry list
    PAS 11 (chr5:6300) -> intersects record 1 only -> named once

The UTR group must therefore hold ``[10, 11]``, not ``[10, 10, 11]``.

Two clusters "A" (cellA1, cellA2) and "B" (cellB1, cellB2).
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
# Fixture data
# ---------------------------------------------------------------------------

ISOFORM_UTRS = {
    "GENE_S": {
        # Spliced 3'UTR: two three_prime_utr records for one transcript, both
        # covering chr5:6100-6200.
        "TRANSCRIPT_S1": [
            ("chr5", 6000, 6200, "+", 0),
            ("chr5", 6100, 6400, "+", 1),
        ],
    },
}

# (chrom, start, end, pas_id, score, strand)
PASBED_RECORDS = [
    ("chr5", 6150, 6152, 10, 0, "+"),   # inside BOTH UTR records
    ("chr5", 6300, 6302, 11, 0, "+"),   # inside record 1 only
]

PAS_GENE_ASSIGNMENT = {"10": "GENE_S", "11": "GENE_S"}

CELLS = ["cellA1", "cellA2", "cellB1", "cellB2"]
CLUSTERS = ["A", "A", "B", "B"]

# counts[pas_id] = [cellA1, cellA2, cellB1, cellB2]
COUNTS = {
    "10": [8, 8, 1, 1],   # cluster totals: A=16, B=2
    "11": [2, 2, 7, 7],   # cluster totals: A=4,  B=14
}
# De-duplicated UTR background (PAS10 + PAS11): A=20, B=16.
# Had PAS 10 been counted twice it would have been A=36, B=18.
UTR_BG_C1 = 20
UTR_BG_C2 = 16


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
def _spliced_utr_env(tmp_path, monkeypatch):
    """Real ``bedtools intersect`` (via ``map_pas_to_isoforms``); only GTF
    parsing is stubbed, mirroring ``tests/test_run_diff_isoform.py``."""
    _write_pasbed(tmp_path)
    adata = _make_adata()
    gtf_path = tmp_path / "genome.gtf"
    gtf_path.write_text("")  # only needs to exist

    monkeypatch.setattr(runner_mod.ad, "read_h5ad", lambda _p: adata)

    import ema.annotate.gtf2isoform_utr as gtf2isoform_mod
    monkeypatch.setattr(gtf2isoform_mod, "parse_isoform_utrs", lambda *a, **kw: ISOFORM_UTRS)

    return tmp_path, adata, gtf_path


def test_spliced_utr_intersects_the_transcript_twice_but_the_map_dedups(
    _spliced_utr_env,
):
    """Guard on the fixture itself: without it the regression below is vacuous.

    The fixture must still make ``bedtools`` report PAS 10 once per
    ``three_prime_utr`` record -- that is the condition the consumer-side guard
    below protects against.  ``map_pas_to_isoforms`` now drops the duplicate at
    the source (so every other consumer is fixed too), which is why the map
    itself names the transcript once.
    """
    from ema.quantification.pas_to_isoform import (
        _build_utr_bed,
        _run_bedtools_intersect,
        map_pas_to_isoforms,
    )

    tmp_path, _adata, _gtf = _spliced_utr_env
    pasbed = tmp_path / "pasbed.bed"

    hits = _run_bedtools_intersect(pasbed, _build_utr_bed(ISOFORM_UTRS))
    assert len(hits[hits["pas_id"] == 10]) == 2, hits

    raw = map_pas_to_isoforms(pasbed, ISOFORM_UTRS)
    named = [(g, t) for g, t, *_ in raw[10]]
    assert named.count(("GENE_S", "TRANSCRIPT_S1")) == 1, named


def test_spliced_utr_does_not_duplicate_a_pas_in_its_group(_spliced_utr_env):
    """The UTR group holds each PAS once, in first-seen order."""
    tmp_path, adata, gtf_path = _spliced_utr_env
    diff_df = pd.DataFrame(
        np.array([[COUNTS[p][ci] for p in ("10", "11")] for ci in range(len(CELLS))], dtype=float),
        index=CELLS,
        columns=[10, 11],
    )

    test_matrix, groups = runner_mod._build_diff_isoform_groups(
        adata=adata,
        diff_df=diff_df,
        h5ad_path=str(tmp_path / "clustered.h5ad"),
        gtf=str(gtf_path),
        pasbed_path=tmp_path / "pasbed.bed",
        isoform_agg="within_utr",
        utr_unmatched="drop",
        n_jobs=1,
    )

    assert [g["group_id"] for g in groups] == ["GENE_S::TRANSCRIPT_S1"]
    group = groups[0]
    assert group["bg_cols"] == [10, 11], (
        f"UTR background is {group['bg_cols']}; a spliced 3'UTR must not append a "
        "PAS twice -- duplicate columns raise TypeError in the fisher strategy."
    )
    assert group["report_cols"] == [10, 11]
    assert test_matrix is diff_df


def test_within_utr_runs_and_uses_the_deduplicated_background(_spliced_utr_env):
    """End-to-end: `within_utr` completes over a spliced 3'UTR (no TypeError)
    and the reported background is the de-duplicated UTR total."""
    tmp_path, _adata, gtf_path = _spliced_utr_env

    pair_results = runner_mod.run_diff(
        h5ad_paths=[str(tmp_path / "clustered.h5ad")],
        pasbed=None,
        gtf=str(gtf_path),
        output_dir=str(tmp_path / "out"),
        cluster_pairs=None,
        cluster_key="leiden",
        marker_top_n=0,
        marker_method="wilcoxon",
        strategy="fisher",
        fdr=0.05,
        threads=1,
        per_worker_mb=300,
        min_cells_per_group=1,
        isoform_agg="within_utr",
        utr_unmatched="drop",
    )

    df = pair_results[("A", "B")]
    assert set(df["diff_group_id"]) == {"GENE_S::TRANSCRIPT_S1"}
    for pas_id in ("10", "11"):
        rows = df[df["pas_id"] == pas_id]
        assert len(rows) == 1, f"PAS {pas_id} got {len(rows)} rows, expected 1"
        assert int(rows["n_reads_gene_cluster1"].iloc[0]) == UTR_BG_C1
        assert int(rows["n_reads_gene_cluster2"].iloc[0]) == UTR_BG_C2

    row10 = df[df["pas_id"] == "10"].iloc[0]
    assert int(row10["n_reads_pas_cluster1"]) == 16
    assert int(row10["n_reads_pas_cluster2"]) == 2


# ---------------------------------------------------------------------------
# between_utr: a spliced UTR was silently double-counted (needs TWO isoforms,
# because `between_utr` compares the UTRs of one gene against each other).
# ---------------------------------------------------------------------------

TWO_ISOFORM_UTRS = {
    "GENE_T": {
        # Spliced: two three_prime_utr records, both covering PAS 20.
        "TRANSCRIPT_T1": [
            ("chr9", 1000, 1200, "+", 0),
            ("chr9", 1100, 1400, "+", 1),
        ],
        # Unspliced control: one record, covering PAS 21 only.
        "TRANSCRIPT_T2": [
            ("chr9", 2000, 2400, "+", 0),
        ],
    },
}

TWO_ISOFORM_PASBED = [
    ("chr9", 1150, 1152, 20, 0, "+"),   # inside BOTH records of T1
    ("chr9", 2100, 2102, 21, 0, "+"),   # inside T2 only
]

TWO_ISOFORM_COUNTS = {
    "20": [8, 8, 1, 1],   # cluster totals: A=16, B=2
    "21": [2, 2, 7, 7],   # cluster totals: A=4,  B=14
}


@pytest.fixture
def _two_isoform_env(tmp_path, monkeypatch):
    (tmp_path / "pasbed.bed").write_text(
        "".join(f"{c}\t{s}\t{e}\t{p}\t{sc}\t{st}\n"
                for c, s, e, p, sc, st in TWO_ISOFORM_PASBED)
    )
    pas_ids = list(TWO_ISOFORM_COUNTS.keys())
    var = pd.DataFrame({"gene_id": ["GENE_T"] * len(pas_ids)}, index=pas_ids)
    obs = pd.DataFrame({"leiden": CLUSTERS}, index=CELLS)
    X = np.array([[TWO_ISOFORM_COUNTS[p][ci] for p in pas_ids]
                  for ci in range(len(CELLS))], dtype=float)
    adata = ad.AnnData(X=X, obs=obs, var=var)

    gtf_path = tmp_path / "genome.gtf"
    gtf_path.write_text("")
    monkeypatch.setattr(runner_mod.ad, "read_h5ad", lambda _p: adata)
    import ema.annotate.gtf2isoform_utr as gtf2isoform_mod
    monkeypatch.setattr(gtf2isoform_mod, "parse_isoform_utrs",
                        lambda *a, **kw: TWO_ISOFORM_UTRS)
    return tmp_path, gtf_path


def test_between_utr_does_not_double_count_a_spliced_utr(_two_isoform_env):
    """`between_utr` never raised on a spliced UTR -- it silently DOUBLE-COUNTED.

    The row *is* the UTR isoform, so a PAS falling inside two
    ``three_prime_utr`` records of one transcript was added to that UTR twice
    and its reads counted twice.  Nothing failed; the totals were just
    inflated, which is the worst kind of bug to ship untested.  PAS 20
    (A=16, B=2) is the duplicated one; TRANSCRIPT_T1's honest totals are
    therefore A=16 / B=2, and before the de-duplication they were A=32 / B=4.
    """
    tmp_path, gtf_path = _two_isoform_env

    pair_results = runner_mod.run_diff(
        h5ad_paths=[str(tmp_path / "clustered.h5ad")],
        pasbed=None,
        gtf=str(gtf_path),
        output_dir=str(tmp_path / "out_between"),
        cluster_pairs=None,
        cluster_key="leiden",
        marker_top_n=0,
        marker_method="wilcoxon",
        strategy="fisher",
        fdr=0.05,
        threads=1,
        per_worker_mb=300,
        min_cells_per_group=1,
        isoform_agg="between_utr",
        utr_unmatched="drop",
    )

    df = pair_results[("A", "B")]
    rows = df[df["pas_id"] == "GENE_T::TRANSCRIPT_T1"]
    assert len(rows) == 1, f"expected one row for the spliced UTR, got {len(rows)}"
    row = rows.iloc[0]

    assert int(row["n_reads_pas_cluster1"]) == 16, (
        f"spliced UTR reads in cluster A = {int(row['n_reads_pas_cluster1'])}, "
        "expected 16; 32 means PAS 20 was counted once per three_prime_utr record"
    )
    assert int(row["n_reads_pas_cluster2"]) == 2, (
        f"spliced UTR reads in cluster B = {int(row['n_reads_pas_cluster2'])}, "
        "expected 2; 4 means PAS 20 was counted twice"
    )
