"""Record-identical regression test for the vectorized per-gene isoform
grouping in ``ema.quantification.strategies.classic``.

Perf pass: ``_gene_per_isoform`` (called once per gene, potentially
thousands of genes) used to build its ``transcript_id -> [(pas_id, rank),
...]`` map with a plain ``for _, row in gene_entries.iterrows(): ...``
loop. Replaced with ``gene_entries.groupby("transcript_id", sort=False)``,
which is order-preserving in the two ways that matter downstream:

  * group (transcript_id) insertion order == first-appearance order in
    ``gene_entries`` — this drives the FINAL output row order, since
    ``_pdui_per_gene_isoform_level`` iterates ``transcript_pas.items()``
    in dict order;
  * within-group row order == original DataFrame row order — this matters
    because ``_pdui_per_gene_isoform_level`` does a *stable*
    ``sorted(ranked_pairs, key=lambda x: x[1])``, so two PAS tied on the
    same ``rank`` must keep their original relative order (their
    proximal/distal choice is order-dependent when ranks tie).

This test builds an interleaved (not pre-grouped) synthetic ``gene_entries``
frame — including a transcript with a tied rank to exercise the
stable-sort dependency — runs both the OLD (iterrows) and NEW (groupby)
``transcript_id -> [(pas_id, rank), ...]`` builders, and asserts they
produce the identical dict (including key order), and that
``_pdui_per_gene_isoform_level`` produces byte-identical rows from each.
It also runs the full public ``ClassicPDUIStrategy.compute(...,
aggregation="per_isoform")`` entry point end-to-end for good measure.
"""
from __future__ import annotations

import pandas as pd

from ema.quantification.strategies.classic import (
    ClassicPDUIStrategy,
    _pdui_per_gene_isoform_level,
)


def _build_t_pas_legacy(gene_entries: pd.DataFrame) -> dict:
    """Verbatim copy of the pre-vectorization per-row loop."""
    t_pas: dict = {}
    for _, row in gene_entries.iterrows():
        tid = row["transcript_id"]
        t_pas.setdefault(tid, []).append((int(row["pas_id"]), int(row["rank"])))
    return t_pas


def _build_t_pas_new(gene_entries: pd.DataFrame) -> dict:
    """The new vectorized builder (copied from classic.py's _gene_per_isoform)."""
    return {
        tid: list(zip(g["pas_id"].astype(int), g["rank"].astype(int)))
        for tid, g in gene_entries.groupby("transcript_id", sort=False)
    }


def _make_gene_entries() -> pd.DataFrame:
    # Deliberately interleaved (not pre-grouped by transcript), and T2 has
    # a tied rank (2) between pas_id 5 and pas_id 7 to exercise stable-sort
    # tie-breaking dependence on original row order.
    return pd.DataFrame(
        [
            {"transcript_id": "T2", "pas_id": 5, "rank": 2},
            {"transcript_id": "T1", "pas_id": 1, "rank": 0},
            {"transcript_id": "T2", "pas_id": 6, "rank": 0},
            {"transcript_id": "T1", "pas_id": 2, "rank": 1},
            {"transcript_id": "T3", "pas_id": 9, "rank": 0},  # single-PAS -> excluded
            {"transcript_id": "T2", "pas_id": 7, "rank": 2},  # ties pas_id 5 on rank=2
        ]
    )


def test_t_pas_builder_matches_legacy_dict_and_order():
    gene_entries = _make_gene_entries()
    legacy = _build_t_pas_legacy(gene_entries)
    new = _build_t_pas_new(gene_entries)

    assert new == legacy
    # Dict key insertion order matters (drives final row order downstream).
    assert list(new.keys()) == list(legacy.keys()) == ["T2", "T1", "T3"]
    # Within-group row order matters for the stable-sort tie-break.
    assert new["T2"] == legacy["T2"] == [(5, 2), (6, 0), (7, 2)]
    assert new["T1"] == legacy["T1"] == [(1, 0), (2, 1)]
    assert new["T3"] == legacy["T3"] == [(9, 0)]


def test_pdui_per_gene_isoform_level_identical_from_both_builders():
    gene_entries = _make_gene_entries()
    t_pas_legacy = _build_t_pas_legacy(gene_entries)
    t_pas_new = _build_t_pas_new(gene_entries)

    count_matrix = pd.DataFrame(
        {
            "cellA": [10.0, 5.0, 3.0, 8.0, 1.0, 4.0, 2.0],
            "cellB": [0.0, 0.0, 6.0, 2.0, 0.0, 9.0, 1.0],
        },
        index=[1, 2, 5, 6, 7, 9, 3],  # pas_ids 1,2,5,6,7,9 valid; 3 unused filler
    )

    rows_legacy = _pdui_per_gene_isoform_level("GENE1", t_pas_legacy, count_matrix, pseudocount=0.0)
    rows_new = _pdui_per_gene_isoform_level("GENE1", t_pas_new, count_matrix, pseudocount=0.0)

    # NaN != NaN under plain dict equality even though it's the same
    # "no signal" value on both sides — use a NaN-aware DataFrame compare.
    pd.testing.assert_frame_equal(pd.DataFrame(rows_new), pd.DataFrame(rows_legacy))
    # T3 excluded (single PAS); T2's tie must resolve proximal=6, distal=7
    # (stable sort on [(5,2),(6,0),(7,2)] -> [(6,0),(5,2),(7,2)]).
    transcripts_seen = {r["transcript_id"] for r in rows_new}
    assert transcripts_seen == {"T1", "T2"}
    t2_rows = [r for r in rows_new if r["transcript_id"] == "T2"]
    assert all(r["proximal_pas_id"] == 6 and r["distal_pas_id"] == 7 for r in t2_rows)
    # Row order: T2 first (first-appearance in gene_entries), then T1.
    assert [r["transcript_id"] for r in rows_new] == ["T2", "T2", "T1", "T1"]


def test_classic_strategy_compute_per_isoform_end_to_end():
    """Full public entry point, small enough to stay on the serial
    (non-Parallel) code path (use_parallel triggers only above 5000 genes).
    """
    # pas_isoform_map: {pas_id: [(gene_id, transcript_id, transcript_pos, rank, total)]}
    pas_isoform_map = {
        1: [("G1", "T1", 0, 0, 2)],
        2: [("G1", "T1", 100, 1, 2)],
        5: [("G2", "T2", 0, 2, 3)],
        6: [("G2", "T2", 10, 0, 3)],
        7: [("G2", "T2", 20, 2, 3)],
        9: [("G3", "T3", 0, 0, 1)],  # single-PAS transcript -> excluded
    }
    count_matrix = pd.DataFrame(
        {
            "cellA": [10.0, 5.0, 3.0, 8.0, 1.0, 4.0],
            "cellB": [0.0, 0.0, 6.0, 2.0, 0.0, 9.0],
        },
        index=[1, 2, 5, 6, 7, 9],
    )

    strategy = ClassicPDUIStrategy()
    out = strategy.compute(
        count_matrix, pas_isoform_map, aggregation="per_isoform", pseudocount=0.0,
    )

    assert set(out["gene_id"]) == {"G1", "G2"}
    assert set(out["transcript_id"]) == {"T1", "T2"}
    g2 = out[out["transcript_id"] == "T2"]
    assert set(g2["proximal_pas_id"]) == {6}
    assert set(g2["distal_pas_id"]) == {7}
