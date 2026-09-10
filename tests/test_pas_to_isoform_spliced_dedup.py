"""Regression tests: a spliced 3'UTR must not duplicate a PAS in the isoform map.

A spliced 3'UTR is annotated by one ``three_prime_utr`` record per exon.  When a
PAS falls inside more than one record of the SAME transcript, ``bedtools
intersect`` reports it once per record and ``map_pas_to_isoforms`` used to emit
the same ``(pas_id, gene_id, transcript_id)`` entry more than once -- giving that
one PAS two different ``rank`` values and an inflated
``total_pas_in_transcript``.

PR #112 de-duplicated this at ONE consumer (``switch diff``'s
``_build_diff_isoform_groups``).  Every other consumer read the map raw and was
silently wrong; the fix now lives in ``map_pas_to_isoforms`` itself, so these
tests exercise the three PDUI strategies through the real map.

Fixture (``chr7``, ``+`` strand)
--------------------------------
``GENE_S / TRANSCRIPT_S1`` has a 3'UTR annotated by TWO ``three_prime_utr``
records that both cover chr7:6100-6200::

    record 0: chr7:6000-6200
    record 1: chr7:6100-6400

    PAS 1 (chr7:6150) -> inside BOTH records -> reported twice by bedtools
    PAS 3 (chr7:6300) -> inside record 1 only -> reported once

One cell, ``cell1``, with 10 reads on PAS 1 and 90 on PAS 3: the honest totals
are 2 distinct PAS / 100 reads, proportions 0.1 / 0.9 and Shannon entropy
0.469 bits.  Counting PAS 1 twice gives 3 "PAS" / 110 reads, proportions
0.0909 / 0.0909 / 0.818 and entropy 0.866 bits.

``GENE_U / TRANSCRIPT_U1`` is the unspliced control: same PAS, one UTR record,
so its output must be identical before and after the fix.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from ema.quantification.pas_to_isoform import map_pas_to_isoforms
from ema.quantification.strategies import get_pdui_strategy
from ema.quantification.strategies.classic import _build_pas_info_from_map


# ---------------------------------------------------------------------------
# Fixture data
# ---------------------------------------------------------------------------

SPLICED_UTRS = {
    "GENE_S": {
        "TRANSCRIPT_S1": [
            ("chr7", 6000, 6200, "+", 0),
            ("chr7", 6100, 6400, "+", 1),
        ],
    },
}

UNSPLICED_UTRS = {
    "GENE_U": {
        "TRANSCRIPT_U1": [
            ("chr7", 6000, 6400, "+", 0),
        ],
    },
}

# (chrom, start, end, pas_id, score, strand)
PASBED_RECORDS = [
    ("chr7", 6150, 6152, 1, 0, "+"),   # inside BOTH spliced records
    ("chr7", 6300, 6302, 3, 0, "+"),   # inside the distal record only
]

COUNTS = pd.DataFrame({"cell1": [10.0, 90.0]}, index=[1, 3])

# -(0.1*log2 0.1 + 0.9*log2 0.9)
CORRECT_ENTROPY = 0.4689955935892812
# What counting PAS 1 twice produced (10/110, 10/110, 90/110).
INFLATED_ENTROPY = 0.8661525353626755


@pytest.fixture
def pasbed(tmp_path) -> Path:
    p = tmp_path / "pas.bed"
    p.write_text(
        "".join(f"{c}\t{s}\t{e}\t{i}\t{sc}\t{st}\n"
                for c, s, e, i, sc, st in PASBED_RECORDS)
    )
    return p


@pytest.fixture
def spliced_map(pasbed) -> dict:
    return map_pas_to_isoforms(pasbed, SPLICED_UTRS)


# ---------------------------------------------------------------------------
# The map itself
# ---------------------------------------------------------------------------

def test_fixture_really_intersects_two_utr_records(pasbed):
    """Guard on the fixture: without the double intersection this is vacuous."""
    from ema.quantification.pas_to_isoform import (
        _build_utr_bed,
        _run_bedtools_intersect,
    )

    raw = _run_bedtools_intersect(pasbed, _build_utr_bed(SPLICED_UTRS))
    assert len(raw[raw["pas_id"] == 1]) == 2, raw


def test_map_emits_each_transcript_once_per_pas(spliced_map):
    named = [(g, t) for g, t, *_ in spliced_map[1]]
    assert named == [("GENE_S", "TRANSCRIPT_S1")]


def test_map_rank_and_total_use_the_distinct_pas_set(spliced_map):
    ranks = {pas_id: entries[0][3] for pas_id, entries in spliced_map.items()}
    totals = {pas_id: entries[0][4] for pas_id, entries in spliced_map.items()}
    assert ranks == {1: 0, 3: 1}          # one rank each, proximal -> distal
    assert totals == {1: 2, 3: 2}         # 2 distinct PAS, not 3


# ---------------------------------------------------------------------------
# Consumers
# ---------------------------------------------------------------------------

def test_shannon_entropy_is_not_inflated_by_the_spliced_utr(spliced_map):
    """Was pas_ids='1;1;3', n_pas=3, reads=110, entropy=0.866."""
    out = get_pdui_strategy("shannon").compute(
        COUNTS, spliced_map, aggregation="per_isoform"
    )
    row = out[out["transcript_id"] == "TRANSCRIPT_S1"].iloc[0]
    assert row["pas_ids"] == "1;3"
    assert row["n_pas"] == 2
    assert row["total_reads_transcript"] == 100.0
    assert row["entropy"] == pytest.approx(CORRECT_ENTROPY, abs=1e-9)
    assert abs(row["entropy"] - INFLATED_ENTROPY) > 1e-3


def test_proportion_counts_each_pas_once(spliced_map):
    """Was two rows for PAS 1, each with proportion 10/110 = 0.0909."""
    out = get_pdui_strategy("proportion").compute(
        COUNTS, spliced_map, aggregation="per_isoform"
    )
    iso = out[out["transcript_id"] == "TRANSCRIPT_S1"]
    assert len(iso) == 2                              # one row per PAS, not 3
    assert iso["pas_id"].tolist() == [1, 3]
    prox = iso[iso["pas_id"] == 1].iloc[0]
    assert prox["proportion"] == pytest.approx(0.1, abs=1e-12)
    assert iso["proportion"].sum() == pytest.approx(1.0, abs=1e-12)


def test_classic_pas_info_rank_and_total_are_not_inflated(spliced_map):
    """``_build_pas_info_from_map`` feeds every classic/PDUI code path.

    Was three rows -- PAS 1 twice, with ranks 0 AND 1 and
    ``total_pas_in_transcript`` = 3, so the transcript looked like a 3-PAS UTR.
    """
    info = _build_pas_info_from_map(spliced_map)
    assert len(info) == 2
    assert not info.duplicated(
        subset=["pas_id", "gene_id", "transcript_id"]
    ).any()
    assert info["total_pas_in_transcript"].tolist() == [2, 2]
    assert info[info["pas_id"] == 1]["rank"].tolist() == [0]

    pdui = get_pdui_strategy("classic").compute(
        COUNTS, spliced_map, aggregation="per_isoform"
    )
    row = pdui.iloc[0]
    assert (row["proximal_pas_id"], row["distal_pas_id"]) == (1, 3)
    assert row["pdui"] == pytest.approx(0.9, abs=1e-12)
    assert row["total_reads"] == 100.0


# ---------------------------------------------------------------------------
# Unspliced control: the fix must be a no-op here
# ---------------------------------------------------------------------------

def test_unspliced_transcript_is_unchanged(pasbed):
    """Same PAS, one UTR record: nothing to de-duplicate."""
    m = map_pas_to_isoforms(pasbed, UNSPLICED_UTRS)
    assert m == {
        1: [("GENE_U", "TRANSCRIPT_U1", 151, 0, 2)],
        3: [("GENE_U", "TRANSCRIPT_U1", 301, 1, 2)],
    }
    out = get_pdui_strategy("shannon").compute(
        COUNTS, m, aggregation="per_isoform"
    )
    row = out.iloc[0]
    assert row["n_pas"] == 2
    assert row["total_reads_transcript"] == 100.0
    assert row["entropy"] == pytest.approx(CORRECT_ENTROPY, abs=1e-9)
