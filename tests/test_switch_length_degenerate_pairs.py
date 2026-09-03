"""`ema switch length --isoform-agg per_isoform` — degenerate proximal==distal
pairs (issue #98).

A per-isoform unit whose proximal and distal endpoints resolve to the SAME
PAS carries no length information: ``pdui = distal/(proximal+distal)`` is
0.5 by construction, ``proximal_reads == distal_reads``, and every
downstream shorten/lengthen call built on that row is meaningless.  Such
units appeared whenever one PAS produced more than one row for a single
(gene, transcript) — one bedtools row per overlapping UTR exon — so the
transcript looked like it had >= 2 PAS while having only one usable site
(18,116 / 13,709,930 rows over 7 genes on GSE104556 mouse1).

Two defects are covered:

(a) PRODUCER.  ``ClassicPDUIStrategy`` (per_isoform) must collapse repeated
    PAS before choosing the endpoints and then EXCLUDE the transcript when a
    single distinct PAS is left — the same treatment any other ``<2`` PAS
    transcript already gets.  Dropping is preferable to flagging: the row has
    no PDUI to report, so keeping it would force every consumer of
    ``pdui_classic.tsv`` to learn about a new column or silently average
    0.5s into their shortening scores.

(b) GUARD.  ``_assert_pdui_strand_convention`` classified equality with
    ``~(proximal_start < distal_start)``, so a degenerate pair landed in the
    "inverted / strand" bucket and the error blamed strand selection — the
    wrong cause (there were ZERO true strand inversions in the reported run).
    Degenerate rows must be counted and named separately from inverted ones.

Both are exercised on BOTH strands, with a normal pair alongside the
degenerate one in the same gene, so the fix cannot work by dropping
everything.
"""
from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from ema.quantification.strategies.classic import ClassicPDUIStrategy
from ema.switch_test.runner import _assert_pdui_strand_convention, run_length

CELLS = ["c0", "c1"]

# Mirror-image loci.  Ranks below are transcript-coordinate ranks (0 =
# proximal), i.e. already strand-aware — that is exactly what the
# per_isoform map hands the strategy.
#
#   GENE_PLUS  (+): pas 1 @1000 proximal, pas 3 @3000 distal
#   GENE_MINUS (-): pas 6 @3000 proximal, pas 4 @1000 distal
#
# Each gene also carries a DEGENERATE transcript: one PAS (2 / 5) that
# produced two rows for the same transcript, so it is ranked as if it were
# both endpoints.
PASBED = [
    ("chr1", 1000, 1001, "1", 0, "+"),
    ("chr1", 2000, 2001, "2", 0, "+"),
    ("chr1", 3000, 3001, "3", 0, "+"),
    ("chr2", 1000, 1001, "4", 0, "-"),
    ("chr2", 2000, 2001, "5", 0, "-"),
    ("chr2", 3000, 3001, "6", 0, "-"),
]

# 10 reads on each true proximal PAS, 90 on each true distal PAS -> the
# correct PDUI of both normal pairs is 0.9.  The degenerate PAS carry 50
# reads so a degenerate row would show the issue's exact signature
# (proximal_reads == distal_reads, pdui == 0.5).
COUNTS = np.array([[10, 50, 90, 90, 50, 10]] * len(CELLS), dtype=np.int64)

PAS_GENE = {
    "1": "GENE_PLUS", "2": "GENE_PLUS", "3": "GENE_PLUS",
    "4": "GENE_MINUS", "5": "GENE_MINUS", "6": "GENE_MINUS",
}


def _map_with_degenerate() -> dict[int, list[tuple[str, str, int, int, int]]]:
    """``{pas_id: [(gene_id, transcript_id, transcript_pos, rank, total)]}``.

    ``T_PLUS_DEG`` / ``T_MINUS_DEG`` hold the reported pathology: the SAME
    PAS listed twice for one transcript (two overlapping UTR-exon rows), so
    it is ranked 0 and 1 and would be picked as both endpoints.
    """
    return {
        1: [("GENE_PLUS", "T_PLUS_OK", 0, 0, 2)],
        3: [("GENE_PLUS", "T_PLUS_OK", 200, 1, 2)],
        2: [("GENE_PLUS", "T_PLUS_DEG", 0, 0, 2),
            ("GENE_PLUS", "T_PLUS_DEG", 0, 1, 2)],
        6: [("GENE_MINUS", "T_MINUS_OK", 0, 0, 2)],
        4: [("GENE_MINUS", "T_MINUS_OK", 200, 1, 2)],
        5: [("GENE_MINUS", "T_MINUS_DEG", 0, 0, 2),
            ("GENE_MINUS", "T_MINUS_DEG", 0, 1, 2)],
    }


def _count_matrix() -> pd.DataFrame:
    return pd.DataFrame(
        COUNTS.T.astype(float),
        index=[1, 2, 3, 4, 5, 6],
        columns=CELLS,
    )


# --- (a) producer: the degenerate unit never reaches the table -------------

def test_per_isoform_drops_degenerate_pairs_on_both_strands():
    out = ClassicPDUIStrategy().compute(
        count_matrix=_count_matrix(),
        pas_isoform_map=_map_with_degenerate(),
        aggregation="per_isoform",
        pseudocount=0.0,
    )

    # The degenerate transcripts are gone, on BOTH strands ...
    assert set(out["transcript_id"]) == {"T_PLUS_OK", "T_MINUS_OK"}
    assert not (out["proximal_pas_id"] == out["distal_pas_id"]).any()
    assert not (out["pdui"] == 0.5).any()
    assert not (out["proximal_reads"] == out["distal_reads"]).any()

    # ... and the normal pair of each gene is untouched: proximal/distal
    # still strand-aware, PDUI still 0.9.
    plus = out[out["transcript_id"] == "T_PLUS_OK"]
    minus = out[out["transcript_id"] == "T_MINUS_OK"]
    assert len(plus) == len(minus) == len(CELLS)
    assert set(plus["proximal_pas_id"]) == {1} and set(plus["distal_pas_id"]) == {3}
    assert set(minus["proximal_pas_id"]) == {6} and set(minus["distal_pas_id"]) == {4}
    assert np.allclose(plus["pdui"], 0.9)
    assert np.allclose(minus["pdui"], 0.9)


def test_duplicate_pas_rows_do_not_hide_a_real_distal_pas():
    """A duplicated PAS alongside a genuine second PAS stays a real pair.

    Collapsing must key on ``pas_id`` (keeping the most proximal rank), not
    just reject transcripts with duplicate rows.
    """
    pas_isoform_map = {
        1: [("GENE_PLUS", "T_DUP", 0, 0, 3),
            ("GENE_PLUS", "T_DUP", 0, 2, 3)],   # same PAS, ranked twice
        3: [("GENE_PLUS", "T_DUP", 200, 1, 3)],
    }
    out = ClassicPDUIStrategy().compute(
        count_matrix=_count_matrix(), pas_isoform_map=pas_isoform_map,
        aggregation="per_isoform", pseudocount=0.0,
    )
    assert set(out["transcript_id"]) == {"T_DUP"}
    assert set(out["proximal_pas_id"]) == {1}
    assert set(out["distal_pas_id"]) == {3}
    assert np.allclose(out["pdui"], 0.9)


# --- (b) guard: degenerate is reported as degenerate, not as a strand bug --

def test_strand_guard_names_degenerate_rows_instead_of_blaming_strand():
    from_both_strands = pd.DataFrame({
        "gene_id": ["GENE_PLUS", "GENE_MINUS"],
        "transcript_id": ["T_PLUS_DEG", "T_MINUS_DEG"],
        "proximal_pas_id": [2, 5],
        "distal_pas_id": [2, 5],
        "proximal_start": [2000, 2000],
        "distal_start": [2000, 2000],
        "proximal_strand": ["+", "-"],
        "distal_strand": ["+", "-"],
    })
    with pytest.raises(RuntimeError) as exc:
        _assert_pdui_strand_convention(from_both_strands, "bad.tsv")
    msg = str(exc.value)
    assert "2 DEGENERATE" in msg
    assert "0 INVERTED" in msg
    assert "0 with proximal and distal on DIFFERENT strands" in msg


def test_strand_guard_still_reports_a_true_inversion_as_inverted():
    inverted = pd.DataFrame({
        "gene_id": ["GENE_MINUS"],
        "proximal_pas_id": [4],
        "distal_pas_id": [6],
        "proximal_start": [1000],   # '-' strand: proximal must be the HIGHER
        "distal_start": [3000],
        "proximal_strand": ["-"],
        "distal_strand": ["-"],
    })
    with pytest.raises(RuntimeError) as exc:
        _assert_pdui_strand_convention(inverted, "bad.tsv")
    msg = str(exc.value)
    assert "1 INVERTED" in msg
    assert "0 DEGENERATE" in msg


# --- end-to-end: `run_length --isoform-agg per_isoform` completes ----------

# UTR model handed to the (stubbed) GTF parser.  ``T_PLUS_DEG`` /
# ``T_MINUS_DEG`` have two OVERLAPPING UTR-exon records, so the real
# `bedtools intersect -wa -wb -s` inside map_pas_to_isoforms emits two rows
# for the single PAS they contain — the exact input that produced the
# degenerate pairs on the testis runs.
ISOFORM_UTRS = {
    "GENE_PLUS": {
        "T_PLUS_OK": [("chr1", 900, 3100, "+", 0)],
        "T_PLUS_DEG": [("chr1", 1900, 2100, "+", 0),
                       ("chr1", 1900, 2100, "+", 1)],
    },
    "GENE_MINUS": {
        "T_MINUS_OK": [("chr2", 900, 3100, "-", 0)],
        "T_MINUS_DEG": [("chr2", 1900, 2100, "-", 0),
                        ("chr2", 1900, 2100, "-", 1)],
    },
}


def test_run_length_per_isoform_completes_and_writes_no_degenerate_rows(
    tmp_path, monkeypatch,
):
    pas_ids = list(PAS_GENE.keys())
    adata = ad.AnnData(
        X=COUNTS.astype(np.float64),
        obs=pd.DataFrame({"leiden": ["0", "1"]}, index=CELLS),
        var=pd.DataFrame({"gene_id": [PAS_GENE[p] for p in pas_ids]},
                         index=pas_ids),
    )
    h5 = tmp_path / "clusters.h5ad"
    adata.write_h5ad(h5)
    pd.DataFrame(PASBED).to_csv(
        tmp_path / "pasbed.bed", sep="\t", header=False, index=False,
    )
    # Keeps run_length's GTF cache inside tmp_path instead of ~/.cache.
    (tmp_path / "run_config.json").write_text("{}")
    gtf = tmp_path / "genome.gtf"
    gtf.write_text("")  # only has to exist; the parser is stubbed

    import ema.annotate.gtf2isoform_utr as gtf2isoform_mod
    monkeypatch.setattr(
        gtf2isoform_mod, "parse_isoform_utrs", lambda *a, **kw: ISOFORM_UTRS,
    )

    out_dir = tmp_path / "out"
    run_length(
        h5ad_paths=[str(h5)], gtf=str(gtf), output_dir=str(out_dir),
        cluster_pairs=None, cluster_key="leiden", strategy="classic",
        isoform_agg="per_isoform", isoform_collapse="none", threads=1,
        utr_unmatched="drop", pasbed=str(tmp_path / "pasbed.bed"),
    )

    df = pd.read_csv(out_dir / "pdui_classic.tsv", sep="\t")
    assert not df.empty
    assert set(df["transcript_id"]) == {"T_PLUS_OK", "T_MINUS_OK"}
    assert not (df["proximal_pas_id"] == df["distal_pas_id"]).any()
    assert not (df["proximal_start"] == df["distal_start"]).any()
    assert not (df["pdui"] == 0.5).any()
    plus = df[df["gene_id"] == "GENE_PLUS"]
    minus = df[df["gene_id"] == "GENE_MINUS"]
    assert (plus["proximal_start"] < plus["distal_start"]).all()
    assert (minus["proximal_start"] > minus["distal_start"]).all()
    assert np.allclose(df["pdui"], 0.9)
