"""E5: findings_long transform conforms to the FindingRow contract."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from ema.switch_test.long_output import (
    FINDING_LONG_COLUMNS,
    classify_direction,
    finding_uid,
    findings_long,
    write_long_table,
)


def test_finding_uid_format_and_sanitization() -> None:
    uid = finding_uid("switch_diff:fisher", "fisher", "clu:A", "12", "chr1:100:+")
    # ':' in fields (arm/strategy/label/var) is replaced; pas_uid kept as tail.
    assert uid == "switch_diff_fisher:fisher:clu_A:12:chr1:100:+"


def test_classify_direction_is_exhaustive_and_never_na() -> None:
    assert classify_direction(0.2, 0.3, fdr=0.05) == "flat"        # not significant
    assert classify_direction(0.01, 0.3, fdr=0.05) == "undetermined"  # sig, polarity unknown
    assert classify_direction(None, 0.3) == "undetermined"
    assert classify_direction(0.01, float("nan")) == "undetermined"


def _aug_df() -> pd.DataFrame:
    return pd.DataFrame({
        "pas_id": ["1", "2"],
        "gene_id": ["ENSG1", "ENSG2"],
        "chrom": ["chr1", "chr2"],
        "start": ["100", "200"],
        "end": ["110", "260"],
        "strand": ["+", "-"],
        "cluster1": ["A", "A"],
        "cluster2": ["B", "B"],
        "pvalue": [0.001, 0.9],
        "qvalue": [0.01, 0.95],
        "delta_proportion": [0.4, 0.01],
        "odds_ratio": [3.0, 1.0],
        "n_cells": [100, 100],
        "n_cells_cluster1": [60, 60],
        "n_cells_cluster2": [40, 40],
        "n_reads_pas_cluster1": [30, 5],
        "n_reads_pas_cluster2": [10, 4],
    })


def test_findings_long_shape_keys_and_direction() -> None:
    out = findings_long({("A", "B"): _aug_df()}, strategy="fisher", arm="switch_diff:fisher")
    assert list(out.columns) == FINDING_LONG_COLUMNS
    assert len(out) == 2

    r0 = out.iloc[0]
    # pas_uid content-addressed, strand-aware 3'
    assert r0["pas_uid"] == "chr1:109:+"
    assert r0["gene_id"] == "ENSG1"
    assert r0["canonical_cluster"] == "A"
    assert r0["comparison_cluster"] == "B"
    assert r0["strategy"] == "fisher"
    # single-PAS genes here → significance-based fallback basis
    assert r0["direction_basis"] == "significance"
    assert r0["direction"] == "undetermined"  # significant, single-PAS gene
    assert r0["n_reads"] == 40  # 30 + 10
    assert r0["n_cells_subject"] == 60
    assert r0["finding_uid"].startswith("switch_diff_fisher:fisher:A:1:chr1:109:+")

    r1 = out.iloc[1]
    assert r1["pas_uid"] == "chr2:200:-"  # - strand summit = start
    assert r1["direction"] == "flat"  # q=0.95 not significant (single-PAS fallback)
    # every direction is one of the exhaustive enum values, never NA
    assert set(out["direction"]) <= {"shorten", "lengthen", "flat", "undetermined"}


def test_canonical_map_applied() -> None:
    out = findings_long(
        {("A", "B"): _aug_df()}, strategy="fisher", arm="arm1",
        canonical_map={"A": "canon3", "B": "canon5"},
    )
    assert set(out["canonical_cluster"]) == {"canon3"}
    assert set(out["comparison_cluster"]) == {"canon5"}


def test_write_long_table_falls_back_to_tsv_without_parquet(tmp_path: Path) -> None:
    out = findings_long({("A", "B"): _aug_df()}, strategy="fisher", arm="arm1")
    written = write_long_table(out, str(tmp_path / "switch_diff_long"))
    # parquet engine may be absent in this env → .tsv fallback; either is fine.
    assert written.endswith((".parquet", ".tsv"))
    assert Path(written).exists()
    if written.endswith(".tsv"):
        back = pd.read_csv(written, sep="\t")
        assert list(back.columns) == FINDING_LONG_COLUMNS


# --- E5 length_long (LengthRow) ---
from ema.switch_test.long_output import LENGTH_LONG_COLUMNS, cell_uid, length_long


def test_cell_uid_namespacing() -> None:
    assert cell_uid("dsA", "AAAA") == "dsA:AAAA"
    # already namespaced -> unchanged
    assert cell_uid("dsA", "dsB:CCCC") == "dsB:CCCC"


def test_length_long_proportion_with_pas_uid_and_rank() -> None:
    df = pd.DataFrame({
        "gene_id": ["ENSG1", "ENSG1"],
        "transcript_id": ["ENST1", "_gene_"],
        "pas_id": ["1", "2"],
        "rank": [0, 1],
        "cell": ["AAAA", "CCCC"],
        "proportion": [0.7, 0.3],
        "cluster": ["A", "A"],
    })
    out = length_long(
        df, strategy="proportion", value_col="proportion", dataset_id="dsA",
        pas_col="pas_id", rank_col="rank", pas_uid_map={"1": "chr1:9:+", "2": "chr1:50:+"},
    )
    assert list(out.columns) == LENGTH_LONG_COLUMNS
    assert list(out["cell_uid"]) == ["dsA:AAAA", "dsA:CCCC"]
    assert list(out["value"]) == [0.7, 0.3]
    assert list(out["pas_uid"]) == ["chr1:9:+", "chr1:50:+"]
    assert list(out["rank"]) == [0, 1]
    # '_gene_' sentinel normalized to a missing value. pandas>=3 stores object
    # missings as NaN (not None), so assert missingness with pd.isna — the check
    # every consumer (and parquet round-trip) actually uses — not identity.
    tids = out["transcript_id"].tolist()
    assert tids[0] == "ENST1"
    assert pd.isna(tids[1])


def test_length_long_shannon_gene_level() -> None:
    df = pd.DataFrame({"gene_id": ["ENSG1"], "cell": ["AAAA"], "entropy": [1.2], "cluster": ["B"]})
    out = length_long(df, strategy="shannon", value_col="entropy", dataset_id="dsA")
    assert out.iloc[0]["value"] == 1.2
    assert pd.isna(out.iloc[0]["pas_uid"])
    # shannon (entropy) has no proximal/distal polarity → undetermined, never NA.
    assert out.iloc[0]["direction"] == "undetermined"
    assert out.iloc[0]["direction_basis"] == "structural"
    assert out.iloc[0]["canonical_cluster"] == "B"


# --- D8 deterministic structural direction ---
from ema.switch_test.long_output import structural_direction_by_gene


def _gene_rows(strand, positions_deltas, qval=0.001):
    """positions_deltas: list of (start,end,delta_proportion) for one gene."""
    rows = []
    for i, (s, e, d) in enumerate(positions_deltas):
        rows.append({
            "gene_id": "G", "chrom": "chr1", "start": str(s), "end": str(e),
            "strand": strand, "delta_proportion": d, "qvalue": qval,
        })
    return rows


def test_structural_plus_strand_distal_up_is_lengthen() -> None:
    # + strand: distal = largest coordinate. Give the distal PAS a positive delta.
    rows = _gene_rows("+", [(100, 110, -0.3), (500, 510, +0.4)])
    assert structural_direction_by_gene(rows)["G"] == "lengthen"


def test_structural_plus_strand_distal_down_is_shorten() -> None:
    rows = _gene_rows("+", [(100, 110, +0.3), (500, 510, -0.4)])
    assert structural_direction_by_gene(rows)["G"] == "shorten"


def test_structural_minus_strand_distal_is_smallest_coord() -> None:
    # - strand: distal = SMALLEST coordinate (farthest in transcription dir).
    # distal PAS at start=100 with +delta → distal usage up → lengthen.
    rows = _gene_rows("-", [(100, 110, +0.4), (500, 510, -0.3)])
    assert structural_direction_by_gene(rows)["G"] == "lengthen"
    # flip the distal delta sign → shorten
    rows2 = _gene_rows("-", [(100, 110, -0.4), (500, 510, +0.3)])
    assert structural_direction_by_gene(rows2)["G"] == "shorten"


def test_structural_not_significant_is_flat() -> None:
    rows = _gene_rows("+", [(100, 110, -0.3), (500, 510, +0.4)], qval=0.9)
    assert structural_direction_by_gene(rows)["G"] == "flat"


def test_structural_single_pas_gene_omitted() -> None:
    rows = _gene_rows("+", [(100, 110, +0.4)])
    assert "G" not in structural_direction_by_gene(rows)  # no structural call


def test_structural_alternative_last_exon_ranks_by_3prime() -> None:
    # Three PAS (tandem + ALE); distal is the farthest 3' one (end 900).
    rows = _gene_rows("+", [(100, 110, +0.1), (500, 510, +0.1), (900, 910, +0.5)])
    assert structural_direction_by_gene(rows)["G"] == "lengthen"


def test_findings_long_two_pas_gene_gets_structural_basis() -> None:
    df = pd.DataFrame({
        "pas_id": ["1", "2"],
        "gene_id": ["G", "G"],
        "chrom": ["chr1", "chr1"],
        "start": ["100", "500"],
        "end": ["110", "510"],
        "strand": ["+", "+"],
        "pvalue": [0.001, 0.001],
        "qvalue": [0.001, 0.001],
        "delta_proportion": [-0.4, 0.4],  # distal (end 510) up → lengthen
        "cluster1": ["A", "A"], "cluster2": ["B", "B"],
    })
    out = findings_long({("A", "B"): df}, strategy="fisher", arm="arm1")
    assert set(out["direction"]) == {"lengthen"}
    assert set(out["direction_basis"]) == {"structural"}


# --- D8 length polarity: structural_length_direction (geneview-critical) ---
from ema.switch_test.long_output import structural_length_direction


def test_length_dir_classic_pdui_lengthen_and_shorten():
    # PDUI = distal fraction. Cluster A high PDUI vs B low → A lengthens, B shortens
    # (one-vs-rest Δdistal-usage sign). Strand handled upstream in PDUI.
    df = pd.DataFrame({
        "gene_id": ["G", "G"],
        "cluster": ["A", "B"],
        "pdui": [0.8, 0.2],
        "cell": ["c1", "c2"],
    })
    dmap = structural_length_direction(df, strategy="classic", value_col="pdui")
    assert dmap[("G", "A")] == "lengthen"
    assert dmap[("G", "B")] == "shorten"


def test_length_dir_classic_flat_when_equal():
    df = pd.DataFrame({
        "gene_id": ["G", "G"], "cluster": ["A", "B"],
        "pdui": [0.5, 0.5], "cell": ["c1", "c2"],
    })
    dmap = structural_length_direction(df, strategy="classic", value_col="pdui")
    assert dmap[("G", "A")] == "flat" and dmap[("G", "B")] == "flat"


def test_length_dir_single_cluster_undetermined():
    df = pd.DataFrame({"gene_id": ["G", "G"], "cluster": ["A", "A"],
                       "pdui": [0.8, 0.6], "cell": ["c1", "c2"]})
    dmap = structural_length_direction(df, strategy="classic", value_col="pdui")
    assert dmap[("G", "A")] == "undetermined"  # no >=2-cluster contrast


def test_length_dir_proportion_uses_distal_max_rank():
    # Tandem 3'UTR: rank 0 = proximal, rank 1 = distal. Distal usage higher in A.
    df = pd.DataFrame({
        "gene_id": ["G", "G", "G", "G"],
        "cluster": ["A", "A", "B", "B"],
        "pas_id": ["p0", "p1", "p0", "p1"],
        "rank": [0, 1, 0, 1],
        "proportion": [0.3, 0.7, 0.8, 0.2],   # distal(rank1): A=0.7 vs B=0.2
        "cell": ["c1", "c1", "c2", "c2"],
    })
    dmap = structural_length_direction(
        df, strategy="proportion", value_col="proportion", rank_col="rank")
    assert dmap[("G", "A")] == "lengthen"   # more distal usage
    assert dmap[("G", "B")] == "shorten"


def test_length_dir_proportion_needs_rank():
    df = pd.DataFrame({"gene_id": ["G", "G"], "cluster": ["A", "B"],
                       "proportion": [0.7, 0.2], "cell": ["c1", "c2"]})
    # no rank_col → cannot identify distal PAS → empty map (rows → undetermined)
    assert structural_length_direction(df, strategy="proportion",
                                       value_col="proportion") == {}


def test_length_dir_shannon_has_no_polarity():
    df = pd.DataFrame({"gene_id": ["G", "G"], "cluster": ["A", "B"],
                       "entropy": [1.2, 0.4], "cell": ["c1", "c2"]})
    assert structural_length_direction(df, strategy="shannon", value_col="entropy") == {}


def test_length_long_classic_stamps_structural_direction():
    # End-to-end through length_long: classic PDUI, two clusters → real polarity.
    df = pd.DataFrame({
        "gene_id": ["G", "G"],
        "transcript_id": ["_gene_", "_gene_"],
        "cluster": ["A", "B"],
        "pdui": [0.9, 0.1],
        "cell": ["c1", "c2"],
    })
    out = length_long(df, strategy="classic", value_col="pdui", dataset_id="dsA")
    by_clu = dict(zip(out["canonical_cluster"], out["direction"]))
    assert by_clu["A"] == "lengthen" and by_clu["B"] == "shorten"
    assert set(out["direction_basis"]) == {"structural"}
    # direction is always a valid contract enum, never NA
    assert set(out["direction"]) <= {"shorten", "lengthen", "flat", "undetermined"}


def test_length_long_direction_never_na_any_strategy():
    for strategy, vcol, extra in [
        ("classic", "pdui", {}),
        ("shannon", "entropy", {}),
    ]:
        df = pd.DataFrame({"gene_id": ["G"], "cluster": ["A"],
                           vcol: [0.5], "cell": ["c1"]})
        out = length_long(df, strategy=strategy, value_col=vcol, dataset_id="d")
        assert out.iloc[0]["direction"] in {"shorten", "lengthen", "flat", "undetermined"}


# ---------------------------------------------------------------------------
# A WITHHELD q-value must never be upgraded into a confident call.
#
# nb_pairwise withholds the q-value of a dispersion-floored test (issue #94)
# precisely because that test is not trustworthy. The structural classifier
# guarded "flat" with `q is not None and q >= fdr`, so a withheld q skipped the
# guard entirely and fell through to lengthen/shorten -- suppressing the
# q-value made the row look MORE biologically confident than an ordinary
# non-significant one, and `findings_long` writes that into the
# switch_diff_long parquet the hub consumes.
# ---------------------------------------------------------------------------

def test_withheld_qvalue_is_not_promoted_to_a_directional_call() -> None:
    """Identical geometry; only the q-value differs."""
    import math

    not_significant = _gene_rows("+", [(100, 110, -0.3), (500, 510, +0.4)], qval=0.9)
    withheld = _gene_rows("+", [(100, 110, -0.3), (500, 510, +0.4)],
                          qval=float("nan"))

    assert structural_direction_by_gene(not_significant)["G"] == "flat"
    got = structural_direction_by_gene(withheld)["G"]
    assert got == "undetermined", (
        "a withheld (NaN) q-value produced the directional call "
        f"{got!r}. Withholding a q-value marks a test as untrustworthy; it "
        "must never make the row look more confident than a plainly "
        "non-significant one."
    )
    assert math.isnan(float("nan"))  # guard: the fixture really is NaN


def test_missing_qvalue_column_is_also_undetermined() -> None:
    """A row with no q at all is not evidence of a direction either."""
    rows = [
        {"gene_id": "G", "chrom": "chr1", "start": "100", "end": "110",
         "strand": "+", "delta_proportion": -0.3},
        {"gene_id": "G", "chrom": "chr1", "start": "500", "end": "510",
         "strand": "+", "delta_proportion": +0.4},
    ]
    assert structural_direction_by_gene(rows)["G"] == "undetermined"
