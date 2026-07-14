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
    assert r0["direction"] == "undetermined"  # significant, polarity unknown
    assert r0["n_reads"] == 40  # 30 + 10
    assert r0["n_cells_subject"] == 60
    assert r0["finding_uid"].startswith("switch_diff_fisher:fisher:A:1:chr1:109:+")

    r1 = out.iloc[1]
    assert r1["pas_uid"] == "chr2:200:-"  # - strand summit = start
    assert r1["direction"] == "flat"  # q=0.95 not significant
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
