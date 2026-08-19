"""Tests for `ema switch trend --metric` / `--combine`.

Covers:
  (a) --metric pdui produces the exact same output as the (unset) default.
  (b) --metric proportion / derive_distal_fraction correctly derives a
      per-(gene,cell) distal-usage scalar from a per-PAS proportion table.
  (c) --metric entropy trends normalized_entropy.
  (d) --combine consensus: agreement + ambiguous-on-conflict logic.
  (e) --combine summary counts (n_shortening/n_lengthening/n_ambiguous +
      pairwise agreement rates).

See ema/switch_test/trend.py and ema/cli/switch_trend.py.
"""
from __future__ import annotations

import json

import pandas as pd
import pytest
from click.testing import CliRunner

from ema.switch_test.trend import (
    consensus_trend_by_gene,
    consensus_trend_summary,
    derive_distal_fraction,
    pdui_trend,
    pdui_trend_by_gene,
)

ORDER = ["Normal", "StageI", "StageII"]


# ---------------------------------------------------------------------------
# Fixtures / small table builders
# ---------------------------------------------------------------------------

def _pdui_table(gene: str = "G1", vals=(0.9, 0.5, 0.1)) -> pd.DataFrame:
    """PDUI shrinks across stages -> 3'UTR shortening (decreasing direction)."""
    rows = []
    for stage, base in zip(ORDER, vals):
        for d in (-0.02, 0.0, 0.02):
            rows.append({"gene_id": gene, "cluster": stage, "pdui": base + d})
    return pd.DataFrame(rows)


def _proportion_table(gene: str = "G1", distal_vals=(0.8, 0.5, 0.2)) -> pd.DataFrame:
    """Two-PAS gene: rank=1 (proximal), rank=2 (distal). Mirrors the real
    ``ema switch length -s proportion`` schema (gene_id/transcript_id/pas_id/
    rank/cell/proportion/cluster + a few extra columns the derivation ignores).
    """
    rows = []
    for stage, distal in zip(ORDER, distal_vals):
        for cell_i in range(2):
            cell = f"{gene}_{stage}_c{cell_i}"
            rows.append({
                "gene_id": gene, "transcript_id": "_gene_", "pas_id": 1, "rank": 1,
                "cell": cell, "proportion": 1.0 - distal, "reads_at_pas": 0.0,
                "total_reads_gene": 1.0, "cluster": stage,
            })
            rows.append({
                "gene_id": gene, "transcript_id": "_gene_", "pas_id": 2, "rank": 2,
                "cell": cell, "proportion": distal, "reads_at_pas": 0.0,
                "total_reads_gene": 1.0, "cluster": stage,
            })
    return pd.DataFrame(rows)


def _entropy_table(gene: str = "G1", vals=(0.1, 0.5, 0.9)) -> pd.DataFrame:
    rows = []
    for stage, ent in zip(ORDER, vals):
        for cell_i in range(2):
            rows.append({
                "gene_id": gene, "transcript_id": "_gene_", "pas_ids": "1,2",
                "cell": f"{gene}_{stage}_c{cell_i}", "entropy": ent,
                "normalized_entropy": ent, "n_pas": 2, "total_reads_gene": 1.0,
                "cluster": stage,
            })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# (a) --metric pdui == default
# ---------------------------------------------------------------------------

def test_metric_pdui_matches_default(tmp_path):
    pdui_tsv = tmp_path / "pdui.tsv"
    _pdui_table().to_csv(pdui_tsv, sep="\t", index=False)

    from ema.cli import main

    out_default = tmp_path / "out_default"
    res_default = CliRunner().invoke(main, [
        "switch", "trend", "--pdui", str(pdui_tsv),
        "--stage-order", ",".join(ORDER), "-o", str(out_default),
    ])
    assert res_default.exit_code == 0, res_default.output

    out_metric = tmp_path / "out_metric"
    res_metric = CliRunner().invoke(main, [
        "switch", "trend", "--pdui", str(pdui_tsv), "--metric", "pdui",
        "--stage-order", ",".join(ORDER), "-o", str(out_metric),
    ])
    assert res_metric.exit_code == 0, res_metric.output

    default_summary = json.loads((out_default / "length_trend.json").read_text())
    metric_summary = json.loads((out_metric / "length_trend.json").read_text())
    assert default_summary == metric_summary
    assert default_summary["direction"] == "decreasing"

    default_by_gene = (out_default / "length_trend_by_gene.tsv").read_text()
    metric_by_gene = (out_metric / "length_trend_by_gene.tsv").read_text()
    assert default_by_gene == metric_by_gene


# ---------------------------------------------------------------------------
# (b) --metric proportion / derive_distal_fraction
# ---------------------------------------------------------------------------

def test_derive_distal_fraction_picks_max_rank_proportion():
    df = _proportion_table(distal_vals=(0.8, 0.5, 0.2))
    derived = derive_distal_fraction(df)
    assert set(derived.columns) == {"gene_id", "cell", "cluster", "value"}
    # one row per (gene, cell); value == the rank=2 (distal) proportion.
    assert len(derived) == 6
    by_stage = derived.groupby("cluster")["value"].mean().round(6).to_dict()
    assert by_stage == {"Normal": 0.8, "StageI": 0.5, "StageII": 0.2}


def test_derive_distal_fraction_missing_column_raises():
    df = pd.DataFrame({"gene_id": ["G1"], "cell": ["c1"]})
    with pytest.raises(ValueError):
        derive_distal_fraction(df)


def test_metric_proportion_trend_is_decreasing(tmp_path):
    prop_tsv = tmp_path / "proportion.tsv"
    _proportion_table(distal_vals=(0.8, 0.5, 0.2)).to_csv(prop_tsv, sep="\t", index=False)

    from ema.cli import main

    out = tmp_path / "out_prop"
    res = CliRunner().invoke(main, [
        "switch", "trend", "--metric", "proportion", "--proportion-tsv", str(prop_tsv),
        "--stage-order", ",".join(ORDER), "-o", str(out),
    ])
    assert res.exit_code == 0, res.output
    summary = json.loads((out / "length_trend.json").read_text())
    assert summary["direction"] == "decreasing"
    assert summary["value_col"] == "value"
    by_gene = pd.read_csv(out / "length_trend_by_gene.tsv", sep="\t", index_col="gene_id")
    assert by_gene.loc["G1", "direction"] == "decreasing"


# ---------------------------------------------------------------------------
# (c) --metric entropy
# ---------------------------------------------------------------------------

def test_metric_entropy_trend(tmp_path):
    ent_tsv = tmp_path / "entropy_shannon.tsv"
    _entropy_table(vals=(0.1, 0.5, 0.9)).to_csv(ent_tsv, sep="\t", index=False)

    from ema.cli import main

    out = tmp_path / "out_entropy"
    res = CliRunner().invoke(main, [
        "switch", "trend", "--metric", "entropy", "--entropy-tsv", str(ent_tsv),
        "--stage-order", ",".join(ORDER), "-o", str(out),
    ])
    assert res.exit_code == 0, res.output
    summary = json.loads((out / "length_trend.json").read_text())
    assert summary["value_col"] == "normalized_entropy"
    assert summary["direction"] == "increasing"


# ---------------------------------------------------------------------------
# (d) --combine consensus: agreement + ambiguous logic (direct function level)
# ---------------------------------------------------------------------------

def _per_gene_by_metric_two_genes():
    """G1: pdui + proportion both decreasing (agree) -> shortening.
    G2: pdui decreasing, proportion INCREASING (conflict) -> ambiguous.
    Entropy is increasing for G1 and flat for G2 (excluded from the vote by
    default; still reported).
    """
    pdui_df = pd.concat([
        _pdui_table("G1", vals=(0.9, 0.5, 0.1)),
        _pdui_table("G2", vals=(0.9, 0.5, 0.1)),
    ], ignore_index=True)
    prop_raw = pd.concat([
        _proportion_table("G1", distal_vals=(0.8, 0.5, 0.2)),
        _proportion_table("G2", distal_vals=(0.2, 0.5, 0.8)),  # conflicts with pdui
    ], ignore_index=True)
    prop_derived = derive_distal_fraction(prop_raw)
    entropy_df = pd.concat([
        _entropy_table("G1", vals=(0.1, 0.5, 0.9)),
        _entropy_table("G2", vals=(0.5, 0.5, 0.5)),  # flat
    ], ignore_index=True)

    return {
        "pdui": pdui_trend_by_gene(pdui_df, ORDER, value_col="pdui"),
        "proportion": pdui_trend_by_gene(prop_derived, ORDER, value_col="value"),
        "entropy": pdui_trend_by_gene(entropy_df, ORDER, value_col="normalized_entropy"),
    }


def test_consensus_trend_agreement_and_ambiguous():
    per_gene = _per_gene_by_metric_two_genes()
    combined = consensus_trend_by_gene(
        per_gene, vote_metrics=["pdui", "proportion"], min_agree=1,
    )
    assert combined.loc["G1", "pdui_direction"] == "decreasing"
    assert combined.loc["G1", "proportion_direction"] == "decreasing"
    assert combined.loc["G1", "entropy_direction"] == "increasing"
    assert combined.loc["G1", "consensus_direction"] == "shortening"
    assert combined.loc["G1", "agree"]
    assert combined.loc["G1", "n_voting"] == 2
    assert combined.loc["G1", "n_agree"] == 2

    assert combined.loc["G2", "pdui_direction"] == "decreasing"
    assert combined.loc["G2", "proportion_direction"] == "increasing"
    assert combined.loc["G2", "consensus_direction"] == "ambiguous"
    assert not combined.loc["G2", "agree"]
    assert combined.loc["G2", "n_voting"] == 2
    assert combined.loc["G2", "n_agree"] == 1


def test_consensus_trend_include_entropy_in_vote_changes_outcome():
    per_gene = _per_gene_by_metric_two_genes()
    # With entropy voting AND requiring all 3 to agree (min_agree=3), G1's
    # entropy (increasing) breaks the pdui+proportion (decreasing) majority
    # that was enough to call "shortening" under the default 2-metric vote.
    combined = consensus_trend_by_gene(
        per_gene, vote_metrics=["pdui", "proportion", "entropy"], min_agree=3,
    )
    assert combined.loc["G1", "consensus_direction"] == "ambiguous"


# ---------------------------------------------------------------------------
# (e) --combine summary counts + pairwise agreement
# ---------------------------------------------------------------------------

def test_consensus_trend_summary_counts_and_pairwise_agreement():
    per_gene = _per_gene_by_metric_two_genes()
    combined = consensus_trend_by_gene(
        per_gene, vote_metrics=["pdui", "proportion"], min_agree=1,
    )
    summary = consensus_trend_summary(
        combined, metrics=["pdui", "proportion", "entropy"],
        vote_metrics=["pdui", "proportion"], min_agree=1,
    )
    assert summary["n_genes"] == 2
    assert summary["n_shortening"] == 1
    assert summary["n_lengthening"] == 0
    assert summary["n_ambiguous"] == 1

    pa = summary["pairwise_agreement"]
    # pdui vs proportion: both determined for both genes; only G1 matches.
    assert pa["pdui_vs_proportion"]["n_genes_compared"] == 2
    assert pa["pdui_vs_proportion"]["agreement_rate"] == pytest.approx(0.5)
    # pdui vs entropy: only G1 has both determined (G2 entropy is flat); mismatch.
    assert pa["pdui_vs_entropy"]["n_genes_compared"] == 1
    assert pa["pdui_vs_entropy"]["agreement_rate"] == pytest.approx(0.0)
    # proportion vs entropy: same as above, only G1 comparable; mismatch.
    assert pa["proportion_vs_entropy"]["n_genes_compared"] == 1
    assert pa["proportion_vs_entropy"]["agreement_rate"] == pytest.approx(0.0)


def test_combine_cli_end_to_end(tmp_path):
    pdui_tsv = tmp_path / "pdui.tsv"
    prop_tsv = tmp_path / "proportion.tsv"
    ent_tsv = tmp_path / "entropy_shannon.tsv"
    _pdui_table("G1", vals=(0.9, 0.5, 0.1)).to_csv(pdui_tsv, sep="\t", index=False)
    _proportion_table("G1", distal_vals=(0.8, 0.5, 0.2)).to_csv(prop_tsv, sep="\t", index=False)
    _entropy_table("G1", vals=(0.1, 0.5, 0.9)).to_csv(ent_tsv, sep="\t", index=False)

    from ema.cli import main

    out = tmp_path / "out_combine"
    res = CliRunner().invoke(main, [
        "switch", "trend", "--combine", "pdui,proportion,entropy",
        "--pdui", str(pdui_tsv), "--proportion-tsv", str(prop_tsv),
        "--entropy-tsv", str(ent_tsv),
        "--stage-order", ",".join(ORDER), "-o", str(out),
    ])
    assert res.exit_code == 0, res.output
    assert (out / "trend_combined.tsv").exists()
    assert (out / "trend_combined_summary.json").exists()

    combined = pd.read_csv(out / "trend_combined.tsv", sep="\t", index_col="gene_id")
    assert combined.loc["G1", "consensus_direction"] == "shortening"
    assert combined.loc["G1", "pdui_direction"] == "decreasing"
    assert combined.loc["G1", "proportion_direction"] == "decreasing"
    assert combined.loc["G1", "entropy_direction"] == "increasing"

    summary = json.loads((out / "trend_combined_summary.json").read_text())
    assert summary["n_genes"] == 1
    assert summary["n_shortening"] == 1
    assert summary["n_ambiguous"] == 0
    assert summary["vote_metrics"] == ["pdui", "proportion"]


def test_combine_requires_at_least_two_metrics(tmp_path):
    from ema.cli import main
    pdui_tsv = tmp_path / "pdui.tsv"
    _pdui_table().to_csv(pdui_tsv, sep="\t", index=False)
    res = CliRunner().invoke(main, [
        "switch", "trend", "--combine", "pdui", "--pdui", str(pdui_tsv),
        "--stage-order", ",".join(ORDER), "-o", str(tmp_path / "out"),
    ])
    assert res.exit_code != 0


def test_combine_unknown_metric_rejected(tmp_path):
    from ema.cli import main
    pdui_tsv = tmp_path / "pdui.tsv"
    _pdui_table().to_csv(pdui_tsv, sep="\t", index=False)
    res = CliRunner().invoke(main, [
        "switch", "trend", "--combine", "pdui,nonsense", "--pdui", str(pdui_tsv),
        "--stage-order", ",".join(ORDER), "-o", str(tmp_path / "out"),
    ])
    assert res.exit_code != 0
