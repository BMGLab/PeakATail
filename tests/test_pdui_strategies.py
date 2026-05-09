"""Tests for ema.quantification.strategies.

Covers:
- Registry: get_pdui_strategy, list_pdui_strategies
- ClassicPDUIStrategy:
    - Output schema [gene_id, transcript_id, cell, pdui]
    - PDUI = distal / (prox + distal) matches hand-computed value
    - Gene with <2 PAS excluded
    - per_gene vs per_isoform aggregation
- ProportionPDUIStrategy:
    - Output schema [gene_id, transcript_id, pas_id, rank, cell, proportion]
    - Proportions sum to 1.0 per (gene/transcript, cell) within 1e-9
    - Zero-count cell → NaN proportions
- ShannonPDUIStrategy:
    - Output schema [gene_id, transcript_id, cell, entropy, normalized_entropy, n_pas]
    - Single-PAS gene → H = 0.0, H_norm = NaN
    - All reads on one PAS (delta) → H = 0.0, H_norm = 0.0
    - Uniform 2-PAS → H = 1.0, H_norm = 1.0
    - Uniform N-PAS → H = log2(N), H_norm = 1.0
    - Zero-count cell → H = NaN, H_norm = NaN
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from ema.quantification.strategies import (
    get_pdui_strategy,
    list_pdui_strategies,
    register_pdui_strategy,
)
from ema.quantification.strategies.base import PDUIStrategy


# ---------------------------------------------------------------------------
# Synthetic fixtures
# ---------------------------------------------------------------------------

def make_count_matrix(
    pas_ids: list[int],
    cells: list[str],
    data: list[list[float]],
) -> pd.DataFrame:
    """Build a count DataFrame (n_pas x n_cells)."""
    return pd.DataFrame(data, index=pas_ids, columns=cells)


# 3 genes with distinct PAS configurations
#
# GENE_A: 2 PAS — PAS 1 (proximal, rank 0) + PAS 2 (distal, rank 1)
# GENE_B: 3 PAS — PAS 3 (rank 0) + PAS 4 (rank 1) + PAS 5 (rank 2)
# GENE_C: 1 PAS — PAS 6 (rank 0 only — should be excluded from classic)

CELLS = ["cell1", "cell2", "cell3"]

# Count matrix: PAS x cells
# PAS_ID | cell1 | cell2 | cell3
#    1   |  10   |   0   |   5
#    2   |   0   |   4   |   5
#    3   |   6   |   0   |   2
#    4   |   2   |   0   |   2
#    5   |   2   |   0   |   2
#    6   |   8   |   0   |   3

COUNTS = pd.DataFrame(
    {
        "cell1": [10, 0, 6, 2, 2, 8],
        "cell2": [ 0, 4, 0, 0, 0, 0],
        "cell3": [ 5, 5, 2, 2, 2, 3],
    },
    index=[1, 2, 3, 4, 5, 6],
)


def make_pas_isoform_map() -> dict:
    """Build a PAS -> isoform map for the synthetic 3-gene setup.

    Per-isoform mapping (single transcript per gene for simplicity):
        GENE_A / T_A : PAS 1 (rank 0), PAS 2 (rank 1) — total 2
        GENE_B / T_B : PAS 3 (rank 0), PAS 4 (rank 1), PAS 5 (rank 2) — total 3
        GENE_C / T_C : PAS 6 (rank 0) — total 1
    """
    return {
        1: [("GENE_A", "T_A", 100, 0, 2)],
        2: [("GENE_A", "T_A", 300, 1, 2)],
        3: [("GENE_B", "T_B", 50,  0, 3)],
        4: [("GENE_B", "T_B", 150, 1, 3)],
        5: [("GENE_B", "T_B", 250, 2, 3)],
        6: [("GENE_C", "T_C", 100, 0, 1)],
    }


PAS_MAP = make_pas_isoform_map()


# ---------------------------------------------------------------------------
# Registry tests
# ---------------------------------------------------------------------------

class TestRegistry:
    def test_list_returns_list(self):
        names = list_pdui_strategies()
        assert isinstance(names, list)

    def test_three_strategies_registered(self):
        names = list_pdui_strategies()
        assert set(names) >= {"classic", "proportion", "shannon"}

    def test_get_returns_strategy_instance(self):
        for name in ("classic", "proportion", "shannon"):
            s = get_pdui_strategy(name)
            assert isinstance(s, PDUIStrategy)

    def test_get_unknown_raises(self):
        with pytest.raises(ValueError, match="Unknown PDUI strategy"):
            get_pdui_strategy("nonexistent_strategy_xyz")

    def test_each_call_returns_new_instance(self):
        s1 = get_pdui_strategy("classic")
        s2 = get_pdui_strategy("classic")
        assert s1 is not s2

    def test_register_requires_name(self):
        with pytest.raises(AttributeError):
            @register_pdui_strategy
            class _Bad(PDUIStrategy):
                name = ""  # empty name should raise
                def compute(self, *a, **kw): ...


# ---------------------------------------------------------------------------
# ClassicPDUIStrategy tests
# ---------------------------------------------------------------------------

class TestClassicStrategy:
    @pytest.fixture
    def strategy(self):
        return get_pdui_strategy("classic")

    def test_output_schema(self, strategy):
        df = strategy.compute(COUNTS, PAS_MAP)
        assert set(df.columns) >= {"gene_id", "transcript_id", "cell", "pdui"}

    def test_returns_dataframe(self, strategy):
        df = strategy.compute(COUNTS, PAS_MAP)
        assert isinstance(df, pd.DataFrame)

    def test_gene_c_excluded_from_per_isoform(self, strategy):
        """GENE_C has 1 PAS → classic excludes it (needs >=2)."""
        df = strategy.compute(COUNTS, PAS_MAP, aggregation="per_isoform")
        genes_in_output = df["gene_id"].unique().tolist()
        assert "GENE_C" not in genes_in_output

    def test_gene_c_excluded_from_per_gene(self, strategy):
        df = strategy.compute(COUNTS, PAS_MAP, aggregation="per_gene")
        assert "GENE_C" not in df["gene_id"].unique()

    def test_gene_a_hand_computed_per_isoform(self, strategy):
        """GENE_A: proximal=PAS1, distal=PAS2.

        cell1: prox=10, dist=0  → PDUI = 0/(10+0) = 0.0
        cell2: prox=0,  dist=4  → PDUI = 4/(0+4)  = 1.0
        cell3: prox=5,  dist=5  → PDUI = 5/(5+5)  = 0.5
        """
        df = strategy.compute(COUNTS, PAS_MAP, aggregation="per_isoform")
        gene_a = df[df["gene_id"] == "GENE_A"].set_index("cell")["pdui"]
        assert abs(gene_a.loc["cell1"] - 0.0) < 1e-9
        assert abs(gene_a.loc["cell2"] - 1.0) < 1e-9
        assert abs(gene_a.loc["cell3"] - 0.5) < 1e-9

    def test_zero_total_cell_is_nan(self, strategy):
        """When proximal + distal = 0 and pseudocount=0 → NaN."""
        # GENE_B, cell2: all PAS have 0 counts
        df = strategy.compute(COUNTS, PAS_MAP, aggregation="per_isoform", pseudocount=0.0)
        gene_b_cell2 = df[(df["gene_id"] == "GENE_B") & (df["cell"] == "cell2")]["pdui"]
        assert gene_b_cell2.isna().all()

    def test_pseudocount_avoids_nan(self, strategy):
        """With pseudocount > 0, PDUI should not be NaN even for zero-count cells."""
        df = strategy.compute(COUNTS, PAS_MAP, aggregation="per_isoform", pseudocount=1.0)
        gene_b_cell2 = df[(df["gene_id"] == "GENE_B") & (df["cell"] == "cell2")]["pdui"]
        assert not gene_b_cell2.isna().all()

    def test_per_gene_transcript_id_sentinel(self, strategy):
        """per_gene mode uses '_gene_' as transcript_id."""
        df = strategy.compute(COUNTS, PAS_MAP, aggregation="per_gene")
        assert (df["transcript_id"] == "_gene_").all()

    def test_per_isoform_has_real_transcript_ids(self, strategy):
        """per_isoform mode uses real transcript IDs."""
        df = strategy.compute(COUNTS, PAS_MAP, aggregation="per_isoform")
        assert "_gene_" not in df["transcript_id"].values

    def test_empty_map_returns_empty_df(self, strategy):
        df = strategy.compute(COUNTS, {})
        assert df.empty
        assert "pdui" in df.columns

    def test_pdui_range_zero_to_one(self, strategy):
        df = strategy.compute(COUNTS, PAS_MAP, aggregation="per_isoform")
        non_nan = df["pdui"].dropna()
        assert (non_nan >= 0.0).all()
        assert (non_nan <= 1.0).all()


# ---------------------------------------------------------------------------
# ProportionPDUIStrategy tests
# ---------------------------------------------------------------------------

class TestProportionStrategy:
    @pytest.fixture
    def strategy(self):
        return get_pdui_strategy("proportion")

    def test_output_schema(self, strategy):
        df = strategy.compute(COUNTS, PAS_MAP)
        for col in ("gene_id", "transcript_id", "pas_id", "rank", "cell", "proportion"):
            assert col in df.columns, f"Missing column: {col}"

    def test_proportions_sum_to_one_per_gene_per_cell(self, strategy):
        """For per_gene aggregation, proportions per (gene, cell) must sum to 1.0.

        Cells where all proportions are NaN (zero-count) are excluded from the
        sum-to-1 check.  pandas groupby sum() treats NaN as 0 by default, so
        we use min_count=1 to propagate NaN, then drop.
        """
        df = strategy.compute(COUNTS, PAS_MAP, aggregation="per_gene")
        grouped = df.groupby(["gene_id", "cell"])["proportion"].sum(min_count=1)
        non_nan = grouped.dropna()
        assert len(non_nan) > 0
        assert (abs(non_nan - 1.0) < 1e-9).all(), \
            f"Proportions don't sum to 1.0: {non_nan[abs(non_nan - 1.0) >= 1e-9]}"

    def test_proportions_sum_to_one_per_isoform_per_cell(self, strategy):
        """For per_isoform, proportions per (transcript, cell) must sum to 1.0.

        Cells where all proportions are NaN (zero-count) are excluded via
        min_count=1 to propagate NaN through the sum.
        """
        df = strategy.compute(COUNTS, PAS_MAP, aggregation="per_isoform")
        grouped = df.groupby(["transcript_id", "cell"])["proportion"].sum(min_count=1)
        non_nan = grouped.dropna()
        assert (abs(non_nan - 1.0) < 1e-9).all(), \
            f"Per-isoform proportions don't sum to 1.0"

    def test_zero_total_cell_nan_proportion(self, strategy):
        """cell2 has 0 counts for GENE_B → proportions should be NaN."""
        df = strategy.compute(COUNTS, PAS_MAP, aggregation="per_isoform")
        gene_b_c2 = df[(df["gene_id"] == "GENE_B") & (df["cell"] == "cell2")]["proportion"]
        assert gene_b_c2.isna().all()

    def test_proportion_dtype(self, strategy):
        df = strategy.compute(COUNTS, PAS_MAP)
        assert df["proportion"].dtype == np.float64

    def test_gene_a_cell1_proportions(self, strategy):
        """GENE_A, cell1: PAS1=10, PAS2=0 → props = [1.0, 0.0]."""
        df = strategy.compute(COUNTS, PAS_MAP, aggregation="per_gene")
        g = df[(df["gene_id"] == "GENE_A") & (df["cell"] == "cell1")].set_index("pas_id")
        assert abs(g.loc[1, "proportion"] - 1.0) < 1e-9
        assert abs(g.loc[2, "proportion"] - 0.0) < 1e-9

    def test_gene_a_cell3_uniform_proportions(self, strategy):
        """GENE_A, cell3: PAS1=5, PAS2=5 → both proportions = 0.5."""
        df = strategy.compute(COUNTS, PAS_MAP, aggregation="per_gene")
        g = df[(df["gene_id"] == "GENE_A") & (df["cell"] == "cell3")].set_index("pas_id")
        assert abs(g.loc[1, "proportion"] - 0.5) < 1e-9
        assert abs(g.loc[2, "proportion"] - 0.5) < 1e-9

    def test_empty_map_returns_empty_df(self, strategy):
        df = strategy.compute(COUNTS, {})
        assert df.empty
        assert "proportion" in df.columns


# ---------------------------------------------------------------------------
# ShannonPDUIStrategy tests
# ---------------------------------------------------------------------------

class TestShannonStrategy:
    @pytest.fixture
    def strategy(self):
        return get_pdui_strategy("shannon")

    def test_output_schema(self, strategy):
        df = strategy.compute(COUNTS, PAS_MAP)
        for col in ("gene_id", "transcript_id", "cell",
                    "entropy", "normalized_entropy", "n_pas"):
            assert col in df.columns, f"Missing column: {col}"

    # --- Boundary cases ---

    def test_single_pas_gene_entropy_zero(self, strategy):
        """GENE_C has 1 PAS → H = 0.0 (all probability on one event)."""
        df = strategy.compute(COUNTS, PAS_MAP, aggregation="per_isoform")
        gene_c = df[(df["gene_id"] == "GENE_C") & (df["cell"] == "cell1")]
        assert not gene_c.empty
        assert abs(gene_c["entropy"].iloc[0] - 0.0) < 1e-9

    def test_single_pas_normalized_entropy_nan(self, strategy):
        """GENE_C has 1 PAS → H_norm = NaN (log2(1) = 0, division undefined)."""
        df = strategy.compute(COUNTS, PAS_MAP, aggregation="per_isoform")
        gene_c = df[(df["gene_id"] == "GENE_C") & (df["cell"] == "cell1")]
        assert gene_c["normalized_entropy"].isna().iloc[0]

    def test_delta_distribution_entropy_zero(self, strategy):
        """GENE_A, cell1: PAS1=10, PAS2=0 → H=0 (delta on PAS1)."""
        df = strategy.compute(COUNTS, PAS_MAP, aggregation="per_isoform")
        gene_a_c1 = df[(df["gene_id"] == "GENE_A") & (df["cell"] == "cell1")]
        assert abs(gene_a_c1["entropy"].iloc[0] - 0.0) < 1e-9
        assert abs(gene_a_c1["normalized_entropy"].iloc[0] - 0.0) < 1e-9

    def test_uniform_2_pas_entropy_equals_1(self, strategy):
        """GENE_A, cell3: PAS1=5, PAS2=5 → H=1.0, H_norm=1.0."""
        df = strategy.compute(COUNTS, PAS_MAP, aggregation="per_isoform")
        gene_a_c3 = df[(df["gene_id"] == "GENE_A") & (df["cell"] == "cell3")]
        assert abs(gene_a_c3["entropy"].iloc[0] - 1.0) < 1e-9
        assert abs(gene_a_c3["normalized_entropy"].iloc[0] - 1.0) < 1e-9

    def test_uniform_n_pas_entropy_equals_log2n(self, strategy):
        """Uniform 3-PAS gene → H = log2(3), H_norm = 1.0."""
        # Build a 3-PAS count matrix with equal counts
        counts = pd.DataFrame(
            {"cell1": [4, 4, 4]},
            index=[10, 11, 12],
        )
        pas_map = {
            10: [("GENEX", "TX", 10, 0, 3)],
            11: [("GENEX", "TX", 20, 1, 3)],
            12: [("GENEX", "TX", 30, 2, 3)],
        }
        df = strategy.compute(counts, pas_map, aggregation="per_isoform")
        row = df[df["gene_id"] == "GENEX"].iloc[0]
        expected_h = math.log2(3)
        assert abs(row["entropy"] - expected_h) < 1e-9
        assert abs(row["normalized_entropy"] - 1.0) < 1e-9

    def test_zero_total_cell_nan_entropy(self, strategy):
        """Zero-count cell → H = NaN."""
        df = strategy.compute(COUNTS, PAS_MAP, aggregation="per_isoform")
        gene_b_c2 = df[(df["gene_id"] == "GENE_B") & (df["cell"] == "cell2")]
        assert gene_b_c2["entropy"].isna().all()
        assert gene_b_c2["normalized_entropy"].isna().all()

    def test_entropy_in_valid_range(self, strategy):
        """H must be in [0, log2(N)] for all non-NaN rows."""
        df = strategy.compute(COUNTS, PAS_MAP, aggregation="per_isoform")
        non_nan = df.dropna(subset=["entropy"])
        assert (non_nan["entropy"] >= 0.0).all()
        max_h = non_nan.apply(
            lambda r: math.log2(r["n_pas"]) if r["n_pas"] > 1 else 0.0,
            axis=1,
        )
        assert (non_nan["entropy"] <= max_h + 1e-9).all()

    def test_normalized_entropy_in_zero_one(self, strategy):
        """H_norm must be in [0, 1] for all non-NaN rows."""
        df = strategy.compute(COUNTS, PAS_MAP, aggregation="per_isoform")
        non_nan = df["normalized_entropy"].dropna()
        assert (non_nan >= 0.0 - 1e-9).all()
        assert (non_nan <= 1.0 + 1e-9).all()

    def test_n_pas_column_correct(self, strategy):
        """n_pas should equal the number of PAS in the transcript."""
        df = strategy.compute(COUNTS, PAS_MAP, aggregation="per_isoform")
        gene_a = df[df["gene_id"] == "GENE_A"]
        assert (gene_a["n_pas"] == 2).all()
        gene_b = df[df["gene_id"] == "GENE_B"]
        assert (gene_b["n_pas"] == 3).all()
        gene_c = df[df["gene_id"] == "GENE_C"]
        assert (gene_c["n_pas"] == 1).all()

    def test_per_gene_vs_per_isoform_same_when_one_transcript(self, strategy):
        """When each gene has one transcript, per_gene and per_isoform H are equal."""
        df_iso = strategy.compute(COUNTS, PAS_MAP, aggregation="per_isoform")
        df_gene = strategy.compute(COUNTS, PAS_MAP, aggregation="per_gene")

        for gene_id in ("GENE_A", "GENE_B", "GENE_C"):
            for cell in CELLS:
                h_iso = df_iso[
                    (df_iso["gene_id"] == gene_id) & (df_iso["cell"] == cell)
                ]["entropy"].values
                h_gene = df_gene[
                    (df_gene["gene_id"] == gene_id) & (df_gene["cell"] == cell)
                ]["entropy"].values
                if len(h_iso) and len(h_gene):
                    both_nan = (
                        (np.isnan(h_iso[0]) and np.isnan(h_gene[0]))
                    )
                    if not both_nan:
                        assert abs(float(h_iso[0]) - float(h_gene[0])) < 1e-9, \
                            f"Entropy mismatch {gene_id}/{cell}"

    def test_empty_map_returns_empty_df(self, strategy):
        df = strategy.compute(COUNTS, {})
        assert df.empty
        assert "entropy" in df.columns
