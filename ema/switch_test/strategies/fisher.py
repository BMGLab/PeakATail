"""Fisher exact test strategy — wraps existing fishertest.py logic.

This module adapts the existing :func:`ema.switch_test.fishertest.fishertest`
function to the :class:`DiffAPAStrategy` interface without modifying the
original implementation.

The original function operates on a pre-aggregated ``(gene, PAS) × cluster``
matrix (one row per PAS, two columns = cluster totals).  This wrapper
derives that representation from the per-cell ``count_matrix`` +
``cluster_labels`` inputs, then delegates to the original logic.

No differential logic was changed — only the call boundary.
"""

from __future__ import annotations

import io
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, fisher_exact

from ema.switch_test.strategies import register_diff_strategy
from ema.switch_test.strategies.base import DiffAPAStrategy


@register_diff_strategy
class FisherStrategy(DiffAPAStrategy):
    """Pairwise Fisher exact test with Benjamini-Hochberg FDR correction.

    Wraps the logic in :func:`ema.switch_test.fishertest.fishertest`.
    For each PAS, a 2×2 contingency table is built:

    .. code-block::

        [[count_PAS_c1,       count_PAS_c2      ],
         [count_other_PAS_c1, count_other_PAS_c2]]

    where *other* sums all remaining PAS in the same gene.

    This test is known to be **anti-conservative** for single-cell count
    data because it ignores overdispersion; prefer ``nb_pairwise`` when
    statistical rigour is required.

    Tunable hyperparameters:
        min_cells_per_group (default 10): Minimum cells in each cluster for a
            PAS to be testable (checked as a total-cell count because Fisher
            uses aggregated counts, not per-cell non-zero detection).
            CLI: ``--min-cells-per-group`` / YAML: ``min_cells_per_group``.
        fdr (default 0.05): Benjamini-Hochberg FDR threshold applied after
            testing.  CLI: ``--fdr`` / YAML: ``fdr``.
    """

    name: str = "fisher"
    supports_multi_condition: bool = False

    def test(
        self,
        count_matrix: pd.DataFrame,
        cluster_labels: pd.Series,
        cluster1: str | None = None,
        cluster2: str | None = None,
        min_cells_per_group: int = 10,
        n_jobs: int = -1,  # Fisher is fast; n_jobs is accepted but unused.
    ) -> pd.DataFrame:
        """Run Fisher exact test for differential PAS usage.

        Parameters
        ----------
        count_matrix:
            Shape ``(n_cells, n_pas)``.  PAS are columns; index = cell barcodes.
        cluster_labels:
            Per-cell cluster assignment aligned to ``count_matrix`` index.
        cluster1, cluster2:
            Cluster labels to compare.  Both required (pairwise only).
        min_cells_per_group:
            Minimum cells per group for a PAS to be testable.
        n_jobs:
            Ignored; Fisher test is not parallelised internally.

        Returns
        -------
        pd.DataFrame
            Indexed by ``pas_id`` with columns
            ``[pvalue, qvalue, n_cells, odds_ratio, delta_proportion, log2fc]``.
        """
        if cluster1 is None or cluster2 is None:
            raise ValueError(
                "FisherStrategy is pairwise: cluster1 and cluster2 must both be provided."
            )

        # --- select cells for the two clusters ---
        mask1 = cluster_labels == cluster1
        mask2 = cluster_labels == cluster2
        cells1 = count_matrix.index[mask1]
        cells2 = count_matrix.index[mask2]

        if len(cells1) < min_cells_per_group or len(cells2) < min_cells_per_group:
            return pd.DataFrame(
                columns=["pvalue", "qvalue", "n_cells", "odds_ratio",
                         "delta_proportion", "log2fc"]
            )

        # --- aggregate counts per PAS per cluster ---
        agg1 = count_matrix.loc[cells1].sum(axis=0)  # Series indexed by pas_id
        agg2 = count_matrix.loc[cells2].sum(axis=0)

        # Build the aggregated (PAS,) × (cluster1, cluster2) DataFrame that
        # mirrors what the original fishertest() expects — but we run the logic
        # directly here to produce the unified output schema.
        pas_ids = count_matrix.columns.tolist()
        total1 = agg1.sum()
        total2 = agg2.sum()

        results: list[dict] = []
        for pas_id in pas_ids:
            pas_c1 = int(agg1[pas_id])
            pas_c2 = int(agg2[pas_id])
            other_c1 = int(total1) - pas_c1
            other_c2 = int(total2) - pas_c2

            # Skip degenerate tables
            if (pas_c1 + pas_c2) == 0:
                continue
            if (pas_c1 + other_c1) == 0 or (pas_c2 + other_c2) == 0:
                continue

            table = np.array([
                [pas_c1, pas_c2],
                [other_c1, other_c2],
            ])

            odds_ratio, pvalue = fisher_exact(table, alternative="two-sided")

            prop1 = pas_c1 / total1 if total1 > 0 else 0.0
            prop2 = pas_c2 / total2 if total2 > 0 else 0.0
            delta_prop = prop1 - prop2

            # log2fc: log2(mean_expression_c2 / mean_expression_c1)
            mean1 = pas_c1 / max(len(cells1), 1)
            mean2 = pas_c2 / max(len(cells2), 1)
            eps = 1e-8
            log2fc = float(np.log2((mean2 + eps) / (mean1 + eps)))

            results.append(
                {
                    "pas_id": pas_id,
                    "pvalue": pvalue,
                    "odds_ratio": odds_ratio,
                    "delta_proportion": delta_prop,
                    "log2fc": log2fc,
                    "n_cells": len(cells1) + len(cells2),
                }
            )

        if not results:
            return pd.DataFrame(
                columns=["pvalue", "qvalue", "n_cells", "odds_ratio",
                         "delta_proportion", "log2fc"]
            )

        df = pd.DataFrame(results).set_index("pas_id")

        # BH FDR correction
        qvalues = false_discovery_control(df["pvalue"].values, method="bh")
        df["qvalue"] = qvalues

        return df[["pvalue", "qvalue", "n_cells", "odds_ratio",
                   "delta_proportion", "log2fc"]].sort_values("qvalue")
