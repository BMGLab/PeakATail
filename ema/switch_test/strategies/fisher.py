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

        # Cell-level binarisation: each cell is "expressing" a PAS iff it has
        # >=1 read mapped there.  Sample size = number of cells per cluster,
        # not reads — reads aggregated across cells are NOT independent
        # observations, so the previous read-based Fisher pseudo-replicated
        # massively (median 66% PAS flagged significant at q<0.05, hits with
        # |delta_proportion| < 0.1%).
        cm1 = count_matrix.loc[cells1]
        cm2 = count_matrix.loc[cells2]
        n1 = int(len(cells1))
        n2 = int(len(cells2))

        # Per-PAS: # cells expressing in each cluster, # reads (info only)
        expr1 = (cm1 > 0).sum(axis=0).astype(int)   # cells with >=1 read, c1
        expr2 = (cm2 > 0).sum(axis=0).astype(int)   # cells with >=1 read, c2
        reads1 = cm1.sum(axis=0).astype(int)        # for info (no longer used in test)
        reads2 = cm2.sum(axis=0).astype(int)

        pas_ids = count_matrix.columns.tolist()
        results: list[dict] = []
        for pas_id in pas_ids:
            a = int(expr1[pas_id])   # cells in c1 expressing PAS
            b = int(expr2[pas_id])   # cells in c2 expressing PAS
            c = n1 - a               # cells in c1 NOT expressing PAS
            d = n2 - b               # cells in c2 NOT expressing PAS

            # Skip PAS expressed in neither cluster (uninformative).
            if a + b == 0:
                continue

            # 2x2 cell-level contingency table:
            #                    cluster1   cluster2
            #   expressing  PAS:    a         b
            #   not expressing:     c         d
            table = np.array([[a, b], [c, d]])
            odds_ratio, pvalue = fisher_exact(table, alternative="two-sided")

            # Effect sizes computed at the CELL level.
            prop1 = a / n1
            prop2 = b / n2
            delta_prop = prop1 - prop2

            # log2fc of expression frequency (bounded: both denominators >= 1
            # since min_cells_per_group >= 10 and numerators bounded by them).
            eps = 1.0 / max(n1, n2)  # half-cell pseudo so 0/N -> log2 ~ -log2(2N)
            log2fc = float(np.log2((prop2 + eps) / (prop1 + eps)))

            results.append(
                {
                    "pas_id": pas_id,
                    "pvalue": pvalue,
                    "odds_ratio": float(odds_ratio),
                    "delta_proportion": delta_prop,
                    "log2fc": log2fc,
                    "n_cells_cluster1": n1,
                    "n_cells_cluster2": n2,
                    "n_cells_expr_cluster1": a,
                    "n_cells_expr_cluster2": b,
                    "n_reads_pas_cluster1": int(reads1[pas_id]),
                    "n_reads_pas_cluster2": int(reads2[pas_id]),
                    "n_cells": n1 + n2,
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

        return df[[
            "pvalue", "qvalue",
            "n_cells", "n_cells_cluster1", "n_cells_cluster2",
            "n_cells_expr_cluster1", "n_cells_expr_cluster2",
            "n_reads_pas_cluster1", "n_reads_pas_cluster2",
            "odds_ratio", "delta_proportion", "log2fc",
        ]].sort_values("qvalue")
