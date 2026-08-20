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
import logging
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, fisher_exact

from ema.switch_test.strategies import register_diff_strategy
from ema.switch_test.strategies.base import DiffAPAStrategy

log = logging.getLogger(__name__)


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
        pas_gene_map: dict[str, str] | None = None,
        n_jobs: int = -1,  # Fisher is fast; n_jobs is accepted but unused.
        count_mode: str = "cells",
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

        # Within-gene Fisher exact test on read counts.  For each gene G and
        # each PAS p of G, the 2x2 table is
        #
        #                          cluster1            cluster2
        #   reads at p           reads_p_c1          reads_p_c2
        #   reads at OTHER PAS   reads_otherG_c1     reads_otherG_c2
        #                        (same gene)         (same gene)
        #
        # This is the canonical APA differential-usage framing (matches the
        # original PeakATail design on develop branch + DEXSeq for splicing):
        # it asks "for this gene, does this PAS get used differentially?".
        # Without pas_gene_map we cannot group within-gene, so fall back to
        # the global "this PAS vs all OTHER PAS" framing with a warning.
        cm1 = count_matrix.loc[cells1]
        cm2 = count_matrix.loc[cells2]
        n1 = int(len(cells1))
        n2 = int(len(cells2))

        agg1 = cm1.sum(axis=0).astype(int)            # reads per PAS in c1
        agg2 = cm2.sum(axis=0).astype(int)            # reads per PAS in c2
        expr1 = (cm1 > 0).sum(axis=0).astype(int)     # cells expressing PAS in c1
        expr2 = (cm2 > 0).sum(axis=0).astype(int)     # cells expressing PAS in c2

        # Group PAS by gene.  When pas_gene_map is missing we treat every
        # PAS as its own (singleton) gene which collapses the within-gene
        # comparison to the legacy global one — still correct but loses
        # APA-specific framing.
        pas_ids = list(count_matrix.columns)
        if pas_gene_map is None:
            log.warning(
                "FisherStrategy: no pas_gene_map supplied; falling back to "
                "global (cross-gene) comparison.  Pass pas_gene_map for the "
                "intended within-gene APA test."
            )
            pas_by_gene: dict[str, list] = {"__global__": list(pas_ids)}
        else:
            pas_by_gene = {}
            for p in pas_ids:
                g = pas_gene_map.get(str(p))
                if not g:
                    continue
                pas_by_gene.setdefault(str(g), []).append(p)

        # D4: read counts within a cell are correlated, so a read-based 2x2
        # table pseudoreplicates (significance scales with sequencing depth:
        # Spearman(cells, %sig)=+0.68 on real data). ``count_mode="cells"`` (the
        # DEFAULT since issue #74) is the de-pseudoreplicated mode — it builds
        # the table from per-cell detection among gene-expressing cells (each
        # cell counted once), which is FDR-calibrated under a permutation null.
        # ``count_mode="reads"`` (legacy, opt-in) is anti-conservative and is
        # retained only for backward comparison / as a ranking screen; the
        # inflation regression tests pass it explicitly.
        if count_mode not in ("cells", "reads"):
            raise ValueError(f"count_mode must be 'cells' or 'reads', got {count_mode!r}")

        results: list[dict] = []
        for gene_id, pas_in_gene in pas_by_gene.items():
            if len(pas_in_gene) < 2:
                continue  # need >=2 PAS in the gene to compare within-gene

            # Reads columns are always reported (for reference), but the TABLE +
            # proportions use the chosen count_mode.
            gene_reads_c1 = int(sum(int(agg1[p]) for p in pas_in_gene))
            gene_reads_c2 = int(sum(int(agg2[p]) for p in pas_in_gene))
            if count_mode == "cells":
                gene_total_c1 = int((cm1[pas_in_gene] > 0).any(axis=1).sum())
                gene_total_c2 = int((cm2[pas_in_gene] > 0).any(axis=1).sum())
            else:
                gene_total_c1, gene_total_c2 = gene_reads_c1, gene_reads_c2
            if gene_total_c1 + gene_total_c2 == 0:
                continue

            for p in pas_in_gene:
                reads_p_c1 = int(agg1[p])
                reads_p_c2 = int(agg2[p])
                if count_mode == "cells":
                    count_p_c1 = int(expr1[p])
                    count_p_c2 = int(expr2[p])
                else:
                    count_p_c1, count_p_c2 = reads_p_c1, reads_p_c2
                # Skip PAS uninformative in BOTH clusters of the gene.
                if count_p_c1 + count_p_c2 == 0:
                    continue
                other_c1 = gene_total_c1 - count_p_c1
                other_c2 = gene_total_c2 - count_p_c2
                table = np.array([[count_p_c1, count_p_c2],
                                  [other_c1,   other_c2]])
                odds_ratio, pvalue = fisher_exact(table, alternative="two-sided")

                # Within-gene proportions in the chosen unit (cells or reads).
                prop1 = count_p_c1 / gene_total_c1 if gene_total_c1 > 0 else 0.0
                prop2 = count_p_c2 / gene_total_c2 if gene_total_c2 > 0 else 0.0
                delta_prop = prop1 - prop2

                # log2fc of the within-gene proportion.  eps prevents log(0)
                # while keeping the value finite when one cluster is empty.
                eps = 1.0 / max(gene_total_c1 + 1, gene_total_c2 + 1)
                log2fc = float(np.log2((prop2 + eps) / (prop1 + eps)))

                results.append({
                    "pas_id": p,
                    "pvalue": pvalue,
                    "odds_ratio": float(odds_ratio),
                    "delta_proportion": delta_prop,
                    "log2fc": log2fc,
                    "n_cells_cluster1": n1,
                    "n_cells_cluster2": n2,
                    "n_cells_expr_cluster1": int(expr1[p]),
                    "n_cells_expr_cluster2": int(expr2[p]),
                    "n_reads_pas_cluster1": reads_p_c1,
                    "n_reads_pas_cluster2": reads_p_c2,
                    "n_reads_gene_cluster1": gene_reads_c1,
                    "n_reads_gene_cluster2": gene_reads_c2,
                    "n_cells": n1 + n2,
                })

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
            "n_reads_gene_cluster1", "n_reads_gene_cluster2",
            "odds_ratio", "delta_proportion", "log2fc",
        ]].sort_values("qvalue")
