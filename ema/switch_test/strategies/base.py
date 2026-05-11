"""Abstract base class for differential APA testing strategies.

All strategies in this package implement this interface so the orchestration
layer can select and invoke them by name without branching on method.
"""

from __future__ import annotations

from abc import ABC, abstractmethod

import pandas as pd


class DiffAPAStrategy(ABC):
    """Abstract strategy for differential alternative polyadenylation testing.

    Subclasses register themselves in the package registry via
    :func:`ema.switch_test.strategies.register_diff_strategy`.

    Class attributes:
        name: Registry key used by :func:`get_diff_strategy`.
        supports_multi_condition: True when the strategy handles >2 groups
            natively (e.g. omnibus LRT).  False for pairwise methods.
    """

    name: str
    supports_multi_condition: bool

    @abstractmethod
    def test(
        self,
        count_matrix: pd.DataFrame,
        cluster_labels: pd.Series,
        cluster1: str | None = None,
        cluster2: str | None = None,
        min_cells_per_group: int = 10,
        pas_gene_map: dict[str, str] | None = None,
        n_jobs: int = -1,
    ) -> pd.DataFrame:
        """Run differential APA test across PAS rows.

        Parameters
        ----------
        count_matrix:
            Shape ``(n_cells, n_pas)``.  Columns are PAS identifiers.
            Values are raw UMI / read counts (integers).
        cluster_labels:
            Series indexed by cell barcode (same index as ``count_matrix``
            rows), values are cluster identifiers (any hashable type).
        cluster1:
            Label of the first group.  Required for pairwise strategies
            (``supports_multi_condition=False``); ignored by multi-condition
            strategies.
        cluster2:
            Label of the second group.  Required for pairwise strategies;
            ignored by multi-condition strategies.
        min_cells_per_group:
            Minimum number of cells with non-zero counts required in each
            group for a PAS to be tested.  PAS that do not meet this
            threshold are dropped silently.
        n_jobs:
            Number of parallel workers for per-PAS fits.
            ``-1`` means use all available CPUs (joblib convention).

        Returns
        -------
        pd.DataFrame
            Indexed by ``pas_id`` (column names of ``count_matrix``).
            Required columns: ``pvalue``, ``qvalue``, ``n_cells``.
            Pairwise strategies also include: ``log2fc``, ``odds_ratio``
            (fisher), ``dispersion`` (NB), ``test_stat``.
            Multi-condition strategies include: ``test_stat``, ``df``,
            ``dispersion``.

        Raises
        ------
        ValueError
            If pairwise strategy called without ``cluster1``/``cluster2``.
        """
        ...
