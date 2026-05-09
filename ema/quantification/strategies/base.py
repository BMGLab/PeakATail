"""Abstract base class for PDUI quantification strategies.

All PDUI strategies implement a single ``compute`` method that accepts a
count matrix (PAS x cells) and a PAS-to-isoform map and returns a
long-format DataFrame.

Aggregation modes (``aggregation`` parameter):
    - ``per_isoform``: keep transcript-level breakdown in the output.
    - ``per_gene``: collapse all isoforms; PAS are pooled per gene.

Isoform-collapse modes (``isoform_collapse`` parameter, used only when
``aggregation="per_gene"`` but isoform information was computed):
    - ``none``: no collapse, return one row per (gene, transcript).
    - ``mean``: mean across isoforms.
    - ``majority``: isoform with the highest expression wins.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Literal

import pandas as pd


AggregationMode = Literal["per_isoform", "per_gene"]
IsoformCollapseMode = Literal["none", "mean", "majority"]


class PDUIStrategy(ABC):
    """Base class for PDUI quantification strategies.

    Subclasses must define the class-level ``name`` attribute (used as the
    registry key) and implement :meth:`compute`.
    """

    name: str  # registry key — defined in each subclass

    @abstractmethod
    def compute(
        self,
        count_matrix: pd.DataFrame,
        pas_isoform_map: dict[int, list[tuple[str, str, int, int, int]]],
        aggregation: AggregationMode = "per_isoform",
        isoform_collapse: IsoformCollapseMode = "none",
        pseudocount: float = 0.0,
    ) -> pd.DataFrame:
        """Compute PDUI values and return a long-format DataFrame.

        Args:
            count_matrix: DataFrame of shape ``(n_pas, n_cells)`` where the
                index contains integer PAS identifiers and columns are cell
                barcodes.
            pas_isoform_map: Output of
                :func:`~ema.quantification.pas_to_isoform.map_pas_to_isoforms`::

                    {
                        pas_id: [
                            (gene_id, transcript_id, transcript_pos,
                             rank, total_pas_in_transcript),
                            ...
                        ]
                    }

            aggregation: Whether to operate at the per-isoform or per-gene
                level.
            isoform_collapse: How to collapse multiple isoforms when
                ``aggregation="per_gene"`` (only applicable when the strategy
                tracks isoform information internally).
            pseudocount: Small non-negative value added to denominators to
                avoid division by zero.  Exact semantics are strategy-specific.

        Returns:
            Long-format ``pd.DataFrame``.  Exact columns depend on the
            strategy — see subclass docstrings for the full schema.
        """
        ...
