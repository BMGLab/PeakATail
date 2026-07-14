"""Cell-typing utilities (A1).

The reusable, dependency-light core of gene-expression (GEX) marker cell typing
and its concordance against PAS-based clusters. The heavy GEX I/O (STARsolo
Solo.out pooling) and scanpy clustering / ``score_genes`` steps that feed this
core are deliberately kept out of the package until validated on real data —
see :mod:`ema.celltype.scoring` for the parts that are unit-tested now.
"""
from ema.celltype.scoring import (
    assign_celltypes,
    concordance,
    zscore_across_cells,
)

__all__ = ["assign_celltypes", "concordance", "zscore_across_cells"]
