"""Label-independent PAS pre-filter for ``switch diff`` (issue #94).

``--marker-top-n`` is a speed knob that ranks the PAS to test with the SAME
``--cluster-key`` labels the differential test then contrasts.  That is a label
double-dip: it makes every strategy anti-conservative (13.0-24.7% of null
p<0.05 under a label-permutation null vs 3.0% at ``--marker-top-n 0``), which
is why its default was flipped to 0 -- and why the speed knob went away with
it.

This module is the replacement knob.  It cuts the tested PAS set using only
properties of the count matrix POOLED OVER ALL CELLS, so it cannot leak group
information into the test:

* :func:`select_expressed_pas` takes a count matrix and an integer.  It has no
  parameter for ``cluster_key``, for a label vector, or for the AnnData object
  that carries them, so there is nothing for it to double-dip on -- the
  independence is structural, not a promise.  ``tests/`` pins that signature.

Why "expressing cells" rather than total reads or variance: the number of cells
in which a PAS is detected is the quantity the per-cell (``--count-mode
cells``) contingency table is actually built from, so the threshold removes
exactly the PAS that could never populate a usable 2xN table in either group,
and nothing else.  A read-total threshold would instead be dominated by
sequencing depth, and a variance threshold would (with two groups) start to
correlate with between-group differences -- i.e. it would drift back towards
selecting on the very signal being tested, even without seeing the labels.
"""
from __future__ import annotations

import pandas as pd

__all__ = ["select_expressed_pas"]


def select_expressed_pas(count_matrix: pd.DataFrame, min_cells: int) -> list:
    """Return the PAS detected in at least *min_cells* cells, labels ignored.

    Args:
        count_matrix: cells x PAS count DataFrame (the ``diff_df`` layout of
            :func:`ema.switch_test.runner.build_count_dfs`).  Every cell in the
            matrix is pooled; no row is ever grouped.
        min_cells: minimum number of cells with a non-zero count for a PAS to
            be kept.  ``<= 0`` keeps everything (the filter is OFF).

    Returns:
        The kept PAS ids, in the matrix's own column order.

    Note:
        There is deliberately no ``cluster_key`` / ``cluster_labels``
        parameter: the criterion is a pooled property of the matrix, so the
        selection is identical under any permutation of the group labels and
        cannot inflate the null (issue #94).
    """
    if min_cells <= 0:
        return list(count_matrix.columns)
    n_expressing = (count_matrix > 0).sum(axis=0)
    return list(count_matrix.columns[n_expressing >= min_cells])
