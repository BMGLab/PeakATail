"""ema.data -- the E4 typed-facade package.

Exposes :class:`Run` (see ``ema/data/run.py``), the manifest-driven,
pandas-based read facade over a PeakATail run directory. Supersedes the
path-hardcoding in ``scripts/harvest_report_data.py`` and
``scripts/harvest_matrix_per_dataset.py``. Engine-self-contained: nothing in
this package imports ``peakatail_contract`` or any ``peakatail-hub`` package.
"""
from __future__ import annotations

from ema.data.run import Run, RunReadError

__all__ = ["Run", "RunReadError"]
