"""B2 regression: annotate() resolves output paths at call time.

Bug B2: annotate() bound ``directory_config.annotated_matrix`` / ``.pasbed`` as
DEFAULT ARGUMENT values, captured at import time before ``set_directory_config``
runs — so a run wrote stray ``emaout/`` shadow files instead of into the
timestamped run dir. Sentinel defaults + call-time resolution fix this.
"""
from __future__ import annotations

import inspect
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ema.annotate.annotate import annotate


def test_annotate_defaults_are_sentinels_not_frozen_paths() -> None:
    """The signature must not freeze directory_config paths at import time."""
    sig = inspect.signature(annotate)
    assert sig.parameters["annotated_matrix"].default is None
    assert sig.parameters["pasbed_dir"].default is None


def test_annotate_writes_into_current_run_dir(tmp_path: Path) -> None:
    import ema.config as cfg

    run = tmp_path / "run_TS"
    (run / "05_annotated_matrix").mkdir(parents=True)
    cfg.set_directory_config(output_dir=run)

    # Minimal inputs: 2 PAS x 2 cells; both PAS have a gene assignment.
    pas_ids = np.array([1, 2])
    matrix = sp.csr_matrix(np.array([[3, 0], [0, 5]], dtype=np.int32))
    collist = ["CB1", "CB2"]
    genes = pd.DataFrame({"gene_id": ["ENSG1", "ENSG2"]}, index=[1, 2])

    # A pasbed on disk at the resolved (run-scoped) path.
    pasbed = Path(cfg.directory_config.pasbed)
    pasbed.parent.mkdir(parents=True, exist_ok=True)
    pasbed.write_text("chr1\t10\t11\t1\t0\t+\nchr1\t20\t21\t2\t0\t-\n")

    annotate(sparse_matrix=matrix, pas_ids=pas_ids, collist=collist, genes=genes)

    # Outputs must land under the run dir, NOT under a stray ./emaout/.
    assert Path(cfg.directory_config.annotated_matrix).exists()
    assert str(run) in str(cfg.directory_config.annotated_matrix)
    assert not Path("emaout/annotated_matrix.mtx").exists()
