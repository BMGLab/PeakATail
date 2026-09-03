"""Issue #94 regression: `--marker-top-n` pre-selection is a label double-dip.

``ema switch diff`` used to default to ``--marker-top-n 200``: before testing,
the PAS matrix was restricted to the union of the top-200 marker PAS per
cluster, ranked by ``scanpy.tl.rank_genes_groups`` on the SAME
``--cluster-key`` labels the differential test then contrasts.  The PAS that
entered the test were therefore exactly the ones that already looked
cluster-associated by chance, and restricting the matrix additionally shrank
the within-gene Fisher denominator ("the rest of the gene" became the same
label-selected subset).  Every strategy came out anti-conservative.

These tests run a **label-permutation null** -- the count matrix is fixed and
carries no group structure whatsoever, only the cluster labels are shuffled,
so every rejection is a false positive by construction -- straight through the
real :func:`ema.switch_test.runner.run_diff` with the real ``fisher``
strategy.  They assert that

* at the DEFAULT settings (``DEFAULTS["marker-top-n"]``) the null p < 0.05
  rate stays near nominal and no permutation produces a q < 0.05 hit, and
* ``--marker-top-n 200`` inflates that rate several-fold.

Everything is deterministic: the counts come from a fixed RandomState, each
permutation from its own fixed seed, and no file is read from disk.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

import ema.switch_test.runner as runner_mod
from ema.cli.defaults import DEFAULTS

# ---------------------------------------------------------------------------
# Synthetic null fixture
# ---------------------------------------------------------------------------
# Sized so that top-200-per-cluster is genuinely selective: with 2 clusters the
# two per-cluster top-N lists are mirror images of one another, so the union is
# ~2 * top_n PAS.  400 of 2000 PAS survive -- a 5x cut.  (With <= 400 PAS the
# union would be the whole matrix and the double-dip would be invisible.)
N_GENES = 500
PAS_PER_GENE = 4
N_PAS = N_GENES * PAS_PER_GENE
N_CELLS = 120
COUNT_SEED = 12345
POISSON_MEAN = 1.0

#: Permutation seeds. Fixed, so the whole test is reproducible.
PERM_SEEDS = (0, 1, 2, 3, 4)

#: The value that used to be the default and is now an explicit opt-in.
LEGACY_TOP_N = 200

#: Nominal alpha of the null.  The Fisher exact test on small per-cell tables
#: is discrete and conservative, so a calibrated run lands at or below this.
ALPHA = 0.05

#: Slack over ALPHA for the calibrated default.  Measured ~3.0-4.2% per seed;
#: the observed rate at top-n 200 is 13.7-24.5%, so the gap is wide.
CALIBRATED_MAX = 0.07


def _null_adata(perm_seed: int) -> ad.AnnData:
    """Cells x PAS counts with NO group structure, plus permuted labels.

    The matrix is identical for every ``perm_seed`` -- only the label vector is
    reshuffled -- so any signal the test finds is purely a label artefact.
    """
    rng = np.random.RandomState(COUNT_SEED)
    X = rng.poisson(POISSON_MEAN, size=(N_CELLS, N_PAS)).astype(float)

    pas_ids = [str(i + 1) for i in range(N_PAS)]
    gene_ids = [f"GENE_{i // PAS_PER_GENE}" for i in range(N_PAS)]
    var = pd.DataFrame({"gene_id": gene_ids}, index=pas_ids)

    labels = np.array(["A"] * (N_CELLS // 2) + ["B"] * (N_CELLS - N_CELLS // 2))
    np.random.RandomState(perm_seed).shuffle(labels)
    obs = pd.DataFrame(
        {"leiden": labels}, index=[f"cell{i}" for i in range(N_CELLS)]
    )
    return ad.AnnData(X=X, obs=obs, var=var)


def _run_diff(adata: ad.AnnData, marker_top_n: int, out_dir: Path) -> dict:
    """One replicate through the real ``run_diff`` (no h5ad ever opened)."""
    real_read_h5ad = runner_mod.ad.read_h5ad
    runner_mod.ad.read_h5ad = lambda _p: adata
    try:
        return runner_mod.run_diff(
            h5ad_paths=["permuted.h5ad"],  # never opened -- read_h5ad is stubbed
            pasbed=None,
            gtf=None,
            output_dir=str(out_dir),
            cluster_pairs=None,
            cluster_key="leiden",
            marker_top_n=marker_top_n,
            marker_method="wilcoxon",
            strategy="fisher",
            count_mode="cells",
            fdr=ALPHA,
            threads=1,
            per_worker_mb=300,
            min_cells_per_group=5,
        )
    finally:
        runner_mod.ad.read_h5ad = real_read_h5ad


def _null_stats(marker_top_n: int, tmp_root: Path) -> list[tuple[float, int]]:
    """``[(fraction of null p<0.05, number of q<0.05 hits), ...]`` per seed."""
    out = []
    for seed in PERM_SEEDS:
        pair_results = _run_diff(
            _null_adata(seed),
            marker_top_n,
            tmp_root / f"perm{seed}_top{marker_top_n}",
        )
        assert list(pair_results) == [("A", "B")]
        df = pair_results[("A", "B")]
        assert not df.empty, "permutation replicate produced no tested PAS"
        out.append(
            (float((df["pvalue"] < ALPHA).mean()), int((df["qvalue"] < ALPHA).sum()))
        )
    return out


@pytest.fixture(scope="module")
def null_stats(tmp_path_factory) -> dict[int, list[tuple[float, int]]]:
    """Run the permutation null once per configuration, share across tests.

    Always covers ``0`` and the legacy ``200`` (so the inflation test states a
    permanent property of the flag) plus whatever the CLI default currently
    is (so the calibration test tracks the default rather than a constant).
    """
    root = tmp_path_factory.mktemp("i94_null")
    wanted = sorted({0, LEGACY_TOP_N, int(DEFAULTS["marker-top-n"])})
    return {top_n: _null_stats(top_n, root) for top_n in wanted}


# ---------------------------------------------------------------------------
# The default itself
# ---------------------------------------------------------------------------


def test_marker_top_n_default_is_zero():
    """A flagless `ema switch diff` must do no label-based pre-selection."""
    from ema.cli.config_schema import RunConfig

    assert DEFAULTS["marker-top-n"] == 0, (
        "--marker-top-n must default to 0: any non-zero default silently "
        "double-dips on the cluster labels (issue #94)"
    )
    assert RunConfig().marker_top_n == 0


# ---------------------------------------------------------------------------
# Permutation null
# ---------------------------------------------------------------------------


def test_permutation_null_is_calibrated_at_defaults(null_stats):
    """Under the null, the DEFAULT settings must not over-reject."""
    stats = null_stats[int(DEFAULTS["marker-top-n"])]
    rates = [r for r, _ in stats]
    hits = [h for _, h in stats]

    assert max(rates) <= CALIBRATED_MAX, (
        f"null p<{ALPHA} rate at the default settings is {rates} -- the "
        f"default is not FDR-controlled (issue #94)"
    )
    assert sum(hits) == 0, (
        f"the label-permutation null produced q<{ALPHA} 'hits' {hits} at the "
        f"default settings -- there is nothing true to find"
    )


def test_marker_top_n_200_inflates_the_null(null_stats):
    """The old default is anti-conservative; that is why it is opt-in now."""
    off_rates = [r for r, _ in null_stats[0]]
    legacy_rates = [r for r, _ in null_stats[LEGACY_TOP_N]]

    assert min(legacy_rates) > max(off_rates), (
        f"--marker-top-n {LEGACY_TOP_N} was expected to inflate the null "
        f"({legacy_rates}) above --marker-top-n 0 ({off_rates})"
    )
    assert np.mean(legacy_rates) > 2.0 * np.mean(off_rates), (
        f"--marker-top-n {LEGACY_TOP_N} null p<{ALPHA} rate "
        f"{np.mean(legacy_rates):.3f} vs {np.mean(off_rates):.3f} at 0 "
        f"-- expected a several-fold inflation"
    )
    # ... and it is well above nominal in absolute terms, too.
    assert np.mean(legacy_rates) > 2.0 * ALPHA


# ---------------------------------------------------------------------------
# Loud warning on the opt-in path
# ---------------------------------------------------------------------------


def _tiny_adata() -> ad.AnnData:
    pas_ids = ["1", "2", "3", "4"]
    var = pd.DataFrame({"gene_id": ["G1", "G1", "G2", "G2"]}, index=pas_ids)
    labels = ["A"] * 6 + ["B"] * 6
    obs = pd.DataFrame({"leiden": labels}, index=[f"c{i}" for i in range(12)])
    X = np.random.RandomState(7).poisson(3.0, size=(12, 4)).astype(float)
    return ad.AnnData(X=X, obs=obs, var=var)


def test_nonzero_marker_top_n_warns_loudly(tmp_path, caplog):
    with caplog.at_level(logging.WARNING, logger="ema.switch_test.runner"):
        _run_diff(_tiny_adata(), 2, tmp_path / "warn")
    text = "\n".join(
        r.getMessage() for r in caplog.records if r.levelno >= logging.WARNING
    )
    assert "double-dip" in text
    assert "#94" in text
    assert "marker_top_n=0" in text


def test_default_marker_top_n_does_not_warn(tmp_path, caplog):
    with caplog.at_level(logging.WARNING, logger="ema.switch_test.runner"):
        _run_diff(_tiny_adata(), int(DEFAULTS["marker-top-n"]), tmp_path / "nowarn")
    text = "\n".join(
        r.getMessage() for r in caplog.records if r.levelno >= logging.WARNING
    )
    assert "double-dip" not in text
