"""Tests for issue #74: expose the D4 per-cell ("cells") count mode on the
``ema switch diff`` CLI, and show it de-pseudoreplicates the fisher test.

Two things are covered:

1. **CLI reachability** — ``--count-mode`` is a real option on ``ema switch
   diff``, accepts ``cells``/``reads``, and rejects anything else.  Before
   this fix the de-pseudoreplicated mode existed on ``FisherStrategy`` but was
   unreachable from the CLI (verified in the issue via ``--help``).

2. **Permutation-null sanity check** — under random (null) group labels, the
   legacy ``reads`` mode is anti-conservative (far too many p<0.05), while the
   ``cells`` mode is roughly calibrated.  This is a compact, deterministic
   reproduction of the miscalibration reported in the issue.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from click.testing import CliRunner

from ema.cli import main
from ema.switch_test.strategies.fisher import FisherStrategy


# --------------------------------------------------------------------------- #
# 1. CLI reachability
# --------------------------------------------------------------------------- #
def test_count_mode_flag_is_on_switch_diff_help():
    """``--count-mode`` (with a 'cells' choice) shows up on the CLI help."""
    result = CliRunner().invoke(main, ["switch", "diff", "--help"])
    assert result.exit_code == 0, result.output
    assert "--count-mode" in result.output
    assert "cells" in result.output


def test_count_mode_rejects_invalid_value(tmp_path):
    """An unknown --count-mode value is rejected by click's Choice."""
    # --h5ad has exists=True, so give it a real (empty) file; click then
    # reaches the --count-mode Choice validation and rejects the bad value
    # before the command body ever tries to read the h5ad.
    h5ad = tmp_path / "dummy.h5ad"
    h5ad.write_bytes(b"")
    result = CliRunner().invoke(
        main,
        ["switch", "diff", "-i", str(h5ad), "--count-mode", "bogus"],
    )
    assert result.exit_code != 0
    assert "bogus" in result.output
    assert "count-mode" in result.output or "count_mode" in result.output


def test_count_mode_default_is_cells():
    """Issue #74: the default was flipped reads->cells so the out-of-the-box
    fisher path is FDR-calibrated. A user can no longer UNKNOWINGLY emit
    miscalibrated q-values; the legacy 'reads' path is an explicit opt-in.

    (This intentionally replaces the earlier test_count_mode_default_is_reads,
    which pinned the pre-fix default.)"""
    from ema.cli.defaults import DEFAULTS
    assert DEFAULTS["count-mode"] == "cells"


def test_run_diff_accepts_count_mode_kwarg():
    """The library entry point threads count_mode (plumbing regression guard)."""
    import inspect
    from ema.switch_test.runner import run_diff
    assert "count_mode" in inspect.signature(run_diff).parameters


# --------------------------------------------------------------------------- #
# 2. Permutation-null sanity check
# --------------------------------------------------------------------------- #
def _make_null_dataset(n_cells: int, n_genes: int, depth: int, seed: int):
    """Build a (cells x PAS) count frame with NO true group signal.

    Each gene has 2 PAS.  Every cell deposits ``depth`` reads at exactly one
    of its gene's two PAS, chosen 50/50 at random.  There is no dependence on
    any group label, so the null is true.  Because each cell dumps ``depth``
    correlated reads, a read-level test sees ``depth``x inflated evidence
    (pseudoreplication) while a per-cell test sees one observation per cell.
    """
    rng = np.random.default_rng(seed)
    pas_cols: list[str] = []
    pas_gene_map: dict[str, str] = {}
    for g in range(n_genes):
        for k in (0, 1):
            pid = f"g{g}_p{k}"
            pas_cols.append(pid)
            pas_gene_map[pid] = f"g{g}"

    mat = np.zeros((n_cells, len(pas_cols)), dtype=int)
    for g in range(n_genes):
        choice = rng.integers(0, 2, size=n_cells)  # which PAS each cell uses
        for k in (0, 1):
            col = 2 * g + k
            mat[choice == k, col] = depth
    cells = [f"cell{i}" for i in range(n_cells)]
    cm = pd.DataFrame(mat, index=cells, columns=pas_cols)
    return cm, pas_gene_map


def _null_reject_rate(count_mode: str, n_perm: int = 12) -> float:
    """Fraction of PAS tests with p<0.05 across permutation-null replicates."""
    strat = FisherStrategy()
    n_cells = 120
    total = 0
    rejected = 0
    for perm in range(n_perm):
        cm, pas_gene_map = _make_null_dataset(
            n_cells=n_cells, n_genes=15, depth=80, seed=1000 + perm
        )
        rng = np.random.default_rng(9000 + perm)
        labels = pd.Series(
            np.where(rng.permutation(n_cells) < n_cells // 2, "A", "B"),
            index=cm.index,
        )
        res = strat.test(
            count_matrix=cm,
            cluster_labels=labels,
            cluster1="A",
            cluster2="B",
            min_cells_per_group=5,
            pas_gene_map=pas_gene_map,
            count_mode=count_mode,
        )
        if res.empty:
            continue
        total += len(res)
        rejected += int((res["pvalue"] < 0.05).sum())
    assert total > 0
    return rejected / total


def test_reads_mode_is_anticonservative_under_null():
    """Legacy reads mode over-rejects under the null (reproduces issue #74)."""
    rate = _null_reject_rate("reads")
    # Nominal would be ~0.05; pseudoreplication pushes it far higher.
    assert rate > 0.20, f"expected reads mode to over-reject, got {rate:.3f}"


def test_cells_mode_is_calibrated_under_null():
    """The per-cell (D4) mode is roughly calibrated under the null."""
    rate = _null_reject_rate("cells")
    assert rate < 0.15, f"cells mode should be ~calibrated, got {rate:.3f}"


def test_cells_mode_less_anticonservative_than_reads():
    """Direct comparison: cells mode rejects the null far less than reads."""
    reads_rate = _null_reject_rate("reads")
    cells_rate = _null_reject_rate("cells")
    assert cells_rate < reads_rate, (
        f"cells ({cells_rate:.3f}) should be < reads ({reads_rate:.3f})"
    )


# --------------------------------------------------------------------------- #
# 3. FDR-calibration under permutation null (issue #74 headline metric)
#
# The issue is stated in terms of q-values, not raw p-values: under label
# permutation, read-level fisher reports q<0.05 hits in *100% of null runs*
# (mean 42 hits/run). The calibrated fix must drive that per-run hit rate to
# ~0. These tests assert exactly that, so they are the evidence the fix works
# (not merely that a flag exists). Small n + fixed seeds keep them fast and
# deterministic.
# --------------------------------------------------------------------------- #
def _null_qhit_stats(count_mode: str, n_perm: int = 12) -> tuple[float, float]:
    """Return (fraction_of_null_runs_with>=1 q<0.05 hit, mean q<0.05 hits/run).

    Each run rebuilds a null dataset AND draws a fresh random group label,
    then BH-corrects within the run -- mirroring the issue's protocol where a
    full permutation is scored end-to-end.
    """
    strat = FisherStrategy()
    n_cells = 120
    runs_with_hit = 0
    hits_per_run: list[int] = []
    for perm in range(n_perm):
        cm, pas_gene_map = _make_null_dataset(
            n_cells=n_cells, n_genes=15, depth=80, seed=1000 + perm
        )
        rng = np.random.default_rng(9000 + perm)
        labels = pd.Series(
            np.where(rng.permutation(n_cells) < n_cells // 2, "A", "B"),
            index=cm.index,
        )
        res = strat.test(
            count_matrix=cm,
            cluster_labels=labels,
            cluster1="A",
            cluster2="B",
            min_cells_per_group=5,
            pas_gene_map=pas_gene_map,
            count_mode=count_mode,
        )
        if res.empty:
            continue
        n_hits = int((res["qvalue"] < 0.05).sum())
        hits_per_run.append(n_hits)
        if n_hits > 0:
            runs_with_hit += 1
    assert hits_per_run, "no non-empty null runs produced"
    return runs_with_hit / len(hits_per_run), float(np.mean(hits_per_run))


def test_reads_mode_reports_qhits_in_most_null_runs():
    """Reproduce issue #74: reads mode reports q<0.05 hits in ~all null runs."""
    frac_runs, mean_hits = _null_qhit_stats("reads")
    assert frac_runs > 0.75, (
        f"expected reads mode to report q<0.05 hits in most null runs, "
        f"got {frac_runs:.2f} of runs (mean {mean_hits:.1f} hits/run)"
    )
    assert mean_hits > 1.0, f"expected many false q-hits/run, got {mean_hits:.1f}"


def test_cells_mode_qhits_near_zero_under_null():
    """The fix: cells mode reports q<0.05 hits in ~no null runs (calibrated)."""
    frac_runs, mean_hits = _null_qhit_stats("cells")
    # Calibrated: a genuine BH-controlled null should almost never emit a hit,
    # and certainly not in 100% of runs the way reads mode does.
    assert frac_runs <= 0.15, (
        f"cells mode should report q<0.05 hits in ~no null runs, got "
        f"{frac_runs:.2f} of runs (mean {mean_hits:.2f} hits/run)"
    )
    assert mean_hits < 1.0, f"cells mode mean q-hits/run too high: {mean_hits:.2f}"


def test_cells_mode_far_fewer_qhits_than_reads_under_null():
    """Head-to-head: the calibrated default rejects the null far less."""
    reads_frac, reads_mean = _null_qhit_stats("reads")
    cells_frac, cells_mean = _null_qhit_stats("cells")
    assert cells_frac < reads_frac, (
        f"cells run-hit rate ({cells_frac:.2f}) must be < reads "
        f"({reads_frac:.2f})"
    )
    assert cells_mean < reads_mean, (
        f"cells mean q-hits/run ({cells_mean:.2f}) must be < reads "
        f"({reads_mean:.2f})"
    )


def test_invalid_count_mode_raises_in_strategy():
    cm, pas_gene_map = _make_null_dataset(60, 5, 40, seed=0)
    labels = pd.Series(["A"] * 30 + ["B"] * 30, index=cm.index)
    with pytest.raises(ValueError, match="count_mode"):
        FisherStrategy().test(
            count_matrix=cm,
            cluster_labels=labels,
            cluster1="A",
            cluster2="B",
            pas_gene_map=pas_gene_map,
            count_mode="nonsense",
        )
