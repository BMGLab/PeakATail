"""Issue #94: `--prefilter-min-cells`, the label-INDEPENDENT speed pre-filter.

``--marker-top-n`` used to be the way to make ``peakatail switch diff`` fast on a
big matrix: test only the top-N marker PAS per cluster.  Those markers are
ranked with the SAME ``--cluster-key`` labels the differential test then
contrasts, so the flag is a label double-dip -- see
``tests/test_marker_top_n_double_dip_i94.py``, which measures the resulting
null inflation.  PR #105 defaulted it to 0, which fixed the statistics and
removed the speed knob.

``--prefilter-min-cells N`` is the replacement knob: keep only the PAS detected
(count > 0) in at least *N* cells, counted over **all cells pooled**.  The
criterion is a property of the count matrix alone, so it cannot leak group
information into the test.  These tests pin the three things that makes true:

1. **Structural label independence.**
   :func:`ema.switch_test.prefilter.select_expressed_pas` takes a matrix and an
   integer -- there is no ``cluster_key`` / label parameter for it to peek at,
   and the module names none.  The consequence is measured too: the selected
   PAS set is identical under every permutation of the labels.
2. **It does not inflate the null.**  Under the same label-permutation null the
   marker flag is measured on (fixed counts, shuffled labels, so every
   rejection is false by construction), the pre-filter leaves the null p < 0.05
   rate at the unfiltered level while ``--marker-top-n 200`` multiplies it
   several-fold.  Stronger still: each surviving p-value is *bit-identical* to
   the one the unfiltered run reports for that PAS, and each q-value is exactly
   Benjamini-Hochberg over that same label-blind subset -- the pre-filter only
   ever removes hypotheses, it never changes a test.
3. **It really cuts the work**, and **OFF changes nothing**: at the default 0
   the result is byte-identical to a run that never heard of the flag.

Everything is deterministic: fixed count seed, fixed permutation seeds, no file
read from disk.
"""

from __future__ import annotations

import inspect
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy.stats import false_discovery_control

sys.path.insert(0, str(Path(__file__).parent.parent))

import ema.switch_test.runner as runner_mod
from ema.cli.defaults import DEFAULTS
from ema.switch_test.prefilter import select_expressed_pas

# ---------------------------------------------------------------------------
# Synthetic null fixture (same shape as tests/test_marker_top_n_double_dip_i94)
# ---------------------------------------------------------------------------
N_GENES = 500
PAS_PER_GENE = 4
N_PAS = N_GENES * PAS_PER_GENE
N_CELLS = 120
COUNT_SEED = 12345
POISSON_MEAN = 1.0

PERM_SEEDS = (0, 1, 2, 3, 4)

#: The value that used to be the default and is now an explicit opt-in.
LEGACY_TOP_N = 200

#: Pre-filter threshold used throughout.  With Poisson(1) counts over 120 cells
#: a PAS is detected in ~76 cells on average, so ">= 80 cells" is genuinely
#: selective (469 of 2000 PAS survive, a 4.3x cut) without emptying the matrix.
PREFILTER_MIN_CELLS = 80

ALPHA = 0.05

#: Slack over ALPHA for a calibrated configuration; identical to the value
#: ``test_marker_top_n_double_dip_i94.py`` uses for the calibrated default.
CALIBRATED_MAX = 0.07


def _null_adata(perm_seed: int) -> ad.AnnData:
    """Cells x PAS counts with NO group structure, plus permuted labels."""
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


def _run_diff(adata: ad.AnnData, out_dir: Path, **kwargs) -> dict:
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
            marker_top_n=kwargs.pop("marker_top_n", 0),
            marker_method="wilcoxon",
            strategy="fisher",
            count_mode="cells",
            fdr=ALPHA,
            threads=1,
            per_worker_mb=300,
            min_cells_per_group=5,
            **kwargs,
        )
    finally:
        runner_mod.ad.read_h5ad = real_read_h5ad


#: The three configurations compared below: unfiltered, label-independent
#: pre-filter, and the label double-dip that used to be the default.
_ARMS: dict[str, dict] = {
    "off": {},
    "prefilter": {"prefilter_min_cells": PREFILTER_MIN_CELLS},
    "marker200": {"marker_top_n": LEGACY_TOP_N},
}


@pytest.fixture(scope="module")
def null_runs(tmp_path_factory) -> dict[str, dict[int, pd.DataFrame]]:
    """``{arm: {perm_seed: result frame}}`` -- run once, shared by every test."""
    root = tmp_path_factory.mktemp("i94_prefilter")
    out: dict[str, dict[int, pd.DataFrame]] = {}
    for arm, kwargs in _ARMS.items():
        per_seed = {}
        for seed in PERM_SEEDS:
            res = _run_diff(
                _null_adata(seed), root / f"{arm}_perm{seed}", **dict(kwargs)
            )
            assert list(res) == [("A", "B")]
            df = res[("A", "B")]
            assert not df.empty, f"{arm}/perm{seed} produced no tested PAS"
            per_seed[seed] = df
        out[arm] = per_seed
    return out


def _p_rate(df: pd.DataFrame) -> float:
    return float((df["pvalue"] < ALPHA).mean())


def _q_hits(df: pd.DataFrame) -> int:
    return int((df["qvalue"] < ALPHA).sum())


# ---------------------------------------------------------------------------
# The flag itself
# ---------------------------------------------------------------------------


def test_prefilter_defaults_to_off():
    """A flagless `peakatail switch diff` must behave exactly as before."""
    from ema.cli.config_schema import RunConfig

    assert DEFAULTS["prefilter-min-cells"] == 0
    assert RunConfig().prefilter_min_cells == 0
    assert (
        inspect.signature(runner_mod.run_diff)
        .parameters["prefilter_min_cells"].default == 0
    )


def test_marker_top_n_help_points_at_the_prefilter():
    """`--marker-top-n`'s help must name the safe speed alternative."""
    from click.testing import CliRunner

    from ema.cli.switch_diff import diff

    text = CliRunner().invoke(diff, ["--help"]).output
    flat = " ".join(text.split())
    assert "--prefilter-min-cells" in flat
    # the pointer lives in --marker-top-n's own help paragraph
    marker_help = flat.split("--marker-top-n INTEGER", 1)[1].split(
        "--prefilter-min-cells INTEGER", 1
    )[0]
    assert "--prefilter-min-cells" in marker_help, (
        "--marker-top-n's help does not point at --prefilter-min-cells as the "
        "label-independent speed alternative (issue #94)"
    )


# ---------------------------------------------------------------------------
# 1. Structural label independence
# ---------------------------------------------------------------------------


def test_selector_has_no_way_to_see_the_labels():
    """The selector's signature makes the double-dip impossible, not unlikely."""
    params = list(inspect.signature(select_expressed_pas).parameters)
    assert params == ["count_matrix", "min_cells"], (
        f"select_expressed_pas{tuple(params)} -- the label-independent "
        "pre-filter must take only the count matrix and a threshold; any "
        "labels/adata/cluster_key parameter reopens issue #94"
    )
    src = Path(runner_mod.__file__).with_name("prefilter.py").read_text()
    body = "\n".join(
        line for line in src.splitlines()
        if not line.lstrip().startswith("#")
    )
    # The docstrings mention cluster_key to say it is absent; the *code* must
    # not name any label-carrying object at all.
    code = body.split('"""')
    executable = "".join(code[::2])
    for forbidden in ("cluster_key", "cluster_labels", "obs", "adata", "groups"):
        assert forbidden not in executable, (
            f"ema/switch_test/prefilter.py mentions {forbidden!r} in code -- "
            "the pre-filter must be computed from the pooled matrix only"
        )


def test_selection_is_invariant_under_label_permutation(null_runs):
    """Measured consequence: shuffling the labels changes nothing selected."""
    tested = {
        seed: set(df["pas_id"]) for seed, df in null_runs["prefilter"].items()
    }
    sizes = {seed: len(s) for seed, s in tested.items()}
    first = tested[PERM_SEEDS[0]]
    assert all(s == first for s in tested.values()), (
        f"the pre-filtered PAS set differs between label permutations "
        f"({sizes}) -- the selection is not label-independent"
    )
    # ... unlike the marker path, which selects a different set every time.
    marker_sets = [set(df["pas_id"]) for df in null_runs["marker200"].values()]
    assert any(s != marker_sets[0] for s in marker_sets), (
        "sanity check: --marker-top-n is supposed to track the labels"
    )


# ---------------------------------------------------------------------------
# 2. The pre-filter does not inflate the null
# ---------------------------------------------------------------------------


def test_permutation_null_is_not_inflated_by_the_prefilter(null_runs):
    """The whole point: speed without the anti-conservatism of #94."""
    off = [_p_rate(df) for df in null_runs["off"].values()]
    pre = [_p_rate(df) for df in null_runs["prefilter"].values()]
    marker = [_p_rate(df) for df in null_runs["marker200"].values()]

    assert max(pre) <= CALIBRATED_MAX, (
        f"null p<{ALPHA} rate with --prefilter-min-cells "
        f"{PREFILTER_MIN_CELLS} is {pre} -- the pre-filter is not calibrated"
    )
    assert np.mean(pre) < 1.5 * np.mean(off), (
        f"the pre-filter moved the null p<{ALPHA} rate from "
        f"{np.mean(off):.4f} (unfiltered) to {np.mean(pre):.4f}"
    )
    # ... whereas the label double-dip is several-fold above both.
    assert np.mean(marker) > 2.0 * np.mean(pre), (
        f"sanity check: --marker-top-n {LEGACY_TOP_N} ({np.mean(marker):.4f}) "
        f"should dwarf the label-independent pre-filter ({np.mean(pre):.4f})"
    )

    # q-values: the double-dip fires under the null; the pre-filter must not
    # fire more often than the unfiltered run does on the same hypotheses.
    q_pre = sum(_q_hits(df) for df in null_runs["prefilter"].values())
    q_bh_of_subset = 0
    for seed, df in null_runs["prefilter"].items():
        full = null_runs["off"][seed].set_index("pas_id")
        p_sub = full.loc[df["pas_id"], "pvalue"].values
        q_bh_of_subset += int(
            (false_discovery_control(p_sub, method="bh") < ALPHA).sum()
        )
    assert q_pre == q_bh_of_subset, (
        f"the pre-filtered run reports {q_pre} null q<{ALPHA} 'hits' but "
        f"Benjamini-Hochberg over the unfiltered p-values of exactly those "
        f"PAS gives {q_bh_of_subset} -- the pre-filter is doing more than "
        "removing hypotheses"
    )


def test_prefilter_only_removes_hypotheses_never_changes_one(null_runs):
    """Every surviving row must equal the unfiltered run's row for that PAS."""
    for seed, df in null_runs["prefilter"].items():
        full = null_runs["off"][seed].set_index("pas_id")
        kept = df.set_index("pas_id")
        assert set(kept.index) <= set(full.index)
        sub = full.loc[kept.index]
        for col in ("pvalue", "n_reads_pas_cluster1", "n_reads_gene_cluster1",
                    "n_reads_gene_cluster2", "odds_ratio"):
            assert np.allclose(
                sub[col].values.astype(float),
                kept[col].values.astype(float),
                equal_nan=True,
            ), (
                f"perm{seed}: --prefilter-min-cells changed {col} -- it must "
                "restrict WHICH PAS are tested, never the test itself "
                "(the within-gene denominator comes from the full matrix)"
            )


# ---------------------------------------------------------------------------
# 3. It cuts work, and OFF changes nothing
# ---------------------------------------------------------------------------


def test_prefilter_reduces_the_number_of_tested_pas(null_runs):
    """The speed knob has to actually shrink the search space."""
    for seed in PERM_SEEDS:
        n_off = len(null_runs["off"][seed])
        n_pre = len(null_runs["prefilter"][seed])
        assert n_pre < n_off, (
            f"perm{seed}: --prefilter-min-cells {PREFILTER_MIN_CELLS} tested "
            f"{n_pre} PAS vs {n_off} unfiltered -- it saved nothing"
        )
        assert n_pre <= 0.5 * n_off, (
            f"perm{seed}: expected a substantial cut, got {n_pre}/{n_off}"
        )


def test_prefilter_off_is_byte_identical(tmp_path):
    """Explicit 0 and the flag never passed must produce the same frame."""
    adata = _null_adata(0)
    absent = _run_diff(adata, tmp_path / "absent")[("A", "B")]
    zero = _run_diff(
        adata, tmp_path / "zero", prefilter_min_cells=0
    )[("A", "B")]
    pd.testing.assert_frame_equal(absent, zero)


# ---------------------------------------------------------------------------
# The selector in isolation
# ---------------------------------------------------------------------------


def test_select_expressed_pas_counts_cells_not_reads():
    """Detection, not depth: one 100-count cell is one cell."""
    m = pd.DataFrame(
        {
            "deep_but_rare": [100.0, 0.0, 0.0, 0.0],
            "broad": [1.0, 1.0, 1.0, 0.0],
            "everywhere": [1.0, 2.0, 1.0, 3.0],
        },
        index=[f"c{i}" for i in range(4)],
    )
    assert select_expressed_pas(m, min_cells=3) == ["broad", "everywhere"]
    assert select_expressed_pas(m, min_cells=4) == ["everywhere"]
    assert select_expressed_pas(m, min_cells=0) == list(m.columns)
    assert select_expressed_pas(m, min_cells=99) == []


# ---------------------------------------------------------------------------
# The invariance promise -- once fisher-only, now general.
#
# The docs, the CHANGELOG and two runtime warnings all used to state, without
# qualification, that a pre-selection "never changes a test" and that every
# surviving p-value is bit-identical to the unrestricted run's. That holds for
# fisher, which is handed ``full_count_matrix`` for its denominator. It is
# false for the NB strategies until 0.3.0: they accepted
# ``full_count_matrix`` via ``**_ignored`` and derived the GLM's per-cell
# library-size offset from the matrix they were actually given, so restricting
# the columns moved every offset, coefficient and p-value. They now take the
# offset from the unrestricted matrix, so the promise holds for them too --
# which is what the second test below pins.
#
# These two tests pin BOTH halves of the real behaviour. If the offset is ever
# fixed to use the full matrix (issue #124), the second test fails loudly --
# which is the point: whoever makes NB invariant must also come back and
# un-qualify the documentation this file's docstring now qualifies.
# ---------------------------------------------------------------------------

def _restriction_fixture():
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(7)
    cells = [f"c{i}" for i in range(80)]
    pas = [f"P{i}" for i in range(8)]
    mat = rng.poisson(6, size=(80, 8))
    mat[:40, 0] = rng.poisson(3, 40)
    mat[40:, 0] = rng.poisson(12, 40)
    counts = pd.DataFrame(mat, index=cells, columns=pas)
    labels = pd.Series(["A"] * 40 + ["B"] * 40, index=cells)
    kept = pas[:4]
    gene_map = {p: "G1" for p in pas}
    return counts, labels, kept, gene_map


def test_fisher_pvalues_survive_a_restriction_unchanged():
    """The documented invariance, for the one strategy that has it."""
    from ema.switch_test.strategies import get_diff_strategy

    counts, labels, kept, gene_map = _restriction_fixture()
    strat = get_diff_strategy("fisher")
    kw = dict(cluster_labels=labels, cluster1="A", cluster2="B",
              min_cells_per_group=1, pas_gene_map=gene_map)
    full = strat.test(count_matrix=counts, **kw)
    restricted = strat.test(count_matrix=counts[kept],
                            full_count_matrix=counts, **kw)
    shared = [p for p in kept if p in full.index and p in restricted.index]
    assert shared, "fixture produced no comparable rows"
    for pas_id in shared:
        assert full.loc[pas_id, "pvalue"] == restricted.loc[pas_id, "pvalue"], (
            f"fisher p-value for {pas_id} changed under a restriction; the "
            "full-matrix denominator fix (issue #94) has regressed"
        )


def test_nb_pairwise_pvalues_survive_a_restriction_unchanged():
    """The NB offset fix: invariance now holds for nb_pairwise too.

    Until 0.3.0 nb_pairwise took ``full_count_matrix`` via ``**_ignored`` and
    built its GLM offset from ``mat_sub.values.sum(axis=1)`` -- the library
    size of whatever matrix it was handed. A pre-selection therefore moved
    every cell's offset and every p-value with it (up to two orders of
    magnitude). The offset is now read from the unrestricted matrix, because a
    cell's sequencing depth cannot depend on which hypotheses were selected.
    """
    from ema.switch_test.strategies import get_diff_strategy

    counts, labels, kept, gene_map = _restriction_fixture()
    strat = get_diff_strategy("nb_pairwise")
    kw = dict(cluster_labels=labels, cluster1="A", cluster2="B",
              min_cells_per_group=1, pas_gene_map=gene_map)
    full = strat.test(count_matrix=counts, **kw)
    restricted = strat.test(count_matrix=counts[kept],
                            full_count_matrix=counts, **kw)
    shared = [p for p in kept if p in full.index and p in restricted.index]
    assert shared, "fixture produced no comparable rows"
    for pas_id in shared:
        assert full.loc[pas_id, "pvalue"] == restricted.loc[pas_id, "pvalue"], (
            f"nb_pairwise p-value for {pas_id} changed under a restriction: "
            f"{full.loc[pas_id, 'pvalue']!r} -> "
            f"{restricted.loc[pas_id, 'pvalue']!r}. The library-size offset "
            "must come from full_count_matrix, not the restricted matrix."
        )


def test_nb_pairwise_without_the_full_matrix_still_shifts():
    """Guard the mechanism, not just the outcome.

    If someone removes the ``full_count_matrix`` plumbing in ``runner.py`` the
    strategy still *accepts* the argument, so the test above would keep passing
    while real runs silently regressed. This pins that the offset genuinely
    comes from the argument: withhold it and the old behaviour returns.
    """
    from ema.switch_test.strategies import get_diff_strategy

    counts, labels, kept, gene_map = _restriction_fixture()
    strat = get_diff_strategy("nb_pairwise")
    kw = dict(cluster_labels=labels, cluster1="A", cluster2="B",
              min_cells_per_group=1, pas_gene_map=gene_map)
    full = strat.test(count_matrix=counts, **kw)
    no_full = strat.test(count_matrix=counts[kept], **kw)
    shared = [p for p in kept if p in full.index and p in no_full.index]
    changed = [p for p in shared
               if full.loc[p, "pvalue"] != no_full.loc[p, "pvalue"]]
    assert changed, (
        "restricting the matrix WITHOUT passing full_count_matrix no longer "
        "changes the p-values -- the offset is evidently not being taken from "
        "the count matrix at all, so this test no longer guards the plumbing."
    )
