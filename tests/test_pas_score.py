"""TASK D -- the calibrated per-site PAS score (``--pas-score``).

What these tests protect, in order of what it would cost to get wrong:

1. **v2 stays reachable.**  ``--pas-score none`` is the branch default and is
   v2: no model is loaded, no column is written, no PAS moves.  The whole-file
   guarantee is ``test_prime_v2_compat_golden.py``; here we pin the narrower,
   more diagnostic statement.
2. **The shipped constants are the fitted model.**  The artefact in
   ``ema/countmatrix/models/`` is a table of numbers, not code, and nothing in
   the test suite can re-fit it -- so a golden probability is pinned on a fixed
   feature row.  If the JSON is ever regenerated, that number moves and this
   test says so.
3. **scikit-learn never enters the import path.**  ``PRIME_PLAN.md``
   non-negotiable 3.  Asserted by importing the scoring module and looking at
   ``sys.modules``, not by reading the source.
4. **The hard gates stay hard.**  The score is a re-ranker inside tier-1 and
   inside the internal-priming veto: ``select`` may only REMOVE a tier-1
   candidate.  It must never promote a tier-2 one and never rescue one the veto
   dropped.
5. **"unscoreable" and "scored badly" are different words.**  A candidate whose
   sequence window could not be read gets ``NA`` and is exempt from the gate.
6. **The model's feature list is a subset of the sidecar's header, verbatim.**
   That is the invariant that makes the offline fitter and the tool the same
   computation; a derived feature would be a place for them to drift.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest

from ema.config import variable_config
from ema.countmatrix import pas_score as ps
from ema.countmatrix.pas_features import SEAM_FEATURE_COLUMNS, FeatureCollector
from ema.countmatrix.paswrite import CALL_FEATURE_COLUMNS, SUPPORT_COLUMNS

SIDECAR_COLUMNS = tuple(SUPPORT_COLUMNS) + tuple(CALL_FEATURE_COLUMNS) + tuple(
    SEAM_FEATURE_COLUMNS)


@pytest.fixture(autouse=True)
def _restore_variable_config():
    keys = ("pas_features", "pas_score", "pas_score_model", "pas_score_min")
    saved = {k: getattr(variable_config, k) for k in keys}
    try:
        yield
    finally:
        for k, v in saved.items():
            setattr(variable_config, k, v)


# ---------------------------------------------------------------------------
# transforms
# ---------------------------------------------------------------------------
def test_rank_pct_is_the_midrank_percentile():
    v = np.array([1.0, 1.0, 2.0, 3.0, 3.0, 3.0, 5.0])
    got = ps.rank_pct(v)
    # ties share the midpoint of the ranks they span
    assert got[0] == got[1] == pytest.approx(1.0 / 7)
    assert got[2] == pytest.approx(2.5 / 7)
    assert got[3] == got[4] == got[5] == pytest.approx(4.5 / 7)
    assert got[6] == pytest.approx(6.5 / 7)


def test_rank_pct_does_not_depend_on_input_order():
    rng = np.random.default_rng(0)
    v = rng.integers(0, 20, 500).astype(float)
    perm = rng.permutation(len(v))
    a = ps.rank_pct(v)[perm]
    b = ps.rank_pct(v[perm])
    assert np.array_equal(a, b)


def test_rank_pct_of_an_empty_or_constant_column():
    assert ps.rank_pct(np.array([])).shape == (0,)
    assert np.allclose(ps.rank_pct(np.zeros(9)), 0.5)


@pytest.mark.parametrize("kind,x,want", [
    ("raw", -3.0, -3.0),
    ("log1p", -3.0, 0.0),                       # negatives clamp to 0 first
    ("log1p", np.e - 1, 1.0),
    ("signlog1p", -(np.e - 1), -1.0),
    ("poslog1p", -5.0, 0.0),
    ("poslog1p", np.e - 1, 1.0),
])
def test_each_transform_does_what_its_name_says(kind, x, want):
    assert ps.apply_transform(np.array([x]), kind)[0] == pytest.approx(want)


def test_nan_becomes_zero_rather_than_propagating():
    assert ps.apply_transform(np.array([np.nan]), "raw")[0] == 0.0


def test_an_unknown_transform_raises_rather_than_passing_through():
    with pytest.raises(ValueError, match="unknown transform"):
        ps.apply_transform(np.array([1.0]), "sqrt")


# ---------------------------------------------------------------------------
# the two model families
# ---------------------------------------------------------------------------
def _linear_spec():
    return {"family": "linear", "name": "t", "features": ["a", "b"],
            "transform": ["raw", "raw"], "coef": [2.0, -1.0], "intercept": 0.5,
            "threshold": 0.5}


def test_a_linear_model_is_a_dot_product_and_a_sigmoid():
    m = ps.PasScoreModel(_linear_spec())
    p = m.predict_proba({"a": np.array([1.0, 0.0]), "b": np.array([0.0, 1.0])})
    assert p[0] == pytest.approx(1 / (1 + np.exp(-2.5)))
    assert p[1] == pytest.approx(1 / (1 + np.exp(0.5)))


def _stump_spec():
    """One tree: `a <= 1.5` -> +1, else -1, on a baseline of 0."""
    return {"family": "hgb", "name": "t", "features": ["a"], "transform": ["raw"],
            "baseline": 0.0,
            "node_feature": [0, 0, 0], "node_threshold": [1.5, 0.0, 0.0],
            "node_left": [1, 0, 0], "node_right": [2, 0, 0],
            "node_is_leaf": [0, 1, 1], "node_value": [0.0, 1.0, -1.0],
            "tree_offset": [0, 3], "threshold": 0.5}


def test_a_tree_stump_splits_where_it_says_it_does():
    m = ps.PasScoreModel(_stump_spec())
    p = m.predict_proba({"a": np.array([1.0, 1.5, 2.0])})
    assert p[0] == pytest.approx(1 / (1 + np.exp(-1.0)))
    assert p[1] == pytest.approx(1 / (1 + np.exp(-1.0)))     # <= is inclusive
    assert p[2] == pytest.approx(1 / (1 + np.exp(1.0)))


def test_traversal_is_chunk_size_invariant(monkeypatch):
    m = ps.PasScoreModel(_stump_spec())
    x = np.linspace(0, 3, 1000)
    full = m.predict_proba({"a": x})
    monkeypatch.setattr(ps.PasScoreModel, "CHUNK", 7)
    assert np.array_equal(full, m.predict_proba({"a": x}))


def test_a_cycle_in_the_node_arrays_fails_loudly_instead_of_hanging():
    spec = _stump_spec()
    spec["node_left"] = [1, 0, 0]
    spec["node_is_leaf"] = [0, 0, 1]        # node 1 points back at the root
    spec["node_right"] = [2, 0, 0]
    with pytest.raises(ValueError, match="did not reach a leaf"):
        ps.PasScoreModel(spec).predict_proba({"a": np.array([0.0])})


def test_a_missing_feature_column_raises_instead_of_scoring_a_zero():
    m = ps.PasScoreModel(_linear_spec())
    with pytest.raises(KeyError, match="missing feature column"):
        m.predict_proba({"a": np.array([1.0])})


def test_an_unknown_family_is_refused():
    spec = _linear_spec()
    spec["family"] = "randomforest"
    with pytest.raises(ValueError, match="unknown family"):
        ps.PasScoreModel(spec)


# ---------------------------------------------------------------------------
# the shipped model
# ---------------------------------------------------------------------------
def test_the_shipped_model_loads_and_declares_its_provenance():
    m = ps.load_model("prime1")
    assert m.name == "prime1"
    assert m.spec["trained_on"] == "mouse1"
    assert m.spec["label"].startswith("PolyASite")
    assert 0.0 <= m.threshold <= 1.0
    assert m.spec["threshold_rule"].startswith("T1") or \
        m.spec["threshold_rule"].startswith("T2")


def test_every_feature_the_shipped_model_names_is_a_sidecar_column():
    """The invariant that makes the offline fitter and the tool one computation.

    ``scripts/prime/taskD_fit_model.py`` builds its design matrix from
    ``pas_support.tsv``; the seam builds it from the caller's columns plus the
    feature collector.  They are the same numbers only while every feature name
    is a sidecar column VERBATIM -- no derived quantities, no renames.
    """
    m = ps.load_model("prime1")
    unknown = [f for f in m.features if f not in SIDECAR_COLUMNS]
    assert unknown == [], f"model features not in pas_support.tsv: {unknown}"


def test_the_shipped_model_reproduces_a_pinned_probability():
    """A golden on the CONSTANTS, so re-fitting cannot land silently.

    The row is deliberately synthetic and unremarkable; the number is whatever
    the shipped JSON returns for it.  If this fails, either the model file was
    regenerated (say so in the commit and update the constant) or the evaluator
    changed (that is a bug).
    """
    m = ps.load_model("prime1")
    row = {f: np.array([_PINNED_ROW[f]], dtype=float) for f in m.features}
    assert float(m.predict_proba(row)[0]) == pytest.approx(_PINNED_PROB, abs=1e-9)


#: One synthetic candidate: a 3-molecule tier-1 site with a strong hexamer.
_PINNED_ROW = {
    "clip_reads": 5, "clip_umis": 3, "clip_reads_f3844": 5, "clip_umis_f3844": 3,
    "tier": 1, "clip_positions": 2, "clip_span": 7, "window_reads": 120,
    "hex_strong": 1, "hex_any12": 1, "hex_n_types": 1, "hex_best_off": -21,
    "hex_strong_off": -21, "seq_ok": 1, "d_prev_cand": 900, "d_next_cand": 1400,
    "n_cand_100": 0, "n_cand_500": 1, "mol_500_sum": 2, "is_local_mol_max": 1,
    "mol_frac_local": 0.6,
}
#: What ema/countmatrix/models/pas_score_model_prime1.json returns for it.
_PINNED_PROB = 0.7208753652231465


def test_scoring_never_imports_scikit_learn(tmp_path):
    """PRIME_PLAN non-negotiable 3, asserted rather than asserted-in-prose.

    In a SUBPROCESS, because ``ema.main`` pulls in scanpy (which imports
    scikit-learn) for the clustering stage, so an in-process check would pass
    or fail depending on which test ran first.  The claim being pinned is the
    one that matters: **loading a model and scoring with it** must not need
    scikit-learn, so the fitted model really is constants and a deployment that
    strips scanpy can still score.
    """
    import subprocess

    prog = (
        "import sys, numpy as np\n"
        "from ema.countmatrix.pas_score import load_model\n"
        "m = load_model('prime1')\n"
        "m.predict_proba({f: np.zeros(3) for f in m.features})\n"
        "bad = [k for k in sys.modules if k == 'sklearn' or k.startswith('sklearn.')]\n"
        "print('SKLEARN' if bad else 'CLEAN')\n"
    )
    root = Path(__file__).resolve().parents[1]
    env = dict(os.environ, PYTHONPATH=str(root))
    out = subprocess.run([sys.executable, "-c", prog], capture_output=True,
                         text=True, env=env, cwd=str(root))
    assert out.returncode == 0, out.stderr[-2000:]
    assert out.stdout.strip() == "CLEAN", (
        "scikit-learn reached the scoring import path: " + out.stdout)


# ---------------------------------------------------------------------------
# the sidecar view
# ---------------------------------------------------------------------------
def test_scored_features_appends_exactly_one_field():
    rows = {"1": "a\tb", "2": "c\td"}
    v = ps.ScoredFeatures(rows, {"1": 0.25})
    assert v.get("1") == "a\tb\t0.250000"
    assert v.get("2") == "c\td\tNA"          # no score -> NA, never 0
    assert v.get("3") is None
    assert sorted(v) == ["1", "2"] and len(v) == 2


def test_scored_features_does_not_mutate_the_collectors_rows():
    rows = {"1": "a\tb"}
    v = ps.ScoredFeatures(rows, {"1": 0.5})
    v.get("1"), v.get("1")
    assert rows["1"] == "a\tb"


def test_a_nan_score_is_NA_and_not_a_zero():
    assert ps.ScoredFeatures.format(float("nan")) == "NA"
    assert ps.ScoredFeatures.format(None) == "NA"
    assert ps.ScoredFeatures.format(0.0) == "0.000000"


# ---------------------------------------------------------------------------
# defaults and configuration
# ---------------------------------------------------------------------------
def test_the_branch_default_is_single_valued():
    """One default, one place -- the defect TASK A shipped a fix for."""
    from ema.cli.config_schema import RunConfig

    rc = RunConfig()
    assert variable_config.pas_score == "none"
    assert rc.pas_score == variable_config.pas_score
    assert rc.pas_score_model == variable_config.pas_score_model
    assert rc.pas_score_min == variable_config.pas_score_min


def test_the_flag_is_validated_not_guessed(tmp_path):
    from ema.config import args
    from ema.main import _validate_pas_score_config

    ns = args._get()
    variable_config.pas_score = "calibrated"
    variable_config.pas_features = "off"
    try:
        with pytest.raises(ValueError, match="requires --pas-features on"):
            _validate_pas_score_config()
        variable_config.pas_features = "on"
        ns.genome_fasta = None
        with pytest.raises(ValueError, match="requires --genome-fasta"):
            _validate_pas_score_config()
    finally:
        try:
            delattr(ns, "genome_fasta")
        except AttributeError:
            pass


def test_pas_score_none_validates_without_a_genome():
    from ema.main import _validate_pas_score_config

    variable_config.pas_score = "none"
    variable_config.pas_features = "off"
    _validate_pas_score_config()             # must not raise


def test_an_unknown_model_name_names_where_it_looked():
    with pytest.raises(FileNotFoundError, match="pas_score_model_nope.json"):
        ps.load_model("nope")


# ---------------------------------------------------------------------------
# the ema.main seam, end to end
# ---------------------------------------------------------------------------
_SEAM_ARG_NAMES = ("ip_filter", "ip_filter_mode", "annot_filter", "genome_fasta",
                   "annotation_bed", "ip_a_stretch", "ip_a_fraction",
                   "ip_window_left", "ip_window_right", "polya_evidence",
                   "polya_mode")


@pytest.fixture
def seam(tmp_path):
    """A run dir with pos/neg PAS BEDs, a full sidecar and a synthetic genome.

    Five candidates, chosen to exercise every branch the gate has:
      1  tier-1, 5 molecules, plain sequence
      2  tier-1, 1 molecule
      3  tier-1, 7 molecules, on '-'
      4  **tier-2** (BED score 0): the score must never promote it and never
         drop it -- it is outside the gate by construction
      5  on a contig the FASTA does not have: unscoreable, so ``NA`` and exempt
      6  tier-1 with an 8-A run immediately downstream: the internal-priming
         veto's own row, used to prove the score cannot rescue it
    """
    pytest.importorskip("pyfaidx")
    from ema.config import args, directory_config, set_directory_config
    from ema.outputs import OutputManager

    ns = args._get()
    for name in _SEAM_ARG_NAMES:
        try:
            delattr(ns, name)
        except AttributeError:
            pass

    run_dir = tmp_path / "run"
    set_directory_config(output_dir=run_dir, gtf_dir=None)
    mgr = OutputManager(base_dir=str(run_dir))
    mgr.setup()

    seq = list(("CG" * 300)[:600])
    seq[320:326] = list("TTTATT")            # revcomp(AATAAA) upstream of PAS 3
    seq[510:518] = list("A" * 8)             # 8-A run downstream of PAS 6
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\n" + "".join(seq) + "\n")

    Path(directory_config.posbed).write_text(
        "chr1\t99\t100\t1\t5\t+\n"
        "chr1\t199\t200\t2\t1\t+\n"
        "chr1\t399\t400\t4\t0\t+\n"
        "chrZ\t50\t51\t5\t4\t+\n"
        "chr1\t499\t500\t6\t9\t+\n"
    )
    Path(directory_config.negbed).write_text("chr1\t300\t301\t3\t7\t-\n")
    support = run_dir / "pas_support.tsv"
    support.write_text(
        "\t".join(SUPPORT_COLUMNS + CALL_FEATURE_COLUMNS) + "\n"
        "1\t9\t5\t9\t5\t120\t1\t2\t7\n"
        "2\t1\t1\t1\t1\t44\t1\t1\t0\n"
        "3\t14\t7\t14\t7\t260\t1\t3\t18\n"
        "4\t0\t0\t0\t0\t61\t2\t0\t0\n"
        "5\t6\t4\t6\t4\t95\t1\t2\t4\n"
        "6\t18\t9\t18\t9\t310\t1\t4\t22\n"
    )
    variable_config.pas_features = "on"
    variable_config.pas_score = "none"
    variable_config.pas_score_model = "prime1"
    variable_config.pas_score_min = -1.0
    yield {"mgr": mgr, "run_dir": run_dir, "genome": genome, "support": support}
    for name in _SEAM_ARG_NAMES:
        try:
            delattr(ns, name)
        except AttributeError:
            pass


def _seam_rows(support: Path):
    lines = support.read_text().splitlines()
    cols = lines[0].split("\t")
    return cols, {ln.split("\t")[0]: dict(zip(cols, ln.split("\t")))
                  for ln in lines[1:]}


def _run_seam(seam, mode, ip_filter=False, score_min=-1.0):
    from ema.config import args, directory_config
    from ema.main import _apply_pas_filters, _write_pas_features

    ns = args._get()
    ns.genome_fasta = str(seam["genome"])
    ns.ip_filter = ip_filter
    ns.ip_filter_mode = "filter"
    variable_config.pas_score = mode
    variable_config.pas_score_min = score_min
    res = _apply_pas_filters(seam["mgr"])
    _write_pas_features(res)
    beds = (Path(directory_config.posbed).read_text()
            + Path(directory_config.negbed).read_text())
    return res, beds, _seam_rows(seam["support"])


def test_calibrated_writes_a_probability_and_moves_no_call(seam):
    from ema.config import directory_config

    before = (Path(directory_config.posbed).read_bytes()
              + Path(directory_config.negbed).read_bytes())
    res, beds, (cols, rows) = _run_seam(seam, "calibrated")
    assert beds.encode() == before, "--pas-score calibrated changed the BED"
    assert cols[-1] == "pas_score"
    assert cols[:-1] == list(SUPPORT_COLUMNS + CALL_FEATURE_COLUMNS
                             + SEAM_FEATURE_COLUMNS)
    for pid in ("1", "2", "3", "4", "6"):
        assert 0.0 <= float(rows[pid]["pas_score"]) <= 1.0
    st = res["pas_score_stats"]
    assert st["mode"] == "calibrated" and st["n_scored"] == 5


def test_none_writes_no_score_column(seam):
    _res, _beds, (cols, _rows) = _run_seam(seam, "none")
    assert "pas_score" not in cols


def test_an_unscoreable_candidate_is_NA_and_survives_the_gate(seam):
    _res, beds, (_cols, rows) = _run_seam(seam, "select", score_min=0.99)
    # PAS 5 sits on a contig the FASTA does not have: seq_ok 0, so it cannot be
    # scored -- and "we could not score it" must not be spelled like "it lost".
    assert rows["5"]["seq_ok"] == "0"
    assert rows["5"]["pas_score"] == "NA"
    assert "chrZ\t50\t51\t5\t4\t+" in beds


def test_select_can_only_remove_tier_one_candidates(seam):
    _res, beds, (_cols, rows) = _run_seam(seam, "select", score_min=0.999)
    kept = {ln.split("\t")[3] for ln in beds.splitlines() if ln.strip()}
    # tier-2 candidate 4 is outside the gate: still there at any threshold
    assert "4" in kept
    # ...and every tier-1 candidate that could be scored is gone at p >= 0.999
    for pid in ("1", "2", "3", "6"):
        assert float(rows[pid]["pas_score"]) < 0.999
        assert pid not in kept


def test_select_at_a_zero_threshold_drops_nothing(seam):
    from ema.config import directory_config

    before = (Path(directory_config.posbed).read_text()
              + Path(directory_config.negbed).read_text())
    _res, beds, _ = _run_seam(seam, "select", score_min=0.0)
    assert beds == before


def test_the_score_cannot_rescue_an_internally_primed_candidate(seam):
    """The veto is a HARD GATE in front of the score (VERIFY §1.4, §7.1)."""
    _res, beds, (_cols, rows) = _run_seam(seam, "select", ip_filter=True,
                                          score_min=0.0)
    kept = {ln.split("\t")[3] for ln in beds.splitlines() if ln.strip()}
    assert rows["6"]["ip_tool_flag"] == "1"     # the 8-A run downstream
    assert "6" not in kept, "an internally-primed PAS survived a 0.0 threshold"
    assert "1" in kept


def test_select_records_what_it_did(seam):
    res, _beds, _ = _run_seam(seam, "select", score_min=0.999)
    st = res["pas_score_stats"]
    assert st["mode"] == "select" and st["threshold"] == pytest.approx(0.999)
    assert st["n_dropped"] == 4 and st["n_unscoreable"] == 1
    stats = json.loads((seam["run_dir"] / "04_pas_gene_assignment"
                        / "pas_score_stats.json").read_text())
    assert stats["n_dropped"] == 4


def test_the_default_threshold_comes_from_the_model(seam):
    res, _beds, _ = _run_seam(seam, "calibrated")
    assert res["pas_score_stats"]["threshold"] == pytest.approx(
        ps.load_model("prime1").threshold)


def test_select_refuses_to_run_without_the_caller_columns(seam):
    """A gate that silently does not apply is worse than one that fails.

    Multi-BAM runs re-key their PAS ids at merge time and have no run-root
    ``pas_support.tsv``; `calibrated` then skips with a warning, but `select`
    must refuse rather than emit a call set the run record says was filtered.
    """
    from ema.config import args
    from ema.main import _apply_pas_filters

    seam["support"].unlink()
    ns = args._get()
    ns.genome_fasta = str(seam["genome"])
    ns.ip_filter = False
    variable_config.pas_score = "select"
    with pytest.raises(ValueError, match="never selected"):
        _apply_pas_filters(seam["mgr"])
    variable_config.pas_score = "calibrated"
    res = _apply_pas_filters(seam["mgr"])          # warns, does not raise
    assert "pas_score" not in (res or {})


def test_a_within_run_transform_is_computed_over_the_whole_run(seam):
    """`rankpct` must see every candidate, not only the scoreable ones.

    The offline fitter computes the percentile over the run's whole
    ``pas_support.tsv``; if the seam computed it over the subset it can score,
    the two would be different numbers with one name.  PAS 5 is on a contig the
    FASTA does not have, so it is unscoreable -- and it must still count towards
    everyone else's percentile.
    """
    from ema.config import args
    from ema.main import _apply_pas_filters

    spec = {"family": "linear", "name": "rank", "features": ["clip_reads"],
            "transform": ["rankpct"], "coef": [8.0], "intercept": -4.0,
            "threshold": 0.5}
    (Path(ps.MODEL_DIR) / "pas_score_model_ranktest.json").write_text(
        json.dumps(spec) + "\n")
    try:
        ns = args._get()
        ns.genome_fasta = str(seam["genome"])
        ns.ip_filter = False
        variable_config.pas_score = "calibrated"
        variable_config.pas_score_model = "ranktest"
        res = _apply_pas_filters(seam["mgr"])
        got = res["pas_score"]
        # sidecar clip_reads: 9, 1, 14, 0, 6, 18 over SIX rows.  PAS 1's value 9
        # is the 4th smallest of six -> midrank percentile 3.5/6, which is only
        # true if PAS 5 (unscoreable, clip_reads 6) was counted.
        assert got["1"] == pytest.approx(1 / (1 + np.exp(-(8 * 3.5 / 6 - 4))))
        assert np.isnan(got["5"])
    finally:
        (Path(ps.MODEL_DIR) / "pas_score_model_ranktest.json").unlink()
