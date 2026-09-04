"""Per-site calibrated PAS score (peakAtail-prime, ``--pas-score``).

WHAT THIS IS
------------
`results/algo_headroom/VERIFY/VERIFICATION_REPORT.md` §1.5/§7.2 found exactly one
change that LIFTS the caller's precision/recall curve rather than sliding along
it (besides the internal-priming veto, which the tool already ships): a
calibrated per-site score used as a **re-ranker inside the existing gates**,
replacing only the ">= 2 clip molecules" threshold.  This module evaluates that
score.  It decides nothing on its own: tier-1 membership and the internal-priming
veto stay HARD GATES in front of it, because every configuration in which a
score was allowed to override the veto looked spectacular on the atlas and no
better on long reads [V §1.4, §7.1].

NO SCIKIT-LEARN AT RUN TIME
---------------------------
`PRIME_PLAN.md` non-negotiable 3 and `manuscript/24_prime_preregistration.md`
§3.3: models are fitted OFFLINE (`scripts/prime/taskD_fit_model.py`) and what
ships is a table of constants evaluated here with numpy alone.  ``ema`` must
never import scikit-learn.  The offline fitter asserts that this module
reproduces scikit-learn's own probabilities to 1e-9 on every training row before
it writes a model out, so "constants" is a checked claim rather than a hope.

THE FEATURE VECTOR IS THE SIDECAR
---------------------------------
Every feature the shipped model uses is a column of ``pas_support.tsv`` -- v2's
seven plus the ``--pas-features on`` columns.  That is deliberate: the score can
therefore never depend on something the caller does not itself write, and any
score in the sidecar can be recomputed from the row next to it.

TRANSFORMS
----------
Four kinds, named in the model spec so the fitter and the tool cannot drift:

``raw``        the value
``log1p``      log(1 + max(v, 0))                 -- counts
``signlog1p``  sign(v) * log(1 + |v|)             -- signed offsets
``poslog1p``   log(1 + clip(v, 0, 1e6))           -- one-sided distances
``rankpct``    the value's mid-rank percentile WITHIN THIS RUN, in [0, 1]

``rankpct`` exists because six covariates (clip_reads, clip_umis, their -F 3844
twins, window_reads, mol_500_sum) carry the library's sequencing depth in their
units, and a model fitted on one library would otherwise learn thresholds in
those units.  It is a within-run statistic over the candidates the seam sees, so
it is computable exactly once, at the seam, from data already in hand.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np

#: Where a model shipped with the package lives.
MODEL_DIR = Path(__file__).resolve().parent / "models"

TRANSFORMS = ("raw", "log1p", "signlog1p", "poslog1p", "rankpct")


def rank_pct(v: np.ndarray) -> np.ndarray:
    """Mid-rank percentile of *v* among its own values, in [0, 1].

    ``(n_below + n_below_or_equal) / (2 n)``.  Ties get the midpoint, so the
    result does not depend on input order.  Three searchsorted lines rather than
    an argsort round-trip so it is O(n log n) with no Python-level loop; the
    offline fitter uses the identical definition
    (``scripts/prime/taskD_transfer.py::rank_pct``).
    """
    n = v.shape[0]
    if n == 0:
        return v.astype(np.float64)
    s = np.sort(v)
    lo = np.searchsorted(s, v, side="left")
    hi = np.searchsorted(s, v, side="right")
    return (lo + hi) / (2.0 * n)


def apply_transform(v: np.ndarray, kind: str) -> np.ndarray:
    v = np.asarray(v, dtype=np.float64)
    v = np.nan_to_num(v, nan=0.0, posinf=0.0, neginf=0.0)
    if kind == "raw":
        return v
    if kind == "log1p":
        return np.log1p(np.maximum(v, 0.0))
    if kind == "signlog1p":
        return np.sign(v) * np.log1p(np.abs(v))
    if kind == "poslog1p":
        return np.log1p(np.minimum(np.maximum(v, 0.0), 1e6))
    if kind == "rankpct":
        return rank_pct(v)
    raise ValueError(f"unknown transform {kind!r}; known: {TRANSFORMS}")


def _sigmoid(z: np.ndarray) -> np.ndarray:
    # branch-free and overflow-free: exp is only ever applied to <= 0
    out = np.empty_like(z)
    pos = z >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-z[pos]))
    ez = np.exp(z[~pos])
    out[~pos] = ez / (1.0 + ez)
    return out


class PasScoreModel:
    """A fitted per-site score, evaluated with numpy only.

    Two families ship-able:

    ``linear``  intercept + coefficients, one dot product.  ~20 numbers a
                reviewer can read, which is why it is the default family.
    ``hgb``     a flattened gradient-boosted tree ensemble (node arrays).
                Supported so the offline comparison between the two is
                measurable rather than asserted, and so a heavier model can be
                supplied with ``--pas-score-model PATH`` without a code change.
    """

    __slots__ = ("spec", "features", "transform", "family", "_coef", "_intercept",
                 "_nf", "_nt", "_nl", "_nr", "_nleaf", "_nv", "_toff", "_baseline")

    def __init__(self, spec: dict):
        self.spec = spec
        self.family = spec["family"]
        self.features = list(spec["features"])
        self.transform = list(spec["transform"])
        if len(self.features) != len(self.transform):
            raise ValueError("model spec: features and transform differ in length")
        for t in self.transform:
            if t not in TRANSFORMS:
                raise ValueError(f"model spec: unknown transform {t!r}")
        if self.family == "linear":
            self._coef = np.asarray(spec["coef"], dtype=np.float64)
            self._intercept = float(spec["intercept"])
            if self._coef.shape[0] != len(self.features):
                raise ValueError("model spec: coef and features differ in length")
        elif self.family == "hgb":
            self._nf = np.asarray(spec["node_feature"], dtype=np.int32)
            self._nt = np.asarray(spec["node_threshold"], dtype=np.float64)
            self._nl = np.asarray(spec["node_left"], dtype=np.int32)
            self._nr = np.asarray(spec["node_right"], dtype=np.int32)
            self._nleaf = np.asarray(spec["node_is_leaf"], dtype=bool)
            self._nv = np.asarray(spec["node_value"], dtype=np.float64)
            self._toff = np.asarray(spec["tree_offset"], dtype=np.int64)
            self._baseline = float(spec["baseline"])
        else:
            raise ValueError(f"model spec: unknown family {self.family!r}")

    # -- properties ---------------------------------------------------------
    @property
    def name(self) -> str:
        return str(self.spec.get("name", "unnamed"))

    @property
    def threshold(self) -> float:
        return float(self.spec.get("threshold", float("nan")))

    def needs_rank(self) -> bool:
        return "rankpct" in self.transform

    # -- evaluation ---------------------------------------------------------
    def build_matrix(self, columns: dict) -> np.ndarray:
        """Assemble the model's design matrix from ``{name: array}``.

        Every feature the model names must be present; a missing one is an
        error, never a silent zero -- a score computed from a column the run did
        not write would be indistinguishable from a real one in the sidecar.
        """
        missing = [f for f in self.features if f not in columns]
        if missing:
            raise KeyError(f"pas_score: missing feature column(s): {missing}")
        n = len(np.asarray(columns[self.features[0]]))
        X = np.empty((n, len(self.features)), dtype=np.float64)
        for j, (f, t) in enumerate(zip(self.features, self.transform)):
            X[:, j] = apply_transform(np.asarray(columns[f], dtype=np.float64), t)
        return X

    #: Rows per traversal chunk.  The ensemble is walked chunk by chunk so the
    #: per-row node index and its gathers stay in cache: 65,536 measured 6.4 s
    #: against 8.1 s unchunked over 195,940 x 21 on 161 trees / 9,821 nodes.
    CHUNK = 65_536

    def raw_matrix(self, X: np.ndarray) -> np.ndarray:
        """Raw (logit-scale) prediction for a pre-transformed design matrix."""
        X = np.ascontiguousarray(X, dtype=np.float64)
        if self.family == "linear":
            return X @ self._coef + self._intercept
        n = X.shape[0]
        out = np.empty(n, dtype=np.float64)
        for s0 in range(0, n, self.CHUNK):
            Xc = X[s0:s0 + self.CHUNK]
            nc = Xc.shape[0]
            rows = np.arange(nc)
            raw = np.full(nc, self._baseline, dtype=np.float64)
            for t in range(self._toff.shape[0] - 1):
                base = int(self._toff[t])
                end = int(self._toff[t + 1])
                nf = self._nf[base:end]
                nt = self._nt[base:end]
                nl = self._nl[base:end]
                nr = self._nr[base:end]
                leaf = self._nleaf[base:end]
                node = np.zeros(nc, dtype=np.int32)
                active = ~leaf[node]
                # bounded by the tree's depth; max_leaf_nodes=31 caps it at 30,
                # and the guard makes a malformed spec fail loudly, not hang.
                for _ in range(1024):
                    if not active.any():
                        break
                    idx = node[active]
                    go_left = Xc[rows[active], nf[idx]] <= nt[idx]
                    node[active] = np.where(go_left, nl[idx], nr[idx])
                    active = ~leaf[node]
                else:
                    raise ValueError("pas_score: tree traversal did not reach a "
                                     "leaf in 1024 steps; the model spec is "
                                     "malformed (cycle in the node arrays?)")
                raw += self._nv[base:end][node]
            out[s0:s0 + self.CHUNK] = raw
        return out

    def predict_proba_matrix(self, X: np.ndarray) -> np.ndarray:
        return _sigmoid(self.raw_matrix(X))

    def predict_proba(self, columns: dict) -> np.ndarray:
        return self.predict_proba_matrix(self.build_matrix(columns))


def load_model(spec_or_path) -> PasScoreModel:
    """Load a model from a name shipped with the package or from a JSON path."""
    if isinstance(spec_or_path, PasScoreModel):
        return spec_or_path
    if isinstance(spec_or_path, dict):
        return PasScoreModel(spec_or_path)
    s = str(spec_or_path)
    p = Path(s)
    if not p.exists():
        p = MODEL_DIR / f"pas_score_model_{s}.json"
    if not p.exists():
        raise FileNotFoundError(
            f"pas_score: no model {s!r}; looked for a file and for "
            f"{MODEL_DIR}/pas_score_model_{s}.json")
    return PasScoreModel(json.loads(p.read_text()))


class ScoredFeatures:
    """The collector's feature rows with ``pas_score`` appended, computed lazily.

    ``FeatureCollector.finish()`` returns ``{pas_id: tab-joined text}``.  This
    view adds one field to each row **without rewriting the dictionary**, so the
    sidecar is still appended to in a single pass and a genome-wide run never
    holds two copies of 650 k rows.  It also makes the append idempotent: the
    stored rows are never mutated, so calling the writer twice cannot produce a
    row with two score fields.

    Presents exactly the mapping surface
    :func:`ema.countmatrix.pas_features.append_columns` and
    :func:`~ema.countmatrix.pas_features.write_features_tsv` use.
    """

    __slots__ = ("_rows", "_scores")

    def __init__(self, rows: dict, scores: dict):
        self._rows = rows
        self._scores = scores

    @staticmethod
    def format(p) -> str:
        # p != p is the NaN test: "we could not score it" must not be spelled
        # the same way as "we scored it and it came out 0".
        return "NA" if p is None or p != p else "%.6f" % p

    def get(self, pas_id, default=None):
        row = self._rows.get(pas_id)
        if row is None:
            return default
        return row + "\t" + self.format(self._scores.get(pas_id))

    def __getitem__(self, pas_id):
        row = self.get(pas_id)
        if row is None:
            raise KeyError(pas_id)
        return row

    def __iter__(self):
        return iter(self._rows)

    def __len__(self):
        return len(self._rows)
