"""`ema switch length` — strand-aware proximal/distal selection + count source.

Two independent defects are covered here (plan `10_caller_fix_plan.md`, 0b/0c):

(a) STRAND.  ``ClassicPDUIStrategy`` never looks at ``strand`` itself; it trusts
    the ``rank`` field of ``pas_isoform_map`` to be proximal(1)->distal(N) in
    *transcription* order.  Three paths must build that rank strand-aware:

      * ``run_length``'s ``per_gene`` branch (fixed at 02f4153);
      * :func:`ema.switch_test.runner._build_gene_fallback_map`, which used to
        hardcode ``rank=1, total=1`` for every PAS — the map used by
        ``--isoform-agg per_isoform --utr-unmatched gene``.  With every rank
        tied, the classic strategy's stable sort fell back to insertion
        (genomic-ascending) order, inverting all 1667 `_gene_` rows shipped;
      * ``run_length`` when no ``pasbed.bed`` can be resolved — there is no
        strand-free answer, so it must REFUSE rather than rank by input order.

(b) COUNT SOURCE.  ``build_count_dfs`` read ``adata.X``.  For a
    ``clusters.h5ad`` written by the ``leiden_tfidf`` strategy that is
    ``log1p(TF * IDF * scale_factor)``, not counts.  It must now prefer
    ``layers['counts']`` and refuse a non-integral matrix.

Plus the cheap permanent guard the plan asks for: every written
``pdui_classic.tsv`` row must satisfy ``'+' => proximal_start < distal_start``
and ``'-' => proximal_start > distal_start``.
"""
from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from ema.switch_test.runner import run_length

CELLS = ["c0", "c1", "c2", "c3"]

# Two mirror-image 3-PAS genes on the same locus:
#   GENE_PLUS  (+):  pas 1 @1000 PROXIMAL ... pas 3 @3000 DISTAL
#   GENE_MINUS (-):  pas 4 @1000 DISTAL   ... pas 6 @3000 PROXIMAL
PASBED = [
    ("chr1", 1000, 1001, "1", 0, "+"),
    ("chr1", 2000, 2001, "2", 0, "+"),
    ("chr1", 3000, 3001, "3", 0, "+"),
    ("chr1", 1000, 1001, "4", 0, "-"),
    ("chr1", 2000, 2001, "5", 0, "-"),
    ("chr1", 3000, 3001, "6", 0, "-"),
]

# 90 reads on the true DISTAL PAS, 10 on the true PROXIMAL PAS of each gene,
# so the correct PDUI is 0.9 for both genes and a strand swap yields 0.1.
COUNTS = np.array([[10, 0, 90, 90, 0, 10]] * 4, dtype=np.int64)

VAR = pd.DataFrame(
    {"gene_id": ["GENE_PLUS"] * 3 + ["GENE_MINUS"] * 3},
    index=[str(i) for i in range(1, 7)],
)
OBS = pd.DataFrame({"leiden": ["0", "0", "1", "1"]}, index=CELLS)


def _write_inputs(tmp_path: Path, X=None, counts_layer=False) -> Path:
    adata = ad.AnnData(
        X=sp.csr_matrix(COUNTS.astype(np.float64) if X is None else X),
        obs=OBS.copy(), var=VAR.copy(),
    )
    if counts_layer:
        adata.layers["counts"] = sp.csr_matrix(COUNTS.astype(np.float64))
    h5 = tmp_path / "clusters.h5ad"
    adata.write_h5ad(h5)
    pd.DataFrame(PASBED).to_csv(
        tmp_path / "pasbed.bed", sep="\t", header=False, index=False,
    )
    return h5


def _run(tmp_path: Path, h5: Path, **kwargs) -> pd.DataFrame:
    out = tmp_path / kwargs.pop("outname", "out")
    run_length(
        h5ad_paths=[str(h5)], gtf=None, output_dir=str(out), cluster_pairs=None,
        cluster_key="leiden", strategy="classic", isoform_agg="per_gene",
        isoform_collapse="none", threads=1, **kwargs,
    )
    return pd.read_csv(out / "pdui_classic.tsv", sep="\t")


# --- (a) strand ------------------------------------------------------------

def test_minus_strand_proximal_distal_not_swapped(tmp_path):
    """Regression guard: minus-strand proximal = HIGHEST genomic coordinate.

    Production output written before the rank fix had 2299/2299 minus-strand
    genes with ``proximal_start < distal_start`` — i.e. every one inverted.
    """
    df = _run(tmp_path, _write_inputs(tmp_path))
    plus = df[df.gene_id == "GENE_PLUS"].iloc[0]
    minus = df[df.gene_id == "GENE_MINUS"].iloc[0]

    assert int(plus.proximal_pas_id) == 1 and int(plus.distal_pas_id) == 3
    assert plus.pdui == pytest.approx(0.9)

    assert int(minus.proximal_pas_id) == 6, "minus strand: proximal = highest coord"
    assert int(minus.distal_pas_id) == 4, "minus strand: distal = lowest coord"
    assert minus.pdui == pytest.approx(0.9), "a strand swap reports 0.1 here"

    # The augmented coordinate columns must agree with the strand convention.
    assert (minus.proximal_start > minus.distal_start)
    assert (plus.proximal_start < plus.distal_start)


def test_unresolvable_pasbed_must_not_silently_invert(tmp_path):
    """No pasbed => REFUSE. Ranking by adata.var order is a plus-strand guess.

    Was xfail(strict): `run_length` silently ranked PAS by input order and
    `ema switch length` had no --pasbed flag at all.
    """
    h5 = _write_inputs(tmp_path)
    moved = tmp_path / "elsewhere" / "pasbed.bed"
    moved.parent.mkdir()
    (tmp_path / "pasbed.bed").rename(moved)

    with pytest.raises(FileNotFoundError, match="--pasbed"):
        _run(tmp_path, h5)

    # ...and --pasbed is the documented way out, with the strand kept right.
    df = _run(tmp_path, h5, pasbed=str(moved), outname="out_explicit")
    minus = df[df.gene_id == "GENE_MINUS"].iloc[0]
    assert int(minus.proximal_pas_id) == 6 and int(minus.distal_pas_id) == 4
    assert minus.pdui == pytest.approx(0.9)


def test_per_isoform_gene_fallback_is_strand_aware():
    """The ``--utr-unmatched=gene`` fallback map must carry strand-aware ranks.

    Was xfail(strict): ``_build_gene_fallback_map`` hardcoded
    ``rank=1, total=1``, so ``ClassicPDUIStrategy``'s stable sort fell back to
    insertion order and reported the LOWEST coordinate as proximal for
    minus-strand genes too (1667/1667 shipped `_gene_` rows inverted).
    """
    from ema.quantification.strategies import get_pdui_strategy
    from ema.switch_test.runner import _build_gene_fallback_map

    adata = ad.AnnData(
        X=sp.csr_matrix(COUNTS[:, 3:].astype(np.float64)),
        obs=pd.DataFrame(index=CELLS),
        var=pd.DataFrame(
            {
                "gene_id": ["GENE_MINUS"] * 3,
                # Coordinates + strand are what make the rank meaningful;
                # run_length passes them from pasbed.bed, and the helper also
                # accepts them on .var (same three PAS as PASBED rows 4-6).
                "start": [1000, 2000, 3000],
                "strand": ["-", "-", "-"],
            },
            index=["4", "5", "6"],
        ),
    )
    fallback = _build_gene_fallback_map(adata)
    # rank 1 == proximal == the HIGHEST coordinate on the minus strand.
    assert fallback[6][0][3] == 1 and fallback[4][0][3] == 3
    assert {entry[0][4] for entry in fallback.values()} == {3}

    cm = pd.DataFrame(
        COUNTS[:, 3:].T.astype(float), index=[4, 5, 6], columns=CELLS,
    )
    df = get_pdui_strategy("classic").compute(
        count_matrix=cm, pas_isoform_map=fallback,
        aggregation="per_isoform", isoform_collapse="none", pseudocount=0.0,
    )
    row = df.iloc[0]
    assert int(row.proximal_pas_id) == 6 and int(row.distal_pas_id) == 4
    assert row.pdui == pytest.approx(0.9)


def test_gene_fallback_map_refuses_to_rank_without_coordinates():
    """No coordinates anywhere => raise, never a silent plus-strand ranking."""
    from ema.switch_test.runner import _build_gene_fallback_map

    adata = ad.AnnData(
        X=sp.csr_matrix(COUNTS[:, 3:].astype(np.float64)),
        obs=pd.DataFrame(index=CELLS),
        var=pd.DataFrame({"gene_id": ["GENE_MINUS"] * 3}, index=["4", "5", "6"]),
    )
    with pytest.raises(ValueError, match="pasbed"):
        _build_gene_fallback_map(adata)


def test_pdui_strand_convention_guard_catches_inversion():
    """The permanent guard: '+' => prox < dist, '-' => prox > dist."""
    from ema.switch_test.runner import _assert_pdui_strand_convention

    good = pd.DataFrame({
        "gene_id": ["GENE_PLUS", "GENE_MINUS"],
        "proximal_start": [1000, 3000],
        "distal_start": [3000, 1000],
        "proximal_strand": ["+", "-"],
        "distal_strand": ["+", "-"],
    })
    _assert_pdui_strand_convention(good, "ok.tsv")  # must not raise

    inverted = good.copy()
    inverted.loc[1, ["proximal_start", "distal_start"]] = [1000, 3000]
    with pytest.raises(RuntimeError, match="convention violated"):
        _assert_pdui_strand_convention(inverted, "bad.tsv")

    mixed = good.copy()
    mixed.loc[1, "distal_strand"] = "+"
    with pytest.raises(RuntimeError, match="DIFFERENT strands"):
        _assert_pdui_strand_convention(mixed, "bad.tsv")


# --- (b) count source ------------------------------------------------------

def _tfidf_like(C: np.ndarray) -> np.ndarray:
    """Signac Method 1, matching ema.clustering.strategies.leiden_tfidf."""
    tf = C / C.sum(axis=1, keepdims=True)
    idf = C.shape[0] / np.maximum((C > 0).sum(axis=0), 1)
    return np.log1p(tf * idf * 1e4)


def test_pdui_prefers_a_raw_counts_layer_over_normalised_X(tmp_path):
    """Was xfail(strict): ``build_count_dfs`` read ``.X`` unconditionally."""
    X = _tfidf_like(COUNTS.astype(np.float64))
    h5 = _write_inputs(tmp_path, X=X, counts_layer=True)
    df = _run(tmp_path, h5)
    plus = df[df.gene_id == "GENE_PLUS"].iloc[0]
    assert plus.pdui == pytest.approx(0.9)
    assert plus.proximal_reads == pytest.approx(10.0)
    assert plus.distal_reads == pytest.approx(90.0)


def test_pdui_on_tfidf_X_raises_without_the_escape_hatch(tmp_path):
    """A TF-IDF .X with no counts layer must abort, not produce a number."""
    X = _tfidf_like(COUNTS.astype(np.float64))
    h5 = _write_inputs(tmp_path, X=X)
    with pytest.raises(ValueError, match="NOT integral"):
        _run(tmp_path, h5)


def test_pdui_on_tfidf_X_is_measurably_wrong(tmp_path):
    """Documents the size of the error the guard prevents.

    Reachable only through the explicit ``--allow-non-count-matrix`` hatch.
    """
    X = _tfidf_like(COUNTS.astype(np.float64))
    df = _run(tmp_path, _write_inputs(tmp_path, X=X), allow_non_count_matrix=True)
    plus = df[df.gene_id == "GENE_PLUS"].iloc[0]
    assert plus.pdui != pytest.approx(0.9, abs=1e-2)
    assert plus.proximal_reads != pytest.approx(10.0)


def test_counts_layer_survives_switch_combine_concat():
    """``ad.concat(join='outer', merge='first')`` must carry layers['counts'].

    ``ema switch combine`` is between clustering (which stashes the layer)
    and ``switch length`` (which reads it); if concat dropped it, the stash
    would never reach the consumer.
    """
    from ema.switch_test.combine import combine_labeled

    def _one(cells):
        a = ad.AnnData(
            X=sp.csr_matrix(_tfidf_like(COUNTS.astype(np.float64))),
            obs=pd.DataFrame(index=cells), var=VAR.copy(),
        )
        a.layers["counts"] = sp.csr_matrix(COUNTS.astype(np.float64))
        return a

    combined = combine_labeled(
        [_one(["a0", "a1", "a2", "a3"]), _one(["b0", "b1", "b2", "b3"])],
        ["A", "B"], group_key="stage",
    )["__all__"]
    assert "counts" in combined.layers
    got = combined.layers["counts"]
    got = got.toarray() if sp.issparse(got) else np.asarray(got)
    np.testing.assert_array_equal(got, np.vstack([COUNTS, COUNTS]))


def test_clustering_stashes_a_counts_layer_before_normalising():
    """The producer half of 0c: raw counts must be stashed before normalize()."""
    from ema.clustering.clustering import _stash_counts_layer

    a = ad.AnnData(
        X=sp.csr_matrix(COUNTS.astype(np.float64)),
        obs=OBS.copy(), var=VAR.copy(),
    )
    _stash_counts_layer(a)
    # Simulate what leiden_tfidf does to .X in place.
    a.X = sp.csr_matrix(_tfidf_like(COUNTS.astype(np.float64)))
    stashed = a.layers["counts"]
    stashed = stashed.toarray() if sp.issparse(stashed) else np.asarray(stashed)
    np.testing.assert_array_equal(stashed, COUNTS)


def test_explicit_pasbed_that_does_not_exist_is_not_silently_replaced(tmp_path):
    """A typo'd --pasbed must abort, not fall back to the walk-up.

    The walk-up would happily find a DIFFERENT run's pasbed.bed sitting next
    to the h5ad and rank against its coordinates.
    """
    h5 = _write_inputs(tmp_path)          # leaves a valid pasbed.bed sibling
    with pytest.raises(FileNotFoundError, match="does not exist"):
        _run(tmp_path, h5, pasbed=str(tmp_path / "typo.bed"))
