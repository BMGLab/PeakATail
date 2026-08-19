"""Regression: the count matrix must stay keyed to its PAS IDs.

``make_dataframe()`` returns ``(sparse_matrix, pas_ids, collist)``.  ``mmread``
gives a FULL-height matrix (rows ``1..header_max``), but ``pas_ids`` was built
as ``np.sort(np.unique(coo.row + 1))`` — only the rows that actually carry
counts.  Any PAS whose reads all belong to barcodes dropped by the ``min_read``
CB filter leaves an EMPTY row behind, so the two drift apart.

``ema.annotate.annotate()`` then does::

    pas_id_to_row = {pid: i for i, pid in enumerate(pas_ids)}   # POSITION
    filtered_sparse = sparse_matrix.tocsr()[keep_rows, :]

i.e. it indexes matrix rows by POSITION IN ``pas_ids`` while the matrix is
indexed by PAS ID.  Every PAS below a gap therefore ships another PAS's counts.

These tests pin the invariant that makes the positional lookup legal:
``sparse_matrix.shape[0] == len(pas_ids)`` and row ``i`` holds PAS
``pas_ids[i]``.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.io as sci
import scipy.sparse as sp

from ema.matrixfilter import make_dataframe


# PAS 2 is empty (its only barcode fell below min_read) -> a gap at row 2.
# PAS 5 is empty too.  Rows that DO carry counts: 1, 3, 4, 6.
_MTX = """\
%%MatrixMarket matrix coordinate integer general
6 2 6
1 1 10
3 1 30
3 2 33
4 2 44
6 1 60
6 2 66
"""

_EXPECTED = {
    1: [10, 0],
    3: [30, 33],
    4: [0, 44],
    6: [60, 66],
}


@pytest.fixture()
def filtered_mtx(tmp_path: Path) -> Path:
    p = tmp_path / "filterdmatrix.mtx"
    p.write_text(_MTX)
    return p


def test_make_dataframe_matrix_is_row_aligned_with_pas_ids(filtered_mtx: Path) -> None:
    """Row i of the returned matrix must be PAS ``pas_ids[i]``."""
    matrix, pas_ids, _ = make_dataframe(matrixpath=str(filtered_mtx), collist=["CB1", "CB2"])

    assert list(pas_ids) == [1, 3, 4, 6]
    assert matrix.shape[0] == len(pas_ids), (
        f"matrix has {matrix.shape[0]} rows but pas_ids has {len(pas_ids)} entries — "
        "positional row lookup in annotate() is then off by the number of empty rows"
    )
    dense = np.asarray(matrix.todense())
    for i, pid in enumerate(pas_ids):
        assert list(dense[i]) == _EXPECTED[int(pid)], (
            f"row {i} should hold PAS {pid} counts {_EXPECTED[int(pid)]}, got {list(dense[i])}"
        )


def test_annotate_ships_each_pas_its_own_counts(filtered_mtx: Path, tmp_path: Path) -> None:
    """End-to-end: annotated_matrix.mtx row for PAS p must hold PAS p's counts."""
    import ema.config as cfg
    from ema.annotate.annotate import annotate

    run = tmp_path / "run_TS"
    (run / "05_annotated_matrix").mkdir(parents=True)
    cfg.set_directory_config(output_dir=run)

    pasbed = Path(cfg.directory_config.pasbed)
    pasbed.parent.mkdir(parents=True, exist_ok=True)
    pasbed.write_text(
        "".join(f"chr1\t{p*100}\t{p*100+50}\t{p}\t0\t+\n" for p in (1, 3, 4, 6))
    )

    matrix, pas_ids, collist = make_dataframe(
        matrixpath=str(filtered_mtx), collist=["CB1", "CB2"]
    )
    genes = pd.DataFrame({"gene_id": [f"ENSG{p}" for p in (1, 3, 4, 6)]}, index=[1, 3, 4, 6])

    result = annotate(sparse_matrix=matrix, pas_ids=pas_ids, collist=collist, genes=genes)

    dense = np.asarray(sp.csr_matrix(result.sparse_matrix).todense())
    for i, pid in enumerate(result.pas_ids):
        assert list(dense[i]) == _EXPECTED[int(pid)], (
            f"annotate() gave PAS {pid} the counts {list(dense[i])}, "
            f"expected {_EXPECTED[int(pid)]}"
        )

    # ...and the same must hold for the file actually written to disk, which is
    # what every downstream consumer reads.
    on_disk = np.asarray(
        sci.mmread(str(cfg.directory_config.annotated_matrix)).todense()
    )
    for i, pid in enumerate(result.pas_ids):
        assert list(on_disk[i]) == _EXPECTED[int(pid)]


def test_annotate_rejects_a_full_height_matrix(filtered_mtx: Path, tmp_path: Path) -> None:
    """The pre-fix shape (full-height matrix + compacted pas_ids) must RAISE.

    This is the exact pairing ``make_dataframe()`` returned between 04e0b3a and
    this fix.  Silently accepting it is what mis-keyed 99.87% of annotated PAS
    on the testis mouse1 run; the guard turns it into a loud failure for any
    caller that builds the two halves itself.
    """
    import ema.config as cfg
    from ema.annotate.annotate import annotate

    run = tmp_path / "run_guard"
    (run / "05_annotated_matrix").mkdir(parents=True)
    cfg.set_directory_config(output_dir=run)

    full_height = sp.csc_matrix(sci.mmread(str(filtered_mtx)))
    _, pas_ids, _ = make_dataframe(matrixpath=str(filtered_mtx), collist=["CB1", "CB2"])
    assert full_height.shape[0] == 6 and len(pas_ids) == 4  # the drift, reproduced

    genes = pd.DataFrame({"gene_id": [f"ENSG{p}" for p in (1, 3, 4, 6)]}, index=[1, 3, 4, 6])
    with pytest.raises(ValueError, match="rows are looked up by position"):
        annotate(sparse_matrix=full_height, pas_ids=pas_ids,
                 collist=["CB1", "CB2"], genes=genes)


def test_preprocessing_rejects_row_count_mismatch() -> None:
    """``pas_ids`` labels rows positionally — a height mismatch must raise."""
    from ema.matrixfilter import preprocessing

    with pytest.raises(ValueError, match="pas_ids labels the rows positionally"):
        preprocessing(
            sparse_matrix=sp.csr_matrix(np.ones((6, 2), dtype=int)),
            pas_ids=np.array([1, 3, 4, 6]),
            collist=["CB1", "CB2"],
            min_cells=0,
            min_genes=0,
        )


def test_write_annotated_matrix_rejects_row_count_mismatch(tmp_path: Path) -> None:
    """annotated_pas_ids.tsv is the row index of annotated_matrix.mtx."""
    import ema.config as cfg
    from ema.outputs import write_annotated_matrix

    run = tmp_path / "run_persist"
    (run / "05_annotated_matrix").mkdir(parents=True)
    cfg.set_directory_config(output_dir=run)

    with pytest.raises(ValueError, match="row index of annotated_matrix.mtx"):
        write_annotated_matrix(
            run, "ds1",
            sp.csr_matrix(np.ones((6, 2), dtype=int)),
            np.array([1, 3, 4, 6]),
            ["CB1", "CB2"],
        )


def test_no_empty_rows_is_unchanged(tmp_path: Path) -> None:
    """The gap-free case (every row non-empty) must behave exactly as before.

    Guards against the re-keying subset perturbing runs that never tripped the
    bug — pas_ids == 1..N there, so the subset is the identity.
    """
    p = tmp_path / "dense.mtx"
    p.write_text(
        "%%MatrixMarket matrix coordinate integer general\n"
        "3 2 4\n1 1 11\n2 2 22\n3 1 31\n3 2 32\n"
    )
    matrix, pas_ids, _ = make_dataframe(matrixpath=str(p), collist=["CB1", "CB2"])
    assert list(pas_ids) == [1, 2, 3]
    assert matrix.shape == (3, 2)
    assert np.asarray(matrix.todense()).tolist() == [[11, 0], [0, 22], [31, 32]]


def test_write_pas_gene_artifacts_rejects_length_mismatch(tmp_path: Path) -> None:
    """pas_ids/gene_ids are a positional pair — a length mismatch must raise."""
    import ema.config as cfg
    from ema.outputs import write_pas_gene_artifacts

    run = tmp_path / "run_pasgene"
    run.mkdir(parents=True)
    cfg.set_directory_config(output_dir=run)

    with pytest.raises(ValueError, match="index the same annotated rows"):
        write_pas_gene_artifacts(
            run, "ds1", np.array([1, 3, 4, 6]), np.array(["A", "B", "C"]),
        )


def test_write_pas_gene_artifacts_does_not_realign_a_series(tmp_path: Path) -> None:
    """A gene_id Series must be paired POSITIONALLY, not by its own index.

    ``annotate()`` builds ``gene_ids`` as ``genes.loc[keep_pas_ids, 'gene_id']``
    — a Series indexed by PAS ID.  Handing that Series straight to the
    DataFrame constructor would re-align it against a fresh RangeIndex and
    write NaN for every PAS whose ID is not also a 0-based position.
    """
    import ema.config as cfg
    from ema.outputs import write_pas_gene_artifacts

    run = tmp_path / "run_series"
    run.mkdir(parents=True)
    cfg.set_directory_config(output_dir=run)

    pas_ids = np.array([1, 3, 4, 6])
    gene_series = pd.Series([f"ENSG{p}" for p in pas_ids], index=pas_ids)
    tsv, _ = write_pas_gene_artifacts(run, "ds1", pas_ids, gene_series)

    got = pd.read_csv(tsv, sep="\t")
    assert list(got["pas_id"]) == [1, 3, 4, 6]
    assert list(got["gene_id"]) == ["ENSG1", "ENSG3", "ENSG4", "ENSG6"]
