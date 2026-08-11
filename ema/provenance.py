"""Append-only provenance ledgers for PeakATail drop-site accounting.

The pipeline computes, at each of several drop points (atlas snapping,
cell-barcode filtering, PAS-to-gene assignment, preprocessing, matrix
concatenation, marker subsetting, ...), which entities survive and which
are dropped. Historically only survivors were kept downstream and the
dropped set + its reason were discarded. This module records *both*
outcomes into two flat, append-only TSV ledgers:

    <run_dir>/provenance/pas_ledger.tsv
    <run_dir>/provenance/cell_ledger.tsv

Design goals:
    * Zero heavy imports at module import time (no pandas/numpy/scanpy).
      Only the standard library (``csv``, ``pathlib``) is used.
    * Cheap/streaming friendly: callers append one row at a time via
      :meth:`ProvenanceLedger.record_pas` / :meth:`record_cell` /
      :meth:`record_drop`; rows are buffered in plain Python lists of
      tuples (fixed-width, no pandas objects) and are only serialized to
      disk on :meth:`flush` / :meth:`write`.
    * Column order is a frozen cross-repo contract — see
      ``PAS_LEDGER_COLUMNS`` / ``CELL_LEDGER_COLUMNS``.

Memory bound: each buffered row is a tuple of ``len(columns)`` short
strings/numbers (no arrays, no references to the original heavy objects
being tracked). For N total records across both ledgers the resident
memory is O(N) small tuples -- roughly a few hundred bytes per row. For
very large N (hundreds of millions of PAS/cell records) callers that
need a hard memory ceiling should call :meth:`flush` periodically (it is
idempotent and safe to call multiple times); subsequent records continue
to buffer and will be written on the next flush. This module intentionally
does not attempt to batch/stream-write incrementally itself, keeping the
implementation simple; periodic ``flush()`` calls are the documented
escape hatch for bounding memory.
"""

from __future__ import annotations

import csv
import io
from pathlib import Path
from typing import Any

__all__ = [
    "PAS_LEDGER_COLUMNS",
    "CELL_LEDGER_COLUMNS",
    "ProvenanceLedger",
    "InvariantError",
    "surviving_count_from_tsv",
    "check_survivor_invariant",
    "record_pas_drops",
    "record_cell_drops",
    "concat_ledger_tsvs",
    "reconcile_dataset_ledgers",
]


class InvariantError(AssertionError):
    """Raised when a provenance survivor-count invariant is violated."""

#: Frozen column order for pas_ledger.tsv. This is a cross-repo contract --
#: do not reorder, rename, or remove columns without a corresponding
#: contract version bump.
PAS_LEDGER_COLUMNS: list[str] = [
    "orig_pas_key",
    "chrom",
    "start",
    "end",
    "strand",
    "unified_pas_id",
    "snap_distance_bp",
    "gene_id",
    "gene_distance_bp",
    "tier",
    "last_stage",
    "dropped_at",
    "drop_reason",
]

#: Frozen column order for cell_ledger.tsv. See ``PAS_LEDGER_COLUMNS``.
CELL_LEDGER_COLUMNS: list[str] = [
    "barcode",
    "dataset_id",
    "total_reads",
    "n_pas",
    "dropped_at",
    "drop_reason",
    "cluster",
]

_KIND_COLUMNS: dict[str, list[str]] = {
    "pas": PAS_LEDGER_COLUMNS,
    "cell": CELL_LEDGER_COLUMNS,
}

_KIND_FILENAMES: dict[str, str] = {
    "pas": "pas_ledger.tsv",
    "cell": "cell_ledger.tsv",
}


def _sanitize(value: Any) -> str:
    """Convert a single ledger field value to a safe TSV cell.

    ``None`` becomes the empty string. Any literal tab/carriage-return/
    newline characters are replaced with a single space so the resulting
    value can never corrupt the tab-separated row structure.
    """
    if value is None:
        return ""
    text = str(value)
    if not text:
        return ""
    for bad in ("\t", "\r", "\n"):
        text = text.replace(bad, " ")
    return text


class ProvenanceLedger:
    """Buffered, append-only recorder for PAS/cell drop provenance.

    Ledgers are written under ``<run_dir>/provenance/``:
        * ``pas_ledger.tsv`` -- one row per PAS/peak entity observed at any
          drop site, columns per :data:`PAS_LEDGER_COLUMNS`.
        * ``cell_ledger.tsv`` -- one row per cell/barcode entity observed
          at any drop site, columns per :data:`CELL_LEDGER_COLUMNS`.

    A surviving entity has ``dropped_at == ""`` and ``drop_reason == ""``.
    A dropped entity has ``dropped_at`` set to the stage id where it died
    (e.g. ``"atlas_snap"``, ``"cb_filter"``, ``"pas_gene"``,
    ``"preprocess"``, ``"matrix_concat"``, ``"marker_subset"``) and a
    human-readable ``drop_reason``.

    Rows are buffered in memory as plain tuples and only serialized to
    disk in :meth:`flush` (alias :meth:`write`), which is idempotent and
    safe to call repeatedly (e.g. once at the end of a run, or
    periodically to bound memory -- see module docstring).
    """

    def __init__(self, run_dir: str | Path) -> None:
        self.run_dir = Path(run_dir)
        self.provenance_dir = self.run_dir / "provenance"
        self._rows: dict[str, list[tuple[str, ...]]] = {"pas": [], "cell": []}

    # -- recording -----------------------------------------------------

    def record_pas(
        self,
        *,
        orig_pas_key: Any,
        chrom: Any = None,
        start: Any = None,
        end: Any = None,
        strand: Any = None,
        unified_pas_id: Any = None,
        snap_distance_bp: Any = None,
        gene_id: Any = None,
        gene_distance_bp: Any = None,
        tier: Any = None,
        last_stage: Any = None,
        dropped_at: Any = "",
        drop_reason: Any = "",
    ) -> None:
        """Append one buffered row to the PAS ledger.

        ``orig_pas_key`` is required; all other fields are optional and
        default to an unknown/empty value. See module docstring for the
        surviving-vs-dropped semantics of ``dropped_at``/``drop_reason``.
        """
        row = (
            _sanitize(orig_pas_key),
            _sanitize(chrom),
            _sanitize(start),
            _sanitize(end),
            _sanitize(strand),
            _sanitize(unified_pas_id),
            _sanitize(snap_distance_bp),
            _sanitize(gene_id),
            _sanitize(gene_distance_bp),
            _sanitize(tier),
            _sanitize(last_stage),
            _sanitize(dropped_at),
            _sanitize(drop_reason),
        )
        self._rows["pas"].append(row)

    def record_cell(
        self,
        *,
        barcode: Any,
        dataset_id: Any,
        total_reads: Any = None,
        n_pas: Any = None,
        dropped_at: Any = "",
        drop_reason: Any = "",
        cluster: Any = None,
    ) -> None:
        """Append one buffered row to the cell ledger.

        ``barcode`` and ``dataset_id`` are required; all other fields are
        optional. See module docstring for surviving-vs-dropped semantics.
        """
        row = (
            _sanitize(barcode),
            _sanitize(dataset_id),
            _sanitize(total_reads),
            _sanitize(n_pas),
            _sanitize(dropped_at),
            _sanitize(drop_reason),
            _sanitize(cluster),
        )
        self._rows["cell"].append(row)

    def record_drop(self, kind: str, *, dropped_at: Any, drop_reason: Any, **fields: Any) -> None:
        """Convenience wrapper for recording a dropped entity.

        Args:
            kind: ``"pas"`` or ``"cell"``.
            dropped_at: required stage id where the entity was dropped.
            drop_reason: required human-readable reason.
            **fields: the remaining keyword args forwarded to
                :meth:`record_pas` / :meth:`record_cell` (e.g.
                ``orig_pas_key=`` for ``kind="pas"`` or ``barcode=``/
                ``dataset_id=`` for ``kind="cell"``).

        Raises:
            ValueError: if ``kind`` is not one of ``{"pas", "cell"}``, or
                if ``dropped_at``/``drop_reason`` are empty (a drop record
                must always carry a stage and a reason).
        """
        if not dropped_at:
            raise ValueError("record_drop requires a non-empty dropped_at")
        if not drop_reason:
            raise ValueError("record_drop requires a non-empty drop_reason")
        if kind == "pas":
            self.record_pas(dropped_at=dropped_at, drop_reason=drop_reason, **fields)
        elif kind == "cell":
            self.record_cell(dropped_at=dropped_at, drop_reason=drop_reason, **fields)
        else:
            raise ValueError(f"record_drop: unknown kind {kind!r}, expected 'pas' or 'cell'")

    # -- counting --------------------------------------------------------

    def count_surviving(self, kind: str) -> int:
        """Return the number of buffered rows of ``kind`` with an empty
        ``dropped_at`` (i.e. entities that survived to their last stage).
        """
        if kind not in _KIND_COLUMNS:
            raise ValueError(f"count_surviving: unknown kind {kind!r}, expected 'pas' or 'cell'")
        dropped_at_index = _KIND_COLUMNS[kind].index("dropped_at")
        return sum(1 for row in self._rows[kind] if row[dropped_at_index] == "")

    @property
    def surviving_pas_count(self) -> int:
        """Number of buffered PAS rows currently marked as surviving."""
        return self.count_surviving("pas")

    # -- persistence -------------------------------------------------------

    def flush(self) -> None:
        """Write both ledgers (header + all buffered rows) to disk.

        Idempotent: re-invoking simply rewrites the current buffer
        contents in full. Creates ``<run_dir>/provenance/`` if needed.
        """
        self.provenance_dir.mkdir(parents=True, exist_ok=True)
        for kind, columns in _KIND_COLUMNS.items():
            path = self.provenance_dir / _KIND_FILENAMES[kind]
            buf = io.StringIO()
            writer = csv.writer(buf, delimiter="\t", lineterminator="\n")
            writer.writerow(columns)
            writer.writerows(self._rows[kind])
            path.write_text(buf.getvalue())

    # alias
    write = flush


def record_pas_drops(
    ledger: "ProvenanceLedger",
    *,
    input_pas_ids,
    annotated_pas_ids,
    final_pas_ids,
    dataset_id: str = "",
    gene_of=None,
    reasons=None,
) -> int:
    """Record one dataset's PAS provenance across the two per-dataset drop sites.

    The per-dataset downstream chain drops PAS at exactly two points after the
    unified/atlas stage:

        input  --annotate(pas_gene)-->  annotated  --preprocess(min_cells)-->  final

    where ``input``  = PAS entering annotate (post cb-filter matrix build),
          ``annotated`` = PAS that got a gene (annotate keeps only the
                          gene-intersection), and
          ``final`` = PAS surviving preprocessing = ``clusters.h5ad`` var_names.

    Records a drop row for every PAS lost at pas_gene (input − annotated) and at
    preprocess (annotated − final), and a surviving row (``dropped_at == ""``,
    ``last_stage="clustering"``) for every final PAS. Returns the surviving
    count, which by construction equals ``len(final_pas_ids)`` — the invariant
    ``surviving == n_vars(clusters.h5ad)``.

    ``gene_of`` optionally maps a surviving pas_id → gene_id; ``reasons`` may
    override the default drop reasons ``{"pas_gene": ..., "preprocess": ...}``.
    """
    gene_of = gene_of or {}
    default_reasons = {
        "pas_gene": "no gene within max_distance (annotate)",
        "preprocess": "below min_cells (preprocessing)",
    }
    if reasons:
        default_reasons.update(reasons)

    annotated = set(annotated_pas_ids)
    final = set(final_pas_ids)
    for p in input_pas_ids:
        if p not in annotated:
            ledger.record_drop(
                "pas", dropped_at="pas_gene", drop_reason=default_reasons["pas_gene"],
                orig_pas_key=p, last_stage="annotate",
            )
    for p in annotated_pas_ids:
        if p not in final:
            ledger.record_drop(
                "pas", dropped_at="preprocess", drop_reason=default_reasons["preprocess"],
                orig_pas_key=p, last_stage="preprocess",
            )
    for p in final_pas_ids:
        # For a surviving PAS the final id IS the unified id (var_name), so set
        # unified_pas_id explicitly — the join key to pasbed.bed / clusters.h5ad
        # (coordinates are sourced there by design, not duplicated onto the
        # ledger). Without this the column was silently NaN for every row.
        ledger.record_pas(
            orig_pas_key=p, unified_pas_id=p, gene_id=gene_of.get(p),
            last_stage="clustering", dropped_at="", drop_reason="",
        )
    return ledger.count_surviving("pas")


def record_cell_drops(
    ledger: "ProvenanceLedger",
    *,
    input_cbs,
    kept_cbs,
    final_cbs=None,
    dataset_id: str = "",
    reasons=None,
) -> int:
    """Record one dataset's cell provenance across cb-filter (+ preprocess).

        input  --cb_filter(min_read)-->  kept  --preprocess(min_genes)-->  final

    Drops at cb_filter (input − kept) and, if ``final_cbs`` is given, at
    preprocess (kept − final). Surviving rows are ``final_cbs`` (or ``kept_cbs``
    when ``final_cbs`` is None). Returns the surviving cell count.
    """
    default_reasons = {
        "cb_filter": "below min_read (cb-filter)",
        "preprocess": "below min_genes (preprocessing)",
    }
    if reasons:
        default_reasons.update(reasons)

    kept = set(kept_cbs)
    for cb in input_cbs:
        if cb not in kept:
            ledger.record_drop(
                "cell", dropped_at="cb_filter", drop_reason=default_reasons["cb_filter"],
                barcode=cb, dataset_id=dataset_id,
            )
    survivors = kept_cbs
    if final_cbs is not None:
        final = set(final_cbs)
        for cb in kept_cbs:
            if cb not in final:
                ledger.record_drop(
                    "cell", dropped_at="preprocess",
                    drop_reason=default_reasons["preprocess"],
                    barcode=cb, dataset_id=dataset_id,
                )
        survivors = final_cbs
    for cb in survivors:
        ledger.record_cell(barcode=cb, dataset_id=dataset_id, dropped_at="", drop_reason="")
    return ledger.count_surviving("cell")


def concat_ledger_tsvs(paths, out_path: str | Path, kind: str = "pas") -> int:
    """Concatenate several same-schema ledger TSVs into one (header written once).

    Missing/empty inputs are skipped. Returns the number of DATA rows written.
    Used by the parent to stitch per-dataset worker sidecars into a cohort
    ledger after the multiprocessing pool joins (sidecar-then-reconcile — no
    shared mutable state across workers).
    """
    if kind not in _KIND_COLUMNS:
        raise ValueError(f"concat_ledger_tsvs: unknown kind {kind!r}")
    columns = _KIND_COLUMNS[kind]
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    n_rows = 0
    buf = io.StringIO()
    writer = csv.writer(buf, delimiter="\t", lineterminator="\n")
    writer.writerow(columns)
    for p in paths:
        p = Path(p)
        if not p.exists():
            continue
        with p.open(newline="") as fh:
            reader = csv.reader(fh, delimiter="\t")
            try:
                next(reader)  # skip that file's header
            except StopIteration:
                continue
            for row in reader:
                if row:
                    writer.writerow(row)
                    n_rows += 1
    out_path.write_text(buf.getvalue())
    return n_rows


def reconcile_dataset_ledgers(datasets, *, cohort_dir=None) -> dict:
    """Reconcile per-dataset worker ledger sidecars after the pool joins.

    Each per-dataset worker writes its own ``<sidecar>/pas_ledger.tsv`` +
    ``cell_ledger.tsv`` (survivor rows = that dataset's final PAS/cells). This
    parent-side step (a) checks the integrity invariant ``surviving PAS ==
    n_vars(clusters.h5ad)`` for EACH dataset, and (b) optionally concatenates
    the per-dataset sidecars into a cohort ledger under ``cohort_dir``.

    IMPORTANT — two-level model: the per-dataset invariant is checked on the
    per-dataset sidecar, NOT on a flat cohort total. A PAS may survive atlas
    snapping (run-level ledger) yet be dropped per-dataset at pas_gene/preprocess,
    so the run-level atlas ledger and the per-dataset ledgers answer different
    questions and are kept SEPARATE (the atlas ledger is not merged in here).

    Args:
        datasets: iterable of ``(ds_id, sidecar_provenance_dir, n_vars)`` where
            ``sidecar_provenance_dir`` contains ``pas_ledger.tsv`` /
            ``cell_ledger.tsv`` and ``n_vars`` is that dataset's
            ``clusters.h5ad`` var count.
        cohort_dir: if given, write concatenated ``pas_ledger.tsv`` +
            ``cell_ledger.tsv`` (per-dataset-scoped survivors) here. Named
            ``cohort_*`` is the CALLER's choice via this dir — this function
            never touches the run-level atlas ledger.

    Returns:
        ``{"datasets": [{ds_id, surviving_pas, n_vars, invariant_ok}...],
           "all_ok": bool, "n_datasets": int}``.
    """
    datasets = list(datasets)
    results = []
    pas_paths, cell_paths = [], []
    for ds_id, sidecar_dir, n_vars in datasets:
        sidecar = Path(sidecar_dir)
        pas_tsv = sidecar / _KIND_FILENAMES["pas"]
        cell_tsv = sidecar / _KIND_FILENAMES["cell"]
        pas_paths.append(pas_tsv)
        cell_paths.append(cell_tsv)
        surviving = surviving_count_from_tsv(pas_tsv, kind="pas") if pas_tsv.exists() else 0
        ok = check_survivor_invariant(surviving, int(n_vars), raise_on_fail=False)
        results.append({
            "ds_id": ds_id,
            "surviving_pas": surviving,
            "n_vars": int(n_vars),
            "invariant_ok": bool(ok),
        })

    if cohort_dir is not None:
        cohort_dir = Path(cohort_dir)
        concat_ledger_tsvs(pas_paths, cohort_dir / _KIND_FILENAMES["pas"], kind="pas")
        concat_ledger_tsvs(cell_paths, cohort_dir / _KIND_FILENAMES["cell"], kind="cell")

    return {
        "datasets": results,
        "all_ok": all(r["invariant_ok"] for r in results),
        "n_datasets": len(results),
    }


def surviving_count_from_tsv(path: str | Path, kind: str = "pas") -> int:
    """Count survivors (``dropped_at == ""``) in an on-disk ledger TSV.

    Reads only the ``dropped_at`` column, so it works for both the pas and
    cell ledgers without loading pandas. A survivor is a data row whose
    ``dropped_at`` field is empty.

    Args:
        path: path to a ``pas_ledger.tsv`` / ``cell_ledger.tsv``.
        kind: ``"pas"`` or ``"cell"`` (selects the expected header).

    Raises:
        ValueError: if ``kind`` is unknown or the header lacks ``dropped_at``.
    """
    if kind not in _KIND_COLUMNS:
        raise ValueError(f"surviving_count_from_tsv: unknown kind {kind!r}")
    p = Path(path)
    with p.open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return 0
        if "dropped_at" not in header:
            raise ValueError(
                f"ledger at {p} has no 'dropped_at' column; header={header!r}"
            )
        idx = header.index("dropped_at")
        return sum(1 for row in reader if len(row) > idx and row[idx] == "")


def check_survivor_invariant(
    surviving: int,
    n_vars: int,
    *,
    label: str = "pas_ledger survivors == n_vars(clusters.h5ad)",
    raise_on_fail: bool = True,
) -> bool:
    """Assert the provenance integrity invariant.

    The contract invariant is::

        rows(pas_ledger where dropped_at == "") == n_vars(clusters.h5ad)

    i.e. every PAS that survives all drop sites is exactly a variable in the
    final clustered AnnData, and vice versa. This function compares the two
    counts. It is the single reusable check used by both the engine (post-run
    self-test) and the hub (contract validation) so the two never drift.

    Args:
        surviving: number of surviving PAS ledger rows (``dropped_at == ""``).
        n_vars: ``clusters.h5ad`` ``n_vars`` (number of PAS in the final matrix).
        label: description used in the error message.
        raise_on_fail: if True (default) raise :class:`InvariantError` on
            mismatch; if False return ``False`` instead.

    Returns:
        ``True`` when ``surviving == n_vars``.

    Raises:
        InvariantError: when the counts differ and ``raise_on_fail`` is True.
    """
    ok = int(surviving) == int(n_vars)
    if not ok and raise_on_fail:
        raise InvariantError(
            f"provenance invariant violated ({label}): "
            f"surviving={surviving} != n_vars={n_vars} "
            f"(delta={int(surviving) - int(n_vars)})"
        )
    return ok
