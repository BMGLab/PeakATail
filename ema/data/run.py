"""E4: the ``ema.data.Run`` loader — a manifest-driven, pandas-based read
facade over one on-disk PeakATail run directory.

This module SUPERSEDES the path-hardcoding in ``scripts/harvest_report_data.py``
and ``scripts/harvest_matrix_per_dataset.py``: both scripts hand-assembled
paths like ``<run>/provenance/pas_ledger.tsv`` or ``<run>/07_clustering/<ds>/
clusters.h5ad`` themselves. :class:`Run` centralizes that path knowledge,
resolving artifacts from ``run_manifest.json`` (see
``ema/outputs.py::OutputManager.write_manifest``) first and falling back to
the conventional on-disk layout (``ema/outputs.py::OutputManager.
_auto_discover_artifacts``) when the manifest doesn't register an artifact.

ENGINE-SELF-CONTAINED: this module must NEVER import ``peakatail_contract``
or any ``peakatail-hub`` package -- the engine repo builds and tests
independently of the hub. It mirrors the METHOD NAMES of the hub's read
facade (``/home/user/D/peakatail-hub/packages/io/src/peakatail_io/run.py``,
a pydantic-typed reader) so callers who know one can guess the other, but
every accessor here returns a plain ``pandas.DataFrame`` / ``dict`` / ``list``
/ ``anndata.AnnData`` -- never a pydantic model.

Lazy + cached: no-argument accessors are ``functools.cached_property`` (each
artifact is read from disk at most once per ``Run`` instance); accessors that
take arguments (``clusters``, ``n_vars``, ``canonical_cluster``) use a manual
per-argument cache dict for the same reason.
"""
from __future__ import annotations

import json
from functools import cached_property
from pathlib import Path
from typing import Any

import pandas as pd

__all__ = ["Run", "RunReadError"]

#: BED6 column names for pasbed.bed / annotatedpas.bed (no header on disk).
_PASBED_COLUMNS = ["chrom", "start", "end", "pas_id", "score", "strand"]

#: Trailing columns ``annotatedpas.bed`` may carry beyond BED6, in the fixed
#: order ``ema/outputs.py::write_pas_gene_artifacts`` appends them. D9 added
#: the last three (atlas-snap / internal-priming become annotate-not-drop by
#: default; see ``ema/datasets/atlas_snap.py`` / ``ema/experimental/
#: internal_priming.py``) -- always present together with gene_id (never
#: gene_id alone anymore going forward), but kept as separate optional
#: columns here so OLD runs' 7-column (BED6 + gene_id) annotatedpas.bed
#: files still parse correctly.
_PASBED_EXTRA_COLUMNS = ["gene_id", "atlas_match", "atlas_distance_bp", "internal_priming"]


def _read_pasbed_like(path: Path) -> pd.DataFrame:
    """Read a BED6(+extra) file with no on-disk header.

    Plain ``pasbed.bed`` is exactly BED6. ``annotatedpas.bed`` extends it
    with up to 4 trailing columns (see :data:`_PASBED_EXTRA_COLUMNS`) --
    older runs may only have ``gene_id`` (7 cols), newer runs have all 4
    (10 cols). Naming columns positionally like this (instead of a fixed
    ``names=_PASBED_COLUMNS`` as before) matters: pandas silently
    misinterprets extra un-named trailing columns as an index level when
    the file has MORE columns than a fixed ``names=`` list, corrupting
    every column's dtype/values -- this reader avoids that by sizing
    ``names`` to the file's actual width.
    """
    raw = pd.read_csv(path, sep="\t", header=None)
    n_cols = raw.shape[1]
    names = list(_PASBED_COLUMNS)
    for extra in _PASBED_EXTRA_COLUMNS:
        if len(names) >= n_cols:
            break
        names.append(extra)
    while len(names) < n_cols:
        names.append(f"extra_{len(names) - len(_PASBED_COLUMNS) + 1}")
    raw.columns = names[:n_cols]
    return raw


class RunReadError(Exception):
    """Raised for any problem reading a run directory or one of its
    artifacts: missing/unparseable ``run_manifest.json``, a missing artifact
    file, an ambiguous ``clusters.h5ad`` lookup, etc. Always carries the
    path(s) that were tried so the caller doesn't have to go spelunking.
    """


class Run:
    """Read-only facade over one on-disk PeakATail pipeline run directory.

    This is the E4 loader: the single place downstream tooling (report
    harvesting, notebooks, the hub's ingest step) should go to read a run's
    artifacts, instead of re-deriving ``<run>/07_clustering/<ds>/
    clusters.h5ad``-style paths by hand. Construct via :meth:`from_dir`; the
    resulting object owns ``root`` (the actual on-disk directory the
    manifest was loaded from -- NOT the ``root`` string recorded *inside*
    the manifest, which may be stale if the run directory was copied/moved)
    and the parsed manifest dict.
    """

    def __init__(self, root: Path, manifest: dict[str, Any]) -> None:
        self._root = root
        self._manifest = manifest
        self._cluster_cache: dict[str | None, Any] = {}
        self._n_vars_cache: dict[str | None, int] = {}

    # -- construction --------------------------------------------------- #

    @classmethod
    def from_dir(cls, path: str | Path) -> "Run":
        """Load ``run_manifest.json`` from ``path``.

        ``path`` may be a run directory (containing ``run_manifest.json``)
        or a direct path to a ``run_manifest.json`` file. Raises
        :class:`RunReadError` (never returns a partially-built ``Run``) if
        the manifest is missing, unreadable, or not valid JSON.
        """
        p = Path(path).expanduser()
        if p.is_file():
            manifest_path = p
            root = p.parent.resolve()
        else:
            root = p.resolve()
            manifest_path = root / "run_manifest.json"

        if not manifest_path.exists():
            raise RunReadError(
                f"no run_manifest.json found at {manifest_path} -- is {root} a "
                "PeakATail run output directory produced with manifest writing "
                "(OutputManager.write_manifest) enabled?"
            )
        try:
            raw = manifest_path.read_text()
        except OSError as exc:
            raise RunReadError(f"could not read {manifest_path}: {exc}") from exc
        try:
            manifest = json.loads(raw)
        except json.JSONDecodeError as exc:
            raise RunReadError(
                f"{manifest_path} exists but is not valid JSON: {exc}"
            ) from exc
        if not isinstance(manifest, dict):
            raise RunReadError(
                f"{manifest_path} did not parse to a JSON object (got {type(manifest).__name__})"
            )
        return cls(root=root, manifest=manifest)

    # -- identity / config ------------------------------------------------ #

    @property
    def manifest(self) -> dict[str, Any]:
        """The parsed ``run_manifest.json`` contents, verbatim."""
        return self._manifest

    @property
    def root(self) -> Path:
        """The actual on-disk run directory (see class docstring)."""
        return self._root

    @property
    def run_id(self) -> str:
        return str(self._manifest.get("run_id") or self._root.name)

    @cached_property
    def config(self) -> dict[str, Any]:
        """The full ``resolved_config`` dict (see
        ``ema/outputs.py::build_resolved_run_config``)."""
        return dict(self._manifest.get("resolved_config") or {})

    @property
    def args(self) -> dict[str, Any]:
        return dict(self.config.get("args") or {})

    @property
    def filters(self) -> dict[str, Any]:
        return dict(self.config.get("filters") or {})

    @property
    def variables(self) -> dict[str, Any]:
        return dict(self.config.get("variables") or {})

    @property
    def directories(self) -> dict[str, Any]:
        return dict(self.config.get("directories") or {})

    @cached_property
    def datasets(self) -> list[str]:
        """``dataset_id`` strings from ``manifest["datasets"]``, in order,
        deduplicated. Tolerant of the dict having ``dataset_id``/``id``/
        ``name`` (the couple of shapes ``OutputManager._datasets_from_config``
        and hand-authored manifests use) or being a bare string already."""
        out: list[str] = []
        for d in self._manifest.get("datasets") or []:
            if isinstance(d, dict):
                ds_id = d.get("dataset_id") or d.get("id") or d.get("name")
            else:
                ds_id = d
            if ds_id and str(ds_id) not in out:
                out.append(str(ds_id))
        return out

    # -- path resolution --------------------------------------------------- #

    def _resolve(self, rel_or_abs: str | Path) -> Path:
        p = Path(rel_or_abs)
        return p if p.is_absolute() else self._root / p

    def _resolve_artifact(self, schema_name: str, *conventional_relpaths: str) -> Path:
        """Resolve one artifact's path: prefer the ``run_manifest.json``
        entry whose ``schema_name`` matches (registered explicitly or by
        ``OutputManager._auto_discover_artifacts``); fall back to each
        ``conventional_relpaths`` under :attr:`root` in order when the
        manifest omits the entry, or when the registered path doesn't
        actually exist on disk (a stale/incomplete manifest is not fatal by
        itself). Raises :class:`RunReadError` naming every path tried if
        nothing resolves to an existing file.

        Multiple conventional paths matter because the SAME logical artifact
        lands in different places by run shape: e.g. the provenance ledgers
        sit at ``provenance/pas_ledger.tsv`` when a run-level (atlas-snap)
        drop ledger exists, but on a multi-sample no-atlas run only the
        reconciled per-dataset ledger at ``provenance/by_dataset/
        pas_ledger.tsv`` is written (verified on a real Laughney run).
        """
        tried: list[Path] = []
        for artifact in self._manifest.get("artifacts") or []:
            if artifact.get("schema_name") == schema_name:
                candidate = self._resolve(artifact.get("path", ""))
                tried.append(candidate)
                if candidate.exists():
                    return candidate
                break  # one manifest entry per schema_name is expected

        for relpath in conventional_relpaths:
            conventional = self._resolve(relpath)
            tried.append(conventional)
            if conventional.exists():
                return conventional

        raise RunReadError(
            f"could not resolve artifact schema_name={schema_name!r} for "
            f"run_id={self.run_id!r}: tried {[str(t) for t in tried]}"
        )

    def _resolve_long_table(self, schema_name: str, base_name: str) -> Path:
        """Like :meth:`_resolve_artifact`, but for the E5 long tables
        (``findings_long`` / ``length_long``) written by
        ``ema/switch_test/long_output.py::write_long_table``, which writes
        ``<base_name>.parquet`` when a parquet engine (pyarrow/fastparquet)
        is installed and falls back to a sibling ``<base_name>.tsv``
        otherwise -- so this MUST try both extensions, not just one.
        """
        for artifact in self._manifest.get("artifacts") or []:
            if artifact.get("schema_name") == schema_name:
                candidate = self._resolve(artifact.get("path", ""))
                if candidate.exists():
                    return candidate

        parquet_path = self._resolve(f"{base_name}.parquet")
        if parquet_path.exists():
            return parquet_path
        tsv_path = self._resolve(f"{base_name}.tsv")
        if tsv_path.exists():
            return tsv_path

        raise RunReadError(
            f"could not resolve {schema_name!r} long table for run_id={self.run_id!r}: "
            f"tried manifest artifacts, {parquet_path}, {tsv_path}"
        )

    # -- provenance ledgers ------------------------------------------------ #

    @cached_property
    def pas_ledger(self) -> pd.DataFrame:
        """``provenance/pas_ledger.tsv`` -- columns per
        ``ema/provenance.py::PAS_LEDGER_COLUMNS``."""
        path = self._resolve_artifact(
            "PasLedgerRow", "provenance/pas_ledger.tsv", "provenance/by_dataset/pas_ledger.tsv"
        )
        return pd.read_csv(path, sep="\t")

    @cached_property
    def cell_ledger(self) -> pd.DataFrame:
        """``provenance/cell_ledger.tsv`` -- columns per
        ``ema/provenance.py::CELL_LEDGER_COLUMNS``."""
        path = self._resolve_artifact(
            "CellLedgerRow", "provenance/cell_ledger.tsv", "provenance/by_dataset/cell_ledger.tsv"
        )
        return pd.read_csv(path, sep="\t")

    @cached_property
    def provenance(self) -> "_Provenance":
        """A small namespace grouping every provenance artifact:
        ``.pas`` / ``.cell`` (same as :attr:`pas_ledger` / :attr:`cell_ledger`),
        ``.reconcile_summary`` (parsed ``provenance/reconcile_summary.json``,
        or ``None`` if absent -- written by ``ema.provenance.
        reconcile_dataset_ledgers`` via ``ema/main.py``), and ``.by_dataset``.

        ``.by_dataset`` shape note: the E3 sidecar-then-reconcile step
        (``reconcile_dataset_ledgers(..., cohort_dir=<run>/provenance/
        by_dataset)``) concatenates every dataset's per-dataset ledger
        *sidecar* into ONE flat ``provenance/by_dataset/pas_ledger.tsv`` +
        ``cell_ledger.tsv`` pair -- there are no per-dataset subdirectories
        in the real engine output today. ``cell_ledger.tsv`` carries a
        ``dataset_id`` column (``CELL_LEDGER_COLUMNS``), so this CAN split
        it into a ``{dataset_id: DataFrame}`` dict; ``pas_ledger.tsv`` has
        NO ``dataset_id`` column (``PAS_LEDGER_COLUMNS``), so it cannot be
        split the same way and is intentionally left out of ``.by_dataset``
        (use the flat, unsplit :attr:`pas_ledger` for that). If a future
        layout instead writes per-dataset subdirectories (``by_dataset/<ds>/
        pas_ledger.tsv``), that shape is read directly, one DataFrame per
        subdirectory. Returns ``{}`` if neither shape is found.
        """
        reconcile_path = self._resolve("provenance/reconcile_summary.json")
        reconcile_summary: dict[str, Any] | None = None
        if reconcile_path.exists():
            try:
                reconcile_summary = json.loads(reconcile_path.read_text())
            except (OSError, json.JSONDecodeError):
                reconcile_summary = None

        by_dataset: dict[str, pd.DataFrame] = {}
        by_dataset_dir = self._resolve("provenance/by_dataset")
        if by_dataset_dir.is_dir():
            subdirs = [d for d in sorted(by_dataset_dir.iterdir()) if d.is_dir()]
            if subdirs:
                # Hypothetical/future per-dataset-subdirectory layout.
                for d in subdirs:
                    pas_p = d / "pas_ledger.tsv"
                    cell_p = d / "cell_ledger.tsv"
                    if pas_p.exists():
                        by_dataset[d.name] = pd.read_csv(pas_p, sep="\t")
                    elif cell_p.exists():
                        by_dataset[d.name] = pd.read_csv(cell_p, sep="\t")
            else:
                # Real (current) engine layout: flat concatenated files.
                cell_p = by_dataset_dir / "cell_ledger.tsv"
                if cell_p.exists():
                    cell_df = pd.read_csv(cell_p, sep="\t")
                    if "dataset_id" in cell_df.columns:
                        for ds_id, sub in cell_df.groupby("dataset_id"):
                            by_dataset[str(ds_id)] = sub.reset_index(drop=True)

        return _Provenance(
            pas=self.pas_ledger,
            cell=self.cell_ledger,
            reconcile_summary=reconcile_summary,
            by_dataset=by_dataset,
        )

    # -- coordinates --------------------------------------------------------- #

    @cached_property
    def pasbed(self) -> pd.DataFrame:
        """``pasbed.bed`` (or ``annotatedpas.bed`` fallback), BED6, no
        header, as a DataFrame with columns
        ``[chrom, start, end, pas_id, score, strand]``."""
        path = self._resolve_artifact("pasbed.bed", "pasbed.bed")
        if not path.exists():
            path = self._resolve_artifact("annotatedpas.bed", "annotatedpas.bed")
        return _read_pasbed_like(path)

    # -- E5 normalized long tables -------------------------------------------- #

    @cached_property
    def findings(self) -> pd.DataFrame:
        """``findings_long.parquet`` (or the ``.tsv`` fallback) -- see
        ``ema/switch_test/long_output.py::write_long_table``."""
        path = self._resolve_long_table("FindingRow", "findings_long")
        return pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path, sep="\t")

    @cached_property
    def length(self) -> pd.DataFrame:
        """``length_long.parquet`` (or the ``.tsv`` fallback)."""
        path = self._resolve_long_table("LengthRow", "length_long")
        return pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path, sep="\t")

    # -- genes --------------------------------------------------------------- #

    @cached_property
    def genes(self) -> list[str]:
        """Sorted, deduplicated ``gene_id`` values, empties/NaN dropped.

        NOTE: ``pasbed.bed`` is plain BED6 (see :attr:`pasbed`) and carries
        no ``gene_id`` column, so despite the name this is derived from
        :attr:`pas_ledger`'s ``gene_id`` column (falling back to
        :attr:`findings`'s ``gene_id`` column if the ledger is unavailable)
        -- either is what ``scripts/harvest_report_data.py`` /
        ``scripts/harvest_matrix_per_dataset.py`` actually read gene ids
        from.
        """
        ids: set[str] = set()
        for accessor in (lambda: self.pas_ledger, lambda: self.findings):
            try:
                df = accessor()
            except RunReadError:
                continue
            if "gene_id" in df.columns:
                ids.update(str(g) for g in df["gene_id"].dropna().tolist())
            if ids:
                break
        ids.discard("")
        ids.discard("nan")
        return sorted(ids)

    # -- clusters.h5ad --------------------------------------------------------- #

    def _cluster_candidates(self) -> list[Path]:
        candidates = [
            self._resolve(a["path"])
            for a in (self._manifest.get("artifacts") or [])
            if a.get("format") == "h5ad"
            and ("cluster" in str(a.get("stage", "")).lower() or a.get("schema_name") == "clusters.h5ad")
        ]
        if not candidates:
            clustering_dir = self._resolve("07_clustering")
            if clustering_dir.is_dir():
                candidates = sorted(clustering_dir.glob("*/clusters.h5ad"))
        return candidates

    def _resolve_cluster_path(self, dataset_id: str | None = None) -> Path:
        candidates = self._cluster_candidates()

        if dataset_id is not None:
            scoped = [p for p in candidates if p.parent.name == dataset_id]
            if scoped:
                if len(scoped) > 1:
                    raise RunReadError(
                        f"ambiguous clusters.h5ad for dataset_id={dataset_id!r}: "
                        f"{[str(p) for p in scoped]!r}"
                    )
                return scoped[0]
            conventional = self._resolve(f"07_clustering/{dataset_id}/clusters.h5ad")
            if conventional.exists():
                return conventional
            raise RunReadError(
                f"no clusters.h5ad found for dataset_id={dataset_id!r} in run_id="
                f"{self.run_id!r}: tried manifest artifacts and {conventional}"
            )

        if len(candidates) == 1:
            return candidates[0]
        if not candidates:
            raise RunReadError(
                f"no clusters.h5ad artifacts registered in run_manifest.json and none "
                f"found under {self._resolve('07_clustering')}/*/clusters.h5ad "
                f"(run_id={self.run_id!r})"
            )
        raise RunReadError(
            f"ambiguous clusters.h5ad for run_id={self.run_id!r}: {len(candidates)} "
            f"candidates {[str(p) for p in candidates]!r}; pass dataset_id to disambiguate"
        )

    def clusters(self, dataset_id: str | None = None):
        """Open the ``clusters.h5ad`` for ``dataset_id`` (or the run's
        sole/unambiguous one) with ``backed='r'`` -- the expression matrix
        stays on disk (h5py-backed access) instead of being fully loaded
        into memory. ``import anndata`` happens here, not at module import
        time, so importing :mod:`ema.data.run` stays cheap for callers that
        never touch a cluster h5ad.

        Raises :class:`RunReadError` if ``dataset_id`` is omitted and more
        than one ``clusters.h5ad`` exists for this run (ambiguous), or if
        no matching file can be found at all.
        """
        if dataset_id in self._cluster_cache:
            return self._cluster_cache[dataset_id]
        import anndata as ad

        path = self._resolve_cluster_path(dataset_id)
        if not path.exists():
            raise RunReadError(f"clusters.h5ad not found at {path} (dataset_id={dataset_id!r})")
        adata = ad.read_h5ad(path, backed="r")
        self._cluster_cache[dataset_id] = adata
        return adata

    def n_vars(self, dataset_id: str | None = None) -> int:
        """``adata.n_vars`` for :meth:`clusters` (``dataset_id``), without
        materializing the full expression matrix."""
        if dataset_id in self._n_vars_cache:
            return self._n_vars_cache[dataset_id]
        n = int(self.clusters(dataset_id).n_vars)
        self._n_vars_cache[dataset_id] = n
        return n

    # -- canonical cluster (B6) ------------------------------------------------ #

    def canonical_cluster(self, dataset_id: str, leiden: str | int) -> str | None:
        """Look up the canonical, cross-sample-comparable label for
        ``(dataset_id, leiden)``.

        Reads ``obs['canonical_cluster']`` off ``clusters(dataset_id)``,
        matched against ``obs['leiden'] == str(leiden)`` -- bug B6 (the
        canonical map was written but never round-tripped into the h5ad's
        ``obs``) is fixed engine-side by
        ``ema/clustering/cross_dataset/roundtrip.py::write_canonical_clusters``,
        which is called from ``ema/main.py`` right after cross-dataset
        cluster matching, so ``obs['canonical_cluster']`` is the current
        source of truth (not a separate ``cross_dataset/
        canonical_cluster_map.tsv`` re-parse). Returns ``None`` -- never
        raises -- when the h5ad has no ``canonical_cluster``/``leiden``
        obs columns (single-dataset run, matching step didn't run, or the
        round-trip was skipped) or when ``leiden`` isn't a known label for
        this dataset.
        """
        adata = self.clusters(dataset_id)
        obs = adata.obs
        if "canonical_cluster" not in obs.columns or "leiden" not in obs.columns:
            return None
        mask = obs["leiden"].astype(str) == str(leiden)
        if not mask.any():
            return None
        values = obs.loc[mask, "canonical_cluster"].dropna()
        if values.empty:
            return None
        return str(values.iloc[0])


class _Provenance:
    """Namespace returned by :attr:`Run.provenance` -- see its docstring."""

    __slots__ = ("pas", "cell", "reconcile_summary", "by_dataset")

    def __init__(
        self,
        pas: pd.DataFrame,
        cell: pd.DataFrame,
        reconcile_summary: dict[str, Any] | None,
        by_dataset: dict[str, pd.DataFrame],
    ) -> None:
        self.pas = pas
        self.cell = cell
        self.reconcile_summary = reconcile_summary
        self.by_dataset = by_dataset
