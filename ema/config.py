from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path


# ---------------------------------------------------------------------------
# Lazy args proxy — replaces the old ``args = cli()`` call at import time.
#
# ``variable_config`` and ``filter_config`` below still reference ``args.*``,
# so we cannot remove the name entirely.  Making it lazy means the CLI parse
# only runs when the first attribute is actually *read*, not at import time.
# This removes the import-time side-effect (argparse invocation, potential
# circular-import risk) while preserving the legacy contract.
# ---------------------------------------------------------------------------

class _LazyArgs:
    """Proxy that defers ``cli()`` until the first attribute access."""

    __slots__ = ("_ns",)

    def __init__(self) -> None:
        object.__setattr__(self, "_ns", None)

    def _get(self):
        ns = object.__getattribute__(self, "_ns")
        if ns is None:
            from ema.cli import cli
            ns = cli()
            object.__setattr__(self, "_ns", ns)
        return ns

    def __getattr__(self, name: str):
        return getattr(self._get(), name)

    def __setattr__(self, name: str, value) -> None:
        if name == "_ns":
            object.__setattr__(self, name, value)
        else:
            setattr(self._get(), name, value)

    def __repr__(self) -> str:  # pragma: no cover
        ns = object.__getattribute__(self, "_ns")
        if ns is None:
            return "<_LazyArgs (not yet initialised)>"
        return f"<_LazyArgs {vars(ns)!r}>"


args = _LazyArgs()


# ---------------------------------------------------------------------------
# DirectoryConfig — frozen dataclass.  A single instance (``directory_config``)
# is exported as the module-level singleton.  All path properties derive from
# ``output_dir`` so changing the run root via ``set_directory_config`` (or
# directly via ``object.__setattr__``) cascades to every path automatically.
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class DirectoryConfig:
    """Single source of truth for every file/dir the pipeline reads or writes.

    Constructed once per run after the timestamped run directory is created
    (see ema/main.py).  All path properties derive from ``output_dir`` so
    changing the run root automatically updates every derived path.
    """

    # ─── injected at construction time ───────────────────────────────────
    output_dir: Path = Path("emaout")
    bam_dir: str | None = None
    gtf_dir: str | None = None
    atlas: str | None = None
    atlas_distance: int = 50
    datasets: list = field(default_factory=list)

    # ─── stage directories (run-root level) ──────────────────────────────
    @property
    def peak_calling_dir(self) -> Path:
        return self.output_dir / "01_peak_calling"

    @property
    def cb_filter_dir(self) -> Path:
        return self.output_dir / "02_cb_filter"

    @property
    def gtf_annotation_dir(self) -> Path:
        return self.output_dir / "03_gtf_annotation"

    @property
    def pas_gene_dir(self) -> Path:
        return self.output_dir / "04_pas_gene_assignment"

    @property
    def annotated_dir(self) -> Path:
        return self.output_dir / "05_annotated_matrix"

    @property
    def preprocessing_dir(self) -> Path:
        return self.output_dir / "06_preprocessing"

    @property
    def clustering_dir(self) -> Path:
        return self.output_dir / "07_clustering"

    @property
    def differential_dir(self) -> Path:
        return self.output_dir / "08_differential"

    @property
    def gtf_cache_dir(self) -> Path:
        return self.output_dir / "gtf_cache"

    @property
    def figures_dir(self) -> Path:
        return self.output_dir / "figures"

    @property
    def per_dataset_dir(self) -> Path:
        return self.output_dir / "per_dataset"

    @property
    def switch_dir(self) -> Path:
        return self.output_dir / "switch_out"

    @property
    def peakcalling_legacy_dir(self) -> Path:
        return self.output_dir / "peakcalling"

    # ─── run-root files ──────────────────────────────────────────────────
    @property
    def run_config_json(self) -> Path:
        return self.output_dir / "run_config.json"

    @property
    def resources_jsonl(self) -> Path:
        return self.output_dir / "resources.jsonl"

    # ─── flat (run-root) data files ──────────────────────────────────────
    # Kept for callsites that still use them until per-dataset migration
    # completes.

    @property
    def posbed(self) -> Path:
        return self.output_dir / "posbed.bed"

    @property
    def negbed(self) -> Path:
        return self.output_dir / "negbed.bed"

    @property
    def pasbed(self) -> Path:
        return self.output_dir / "pasbed.bed"

    @property
    def filtered_cb(self) -> Path:
        return self.output_dir / "filtered_cb.tsv"

    @property
    def annotated_matrix(self) -> Path:
        return self.output_dir / "annotated_matrix.mtx"

    @property
    def annotatedbed(self) -> Path:
        return self.output_dir / "annotatedpas.bed"

    @property
    def pas_geneid(self) -> Path:
        return self.output_dir / "pas_gene.tsv"

    @property
    def endbed(self) -> Path:
        return self.output_dir / "gene_end.bed"

    @property
    def raw_features(self) -> Path:
        return self.output_dir / "raw_feature.tsv"

    @property
    def utr_lengths(self) -> Path:
        return self.output_dir / "utr_lengths.tsv"

    @property
    def cluster_labels(self) -> Path:
        return self.output_dir / "clusterlabels.csv"

    # ─── legacy-named properties (stale filenames, kept for callsite compat)
    # These are preserved until all callsites migrate to per-dataset accessors.

    @property
    def posmatrixpath(self) -> Path:
        return self.output_dir / "posmatrix.mtx"

    @property
    def negmatrixpath(self) -> Path:
        return self.output_dir / "negmatrix.mtx"

    @property
    def matrixpath(self) -> Path:
        return self.output_dir / "pascountmatrix.mtx"

    @property
    def filterd_matrix(self) -> Path:
        return self.output_dir / "filterdmatrix.mtx"

    @property
    def cluster_output(self) -> Path:
        return self.output_dir / "write" / "pbmc3k.h5ad"

    # ─── per-dataset accessors ───────────────────────────────────────────

    def dataset_dir(self, ds: str) -> Path:
        return self.per_dataset_dir / ds

    def pasbed_for(self, ds: str) -> Path:
        return self.dataset_dir(ds) / "pasbed.bed"

    def posbed_for(self, ds: str) -> Path:
        return self.dataset_dir(ds) / "posbed.bed"

    def negbed_for(self, ds: str) -> Path:
        return self.dataset_dir(ds) / "negbed.bed"

    def clusters_h5ad_for(self, ds: str) -> Path:
        return self.dataset_dir(ds) / "clusters.h5ad"

    def preprocessed_h5ad_for(self, ds: str) -> Path:
        return self.dataset_dir(ds) / "preprocessed.h5ad"

    def filtered_cb_for(self, ds: str) -> Path:
        return self.dataset_dir(ds) / "filtered_cb.tsv"

    def pas_gene_for(self, ds: str) -> Path:
        return self.dataset_dir(ds) / "pas_gene.tsv"

    def annotatedpas_for(self, ds: str) -> Path:
        return self.dataset_dir(ds) / "annotatedpas.bed"

    def annotated_matrix_for(self, ds: str) -> Path:
        return self.dataset_dir(ds) / "annotated_matrix.mtx"

    def annotated_pas_ids_for(self, ds: str) -> Path:
        return self.dataset_dir(ds) / "annotated_pas_ids.tsv"

    def annotated_cells_for(self, ds: str) -> Path:
        return self.dataset_dir(ds) / "annotated_cells.tsv"

    def raw_dir_for(self, ds: str) -> Path:
        return self.dataset_dir(ds) / "raw"

    # ─── helpers ─────────────────────────────────────────────────────────

    def stage_stats(self, stage: str) -> Path:
        return getattr(self, f"{stage}_dir") / f"{stage}_stats.json"

    def get_fisher_dir(self, c1, c2) -> Path:
        return self.switch_dir / f"fisherresults_{c1}_{c2}.tsv"

    def setup(self) -> None:
        """Create all standard pipeline directories under output_dir."""
        for p in (
            self.peak_calling_dir,
            self.cb_filter_dir,
            self.gtf_annotation_dir,
            self.pas_gene_dir,
            self.annotated_dir,
            self.preprocessing_dir,
            self.clustering_dir,
            self.differential_dir,
            self.gtf_cache_dir,
            self.per_dataset_dir,
            self.figures_dir,
            self.switch_dir,
        ):
            p.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Module-level singleton
#
# Importers keep writing:
#   from ema.config import directory_config
# but ``directory_config`` is now an *instance* of DirectoryConfig, not the
# class itself.  Because the dataclass is frozen, ``set_directory_config``
# uses ``object.__setattr__`` to mutate fields on the *same* object in-place.
# Every ``from ema.config import directory_config`` binding still points at
# this object, so mutations are immediately visible through existing bindings
# (no re-import needed).
# ---------------------------------------------------------------------------

directory_config = DirectoryConfig()


def set_directory_config(**kwargs) -> None:
    """Mutate the module-level singleton's fields in-place.

    Used by main.py to update ``output_dir`` / ``gtf_dir`` / etc. once the
    timestamped run root is known.  Bypasses ``frozen=True`` via
    ``object.__setattr__`` — this is the one place mutation is allowed.
    Derived properties recompute on every access, so changing ``output_dir``
    cascades to every path attribute automatically.

    Example::

        from ema.config import set_directory_config
        from pathlib import Path
        set_directory_config(output_dir=Path("/path/to/run"), datasets=[...])
    """
    for k, v in kwargs.items():
        object.__setattr__(directory_config, k, v)


# ---------------------------------------------------------------------------
# Legacy dataclasses — kept unchanged (variable_config, filter_config).
# They still reference ``args`` which is now a lazy proxy; the first attribute
# access on ``args`` triggers cli() and populates the namespace.
# ---------------------------------------------------------------------------

@dataclass
class variable_config:
    seqlen = args.seqlen
    cb_len = args.cb_len
    default_threshold = 5
    merge_len = 100
    ignore_chro = ["MT", "mt"]
    barcode_tag = args.barcode_tag
    time = 0
    matrixmarketheader = f"%%MatrixMarket matrix coordinate integer general\n"

@dataclass
class filter_config:
    # Class defaults must match ema/cli/defaults.py::DEFAULTS so a bare
    # `ema run` invocation (no --config and no --min-* flags) produces the
    # same numbers regardless of which surface set them.  Pre-fix:
    # DEFAULTS["min-read"] was 1500 but the class attribute was 2000, so
    # the result you got depended on whether the YAML loader had run yet.
    # Centralisation in Part B will derive these from the schema.
    min_read = 1500
    min_cells = 3
    min_genes = 50
    min_pas_per_cell = args.min_pas_per_cell
