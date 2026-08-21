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
    filenames: dict = field(default_factory=dict)

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
        return self.output_dir / self.filenames.get("posbed", "posbed.bed")

    @property
    def negbed(self) -> Path:
        return self.output_dir / self.filenames.get("negbed", "negbed.bed")

    @property
    def pasbed(self) -> Path:
        return self.output_dir / self.filenames.get("pasbed", "pasbed.bed")

    @property
    def filtered_cb(self) -> Path:
        return self.output_dir / self.filenames.get("filtered_cb", "filtered_cb.tsv")

    @property
    def annotated_matrix(self) -> Path:
        return self.output_dir / self.filenames.get("annotated_matrix", "annotated_matrix.mtx")

    @property
    def annotatedbed(self) -> Path:
        return self.output_dir / self.filenames.get("annotatedbed", "annotatedpas.bed")

    @property
    def pas_geneid(self) -> Path:
        return self.output_dir / self.filenames.get("pas_geneid", "pas_gene.tsv")

    @property
    def endbed(self) -> Path:
        return self.output_dir / self.filenames.get("endbed", "gene_end.bed")

    @property
    def raw_features(self) -> Path:
        return self.output_dir / self.filenames.get("raw_features", "raw_feature.tsv")

    @property
    def utr_lengths(self) -> Path:
        return self.output_dir / self.filenames.get("utr_lengths", "utr_lengths.tsv")

    @property
    def cluster_labels(self) -> Path:
        return self.output_dir / self.filenames.get("cluster_labels", "clusterlabels.csv")

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
    # Data files live inside their numbered stage directory, nested by ds:
    #
    #   01_peak_calling/<ds>/pasbed.bed, posbed.bed, negbed.bed, raw/
    #   02_cb_filter/<ds>/filtered_cb.tsv
    #   03_gtf_annotation/<ds>/annotatedpas.bed
    #   04_pas_gene_assignment/<ds>/pas_gene.tsv
    #   05_annotated_matrix/<ds>/annotated_matrix.mtx + index files
    #   06_preprocessing/<ds>/preprocessed.h5ad
    #   07_clustering/<ds>/clusters.h5ad

    def pasbed_for(self, ds: str) -> Path:
        return self.peak_calling_dir / ds / self.filenames.get("pasbed", "pasbed.bed")

    def posbed_for(self, ds: str) -> Path:
        return self.peak_calling_dir / ds / self.filenames.get("posbed", "posbed.bed")

    def negbed_for(self, ds: str) -> Path:
        return self.peak_calling_dir / ds / self.filenames.get("negbed", "negbed.bed")

    def clusters_h5ad_for(self, ds: str) -> Path:
        return self.clustering_dir / ds / self.filenames.get("clusters_h5ad", "clusters.h5ad")

    def preprocessed_h5ad_for(self, ds: str) -> Path:
        return self.preprocessing_dir / ds / self.filenames.get("preprocessed_h5ad", "preprocessed.h5ad")

    def filtered_cb_for(self, ds: str) -> Path:
        return self.cb_filter_dir / ds / self.filenames.get("filtered_cb", "filtered_cb.tsv")

    def pas_gene_for(self, ds: str) -> Path:
        return self.pas_gene_dir / ds / self.filenames.get("pas_geneid", "pas_gene.tsv")

    def annotatedpas_for(self, ds: str) -> Path:
        return self.gtf_annotation_dir / ds / self.filenames.get("annotatedbed", "annotatedpas.bed")

    def annotated_matrix_for(self, ds: str) -> Path:
        return self.annotated_dir / ds / self.filenames.get("annotated_matrix", "annotated_matrix.mtx")

    def annotated_pas_ids_for(self, ds: str) -> Path:
        return self.annotated_dir / ds / "annotated_pas_ids.tsv"

    def annotated_cells_for(self, ds: str) -> Path:
        return self.annotated_dir / ds / "annotated_cells.tsv"

    def raw_dir_for(self, ds: str) -> Path:
        return self.peak_calling_dir / ds / "raw"

    # ─── helpers ─────────────────────────────────────────────────────────

    def stage_stats(self, stage: str) -> Path:
        return getattr(self, f"{stage}_dir") / f"{stage}_stats.json"

    def get_fisher_dir(self, c1, c2) -> Path:
        return self.switch_dir / f"fisherresults_{c1}_{c2}.tsv"

    def setup(self) -> None:
        """Create all standard pipeline stage directories under output_dir.

        Per-dataset subdirectories inside each stage (01_peak_calling/<ds>/,
        07_clustering/<ds>/, etc.) are created on-demand by writers via
        ``path.parent.mkdir(parents=True, exist_ok=True)`` — no pre-creation
        needed here.
        """
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
    # Post-detection PAS-summit merger (strategy-agnostic; applied at the
    # peak_calling caller layer after every strategy.find_pas() call).
    # Tier 1: hard distance floor in bp.  -1 = auto-detect median read length
    # per BAM; 0 disables the distance tier.
    min_pas_spacing = -1
    # Tier 2: static fallback valley threshold (coverage units) used by
    # non-lambda strategies (original, sierra_iterative).  Lambda strategies
    # ignore this and use compute_lambda(heights) instead.  Negative
    # disables Tier 2 entirely.
    min_pas_prominence = 5.0
    # 3' cleavage-site offset correction (issue #72).  Called peak 3' ends
    # stop ~90-105 nt short of the true cleavage site (10x R2 coverage runs
    # out before the poly(A) junction).  When > 0, the reported PAS 3' end is
    # shifted downstream by this many bp after peak calling so tight-cutoff
    # benchmarks and atlas annotation use the inferred cleavage position.
    # 0 (default) = legacy behaviour (no shift).  See
    # ema/countmatrix/cleavage_offset.py.
    cleavage_offset = 0
    # When True, cleavage_offset is estimated per run from the called peaks +
    # genome FASTA instead of the fixed constant (issue #72; --auto-cleavage-
    # offset).  False (default) = use the cleavage_offset constant as-is.
    auto_cleavage_offset = False
    # Cache populated once per BAM by peak_calling when min_pas_spacing == -1.
    # Keyed by str(bam_path) -> int median read length.  Plain class-level
    # dict (not a dataclass field) so it's accessible on the class itself,
    # matching the access pattern variable_config.dataset_read_lengths.
    dataset_read_lengths = {}
    # --- peakAtail-prime: read acceptance geometry (--read-geometry) --------
    # "fixed" (v2) | "keep" | "true"; see ema.countmatrix.read.READ_GEOMETRIES.
    #
    # THIS IS THE BRANCH DEFAULT AND IT IS DELIBERATELY THE ONLY COPY OF IT.
    # An earlier revision left the legacy global at "fixed" and put the branch
    # default only in RunConfig, so `ema run` used one geometry and a direct
    # library call the other -- and `RunConfig.apply_to_legacy_globals()` then
    # leaked "true" into the process globals halfway through a pytest session,
    # making 15 tests order-dependent.  One default, one place.
    # RunConfig.read_geometry must carry the same literal; that is asserted by
    # tests/test_read_geometry.py::test_the_branch_default_is_single_valued.
    #
    # THE DEFAULT IS "fixed" BECAUSE THE MEASUREMENT SAID SO, not because the
    # change is unfinished.  `true` was implemented as the intended branch
    # default and then measured on two dev slices against
    # manuscript/24 3.1, which requires (i) dP@100 >= -0.005, (ii) dR_det
    # >= +0.010 and (iii) dF1 > 0 on every dataset:
    #   PBMC chr19+21   dP -0.0026 (i ok)   dR +0.0031 (ii FAILS)  dF1 +0.0030
    #   mouse1 ch18+19  dP -0.0292 (i FAILS by 5.9x)  dR +0.0063 (ii FAILS)
    # The mouse precision loss survives --read-exclude-flags 256 (-0.0117) and
    # the ablation puts essentially all of it on "stop discarding", not on
    # "stop fabricating the 3' end".  Per 24 3.1 a FAIL means the behaviour
    # stays behind a non-default flag.
    #
    # This is NOT the whole story and the flag is not dead: `true` is worth
    # +31.2 % (PBMC) / +23.5 % (mouse) of RAW COUNT-MATRIX MASS on those
    # slices -- 88.8 M reads genome-wide -- and 24 3.1 measures detection
    # only.  See results/prime/taskA_read_geometry_slice.tsv and the
    # CHANGELOG; the decision to move this literal belongs to whoever owns
    # the pre-registration, and needs a quantification criterion first.
    read_geometry = "fixed"
    # SAM flag mask vetoed on the COVERAGE channel (--read-exclude-flags).
    # 0 (default, = v2) means no filtering: read_check applies no `-F`, so a
    # multimapper contributes one coverage read PER ALIGNMENT.  3844 is
    # samtools' unmapped+secondary+qcfail+duplicate+supplementary.
    read_exclude_flags = 0

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
