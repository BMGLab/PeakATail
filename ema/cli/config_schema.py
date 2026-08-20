"""Single source of truth for every parameter ``ema run`` consumes.

This module defines :class:`RunConfig`, a dataclass whose fields hold the
canonical defaults, types, choices, CLI flag names, and YAML keys for
*every* parameter the pipeline accepts.  Click options, the wizard,
the YAML loader, and the pipeline body all derive from this single
source so the surfaces cannot drift.

Adding a new flag is a one-liner: append a new ``field()`` entry to
:class:`RunConfig` with the right :class:`FieldSpec` metadata, and
Click options + YAML keys + the legacy bridge pick it up automatically.

Design constraints
------------------
* Fields are organised by section comments so the auto-generated --help
  groups stay readable.
* Field defaults must match the previously-shipped behaviour exactly so
  Part B is a refactor, not a behaviour change.  The atlas baseline in
  reports/baseline_apa_completeness/ is computed against these defaults;
  changing any default requires regenerating the baseline.
* The schema does NOT validate dataset structure (that is the YAML
  loader's job); it only knows about scalar/path parameters.

Why a dataclass and not a TypedDict?
------------------------------------
* Dataclasses give us ``fields(cls)`` reflection for free.
* :func:`dataclasses.field` accepts ``metadata=`` to attach the FieldSpec
  without bloating the type annotation.
* mypy / IDE autocompletion works on attribute access.
"""
from __future__ import annotations

from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Any, Callable, Iterable, Optional

import click


# ---------------------------------------------------------------------------
# Sentinel: marks "no default value supplied" without conflicting with None
# (which IS a valid default for optional path / str fields).
# ---------------------------------------------------------------------------
class _Unset:
    def __repr__(self) -> str:  # pragma: no cover -- diagnostic only
        return "<UNSET>"


UNSET = _Unset()


# ---------------------------------------------------------------------------
# Per-field metadata: tells the generators (Click, wizard, YAML, bridge)
# everything they need to know about a field.
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class FieldSpec:
    """Per-field metadata held in dataclass field() metadata.

    Attributes:
        cli_flag: CLI flag (e.g. ``"--min-read"``).  ``None`` means the
            field is not exposed on the CLI.
        yaml_key: YAML key (e.g. ``"min_read"``).  ``None`` means the
            field is not loadable from YAML.
        description: Human-readable help text used by Click ``--help``
            and the wizard.
        choice: Tuple of valid string values, or ``None`` for free-form.
        is_flag: True if the CLI flag is a boolean ``--xxx/--no-xxx``.
        legacy_alias: Optional legacy YAML key/CLI flag to accept for
            backward compatibility (e.g. ``min_genes`` → ``min_pas_per_cell``).
        click_type: Override the auto-derived Click type.  Use a
            ``click.Path`` instance for path fields that need
            ``exists=True``.
        click_kwargs: Extra keyword args forwarded to ``click.option``.
        legacy_args_attr: Name of the attribute on the legacy ``args``
            namespace this field should be bridged to (when the names
            differ -- e.g. ``peak_strategy`` → ``args.strategy``).
        legacy_dataclass_attr: ``"<which_dataclass>.<attr>"`` for fields
            that bridge into ``directory_config`` / ``variable_config`` /
            ``filter_config`` rather than the catch-all ``args`` namespace.
        skip_legacy_bridge: If True, the bridge does not touch the legacy
            globals (used for purely CLI-only fields like ``--config``).
        applies_to: Optional frozenset of subcommand names this field is
            relevant to (e.g. ``frozenset({"switch_diff"})``).  ``None``
            means the field applies everywhere (i.e. ``ema run`` + all
            subcommands).  Fields that don't apply to ``ema run`` should
            set ``skip_legacy_bridge=True`` so the bridge ignores them.
            This attribute is purely informational — the YAML loader and
            Click generators do NOT filter on it, so switch-only fields
            are harmless on ``ema run``.
    """
    cli_flag: Optional[str] = None
    yaml_key: Optional[str] = None
    description: str = ""
    choice: Optional[tuple[str, ...]] = None
    is_flag: bool = False
    legacy_alias: Optional[str] = None
    click_type: Any = None
    click_kwargs: Optional[dict] = None
    legacy_args_attr: Optional[str] = None
    legacy_dataclass_attr: Optional[str] = None
    skip_legacy_bridge: bool = False
    applies_to: Optional[frozenset] = None


def _spec(**kwargs: Any) -> dict[str, FieldSpec]:
    """Convenience factory: ``metadata={"spec": FieldSpec(**kwargs)}``."""
    return {"spec": FieldSpec(**kwargs)}


# ---------------------------------------------------------------------------
# RunConfig: the schema.  ORDER MATTERS for the auto-generated --help
# layout (Click renders flags top-to-bottom in declaration order).
# ---------------------------------------------------------------------------
@dataclass
class RunConfig:
    """Canonical parameter container for ``ema run``.

    Use :meth:`from_click_kwargs` to construct from a Click invocation,
    :meth:`from_yaml` to construct from a YAML file, and
    :meth:`apply_to_legacy_globals` to push the values into the legacy
    ``ema.config`` module-level dataclasses the pipeline body reads.

    Filename overrides
    ------------------
    The ``filenames`` dict lets users rename specific pipeline output files
    without touching ``DirectoryConfig``.  Supply it under the ``filenames:``
    YAML key; CLI override is not supported (use the YAML instead).

    Overrideable keys (all optional; shown with their defaults)::

        filenames:
          pasbed:             pasbed.bed
          posbed:             posbed.bed
          negbed:             negbed.bed
          filtered_cb:        filtered_cb.tsv
          annotated_matrix:   annotated_matrix.mtx
          annotatedbed:       annotatedpas.bed
          pas_geneid:         pas_gene.tsv
          endbed:             gene_end.bed
          raw_features:       raw_feature.tsv
          utr_lengths:        utr_lengths.tsv
          cluster_labels:     clusterlabels.csv
          clusters_h5ad:      clusters.h5ad
          preprocessed_h5ad:  preprocessed.h5ad
    """

    # ─── inputs / outputs ────────────────────────────────────────────────
    config: Optional[Path] = field(
        default=None,
        metadata=_spec(
            cli_flag="--config", yaml_key=None, skip_legacy_bridge=True,
            click_type=click.Path(exists=True, dir_okay=False, resolve_path=True),
            description="YAML config; CLI flags override individual keys.",
        ),
    )
    output: str = field(
        default="emaout",
        metadata=_spec(
            cli_flag="--output", yaml_key="output_dir",
            # IMPORTANT: do NOT bridge to directory_config.output_dir from
            # the schema -- run.py resolves the timestamp suffix and sets
            # directory_config.output_dir explicitly *before*
            # apply_to_legacy_globals runs.  Bridging here would clobber
            # the timestamped path back to the bare default.
            skip_legacy_bridge=True,
            click_type=click.Path(file_okay=False, resolve_path=True),
            description="Output directory (timestamp suffix added automatically).",
        ),
    )
    bam_dir: Optional[str] = field(
        default=None,
        metadata=_spec(
            cli_flag="--bam-dir", yaml_key=None,
            legacy_dataclass_attr="directory_config.bam_dir",
            click_type=click.Path(exists=True),
            description="Single-BAM convenience.",
        ),
    )
    bam_files: Optional[str] = field(
        default=None,
        metadata=_spec(
            cli_flag="--bam-files", yaml_key=None, skip_legacy_bridge=True,
            description="Comma-separated multi-BAM list.",
        ),
    )
    gtf: Optional[str] = field(
        default=None,
        metadata=_spec(
            cli_flag="--gtf", yaml_key="gtf",
            legacy_dataclass_attr="directory_config.gtf_dir",
            click_type=click.Path(exists=True, dir_okay=False),
            description="GTF annotation file.",
        ),
    )
    atlas: Optional[str] = field(
        default=None,
        metadata=_spec(
            cli_flag="--atlas", yaml_key="atlas",
            legacy_dataclass_attr="directory_config.atlas",
            click_type=click.Path(exists=True, dir_okay=False),
            description="Reference PAS atlas BED.",
        ),
    )
    atlas_distance: int = field(
        default=50,
        metadata=_spec(
            cli_flag="--atlas-distance", yaml_key="atlas_distance",
            legacy_dataclass_attr="directory_config.atlas_distance",
            description="Atlas snap distance (bp).",
        ),
    )
    atlas_mode: str = field(
        default="annotate",
        metadata=_spec(
            cli_flag="--atlas-mode", yaml_key="atlas_mode",
            choice=("annotate", "filter"),
            legacy_args_attr="atlas_mode",
            description=(
                "How atlas snapping handles PAS with no atlas match (only "
                "relevant when --atlas is set). 'annotate' (default) KEEPS "
                "every PAS -- non-matching PAS are alternative-PAS signal, "
                "not noise -- and records atlas_match/atlas_distance_bp on "
                "the PAS ledger + pasbed. 'filter' restores the pre-D9 "
                "behaviour of dropping PAS with no atlas hit within "
                "--atlas-distance."
            ),
        ),
    )
    filenames: Optional[dict] = field(
        default=None,
        metadata=_spec(
            cli_flag=None, yaml_key="filenames",
            skip_legacy_bridge=True,  # wired manually in main.py via set_directory_config
            description=(
                "Optional dict of output filename overrides "
                "(YAML only — see RunConfig docstring for valid keys)."
            ),
        ),
    )

    # ─── read processing ─────────────────────────────────────────────────
    seq_len: Optional[int] = field(
        default=None,
        metadata=_spec(
            cli_flag="--seq-len", yaml_key="seqlen",
            legacy_dataclass_attr="variable_config.seqlen",
            description="Sequencing read length.",
        ),
    )
    cb_len: Optional[int] = field(
        default=None,
        metadata=_spec(
            cli_flag="--cb-len", yaml_key="cb_len",
            legacy_dataclass_attr="variable_config.cb_len",
            description="Cell-barcode length (bp).",
        ),
    )
    barcode_tag: Optional[str] = field(
        default=None,
        metadata=_spec(
            cli_flag="--barcode-tag", yaml_key="barcode_tag",
            legacy_dataclass_attr="variable_config.barcode_tag",
            description="BAM tag holding the cell barcode (default CB).",
        ),
    )
    default_threshold: int = field(
        default=5,
        metadata=_spec(
            cli_flag="--default-threshold", yaml_key="default_threshold",
            legacy_dataclass_attr="variable_config.default_threshold",
            description="Minimum read support per peak (peak-calling threshold).",
        ),
    )
    merge_len: int = field(
        default=100,
        metadata=_spec(
            cli_flag="--merge-len", yaml_key="merge_len",
            legacy_dataclass_attr="variable_config.merge_len",
            description="Maximum gap (bp) to merge adjacent peaks.",
        ),
    )
    ignore_chro: str = field(
        default="MT,mt",
        metadata=_spec(
            cli_flag="--ignore-chro", yaml_key="ignore_chro",
            legacy_dataclass_attr="variable_config.ignore_chro",
            click_type=click.STRING,
            description=(
                "Chromosomes excluded from peak-calling. "
                "Comma-separated string for CLI (e.g. 'MT,mt,chrM'). "
                "YAML accepts a list or a comma-separated string."
            ),
        ),
    )

    # ─── concurrency / runtime ──────────────────────────────────────────
    threads: Optional[int] = field(
        default=None,
        metadata=_spec(
            cli_flag="--threads", yaml_key="threads", skip_legacy_bridge=True,
            description="Max parallel workers (auto-detected if not set).",
        ),
    )
    bam_threads: int = field(
        default=4,
        metadata=_spec(
            cli_flag="--bam-threads", yaml_key="bam_threads",
            legacy_args_attr="bam_threads",
            description="pysam decompression threads per BAM.",
        ),
    )
    pipeline: bool = field(
        default=False,
        metadata=_spec(
            cli_flag="--pipeline", yaml_key="pipeline", is_flag=True,
            legacy_args_attr="pipeline",
            description="Pipeline mode (currently only used by --tiles).",
        ),
    )
    batch_size: int = field(
        default=10000,
        metadata=_spec(
            cli_flag="--batch-size", yaml_key="batch_size",
            legacy_args_attr="batch_size",
            description="Worker batch size for streaming reads.",
        ),
    )
    tiles: bool = field(
        default=False,
        metadata=_spec(
            cli_flag="--tiles", yaml_key="tiles", is_flag=True,
            legacy_args_attr="tiles",
            description="Enable tile-based peak calling (parallel pool).",
        ),
    )
    tile_size: Optional[int] = field(
        default=None,
        metadata=_spec(
            cli_flag="--tile-size", yaml_key="tile_size",
            legacy_args_attr="tile_size",
            description="Tile size in bp (auto if unset).",
        ),
    )
    tile_overlap: int = field(
        default=10000,
        metadata=_spec(
            cli_flag="--tile-overlap", yaml_key="tile_overlap",
            legacy_args_attr="tile_overlap",
            description="Tile overlap in bp.",
        ),
    )

    # ─── peak calling ────────────────────────────────────────────────────
    peak_strategy: str = field(
        # lambda_gradient is the recommended production strategy
        # (40.6% precision @50bp on the benchmark; see technical report
        # Sections 3 + 5).  `original` remains available as a baseline but
        # is no longer the default — its lack of statistical filtering
        # inflates PAS counts (~2x on the v9 BAM) and trips users who
        # don't know to override.
        default="lambda_gradient",
        metadata=_spec(
            cli_flag="--peak-strategy", yaml_key="peak_strategy",
            legacy_args_attr="strategy",
            description="Peak-calling strategy (run --list-strategies to see).",
        ),
    )
    lambda_window: int = field(
        default=5000,
        metadata=_spec(
            cli_flag="--lambda-window", yaml_key="lambda_window",
            legacy_args_attr="lambda_window",
            description="Background lambda estimation window (bp).",
        ),
    )
    lambda_method: str = field(
        default="median",
        metadata=_spec(
            cli_flag="--lambda-method", yaml_key="lambda_method",
            legacy_args_attr="lambda_method",
            description="Lambda estimator (median / mean / ...).",
        ),
    )
    lambda_fold_change: float = field(
        default=2.0,
        metadata=_spec(
            cli_flag="--lambda-fold-change", yaml_key="lambda_fold_change",
            legacy_args_attr="lambda_fold_change",
            description="Lambda fold-change cutoff for peak detection.",
        ),
    )
    max_pas: int = field(
        default=5,
        metadata=_spec(
            cli_flag="--max-pas", yaml_key="max_pas",
            legacy_args_attr="max_pas",
            description="Maximum PAS sites kept per peak.",
        ),
    )
    smoothing_window: int = field(
        default=50,
        metadata=_spec(
            cli_flag="--smoothing-window", yaml_key="smoothing_window",
            legacy_args_attr="smoothing_window",
            description="Coverage smoothing window (bp).",
        ),
    )
    min_prominence: float = field(
        default=5.0,
        metadata=_spec(
            cli_flag="--min-prominence", yaml_key="min_prominence",
            legacy_args_attr="min_prominence",
            description="scipy.signal.find_peaks prominence threshold.",
        ),
    )
    dynamic_threshold: bool = field(
        default=False,
        metadata=_spec(
            cli_flag="--dynamic-threshold", yaml_key="dynamic_threshold",
            is_flag=True, legacy_args_attr="dynamic_threshold",
            description="Use a per-window dynamic peak threshold.",
        ),
    )
    floor_threshold: int = field(
        default=3,
        metadata=_spec(
            cli_flag="--floor-threshold", yaml_key="floor_threshold",
            legacy_args_attr="floor_threshold",
            description="Minimum peak height (clamps dynamic threshold).",
        ),
    )
    pas_gap: int = field(
        default=100,
        metadata=_spec(
            cli_flag="--pas-gap", yaml_key="pas_gap",
            legacy_args_attr="pas_gap",
            description="Minimum gap between PAS within a peak (bp).",
        ),
    )
    min_pas_spacing: int = field(
        default=-1,
        metadata=_spec(
            cli_flag="--min-pas-spacing", yaml_key="min_pas_spacing",
            legacy_dataclass_attr="variable_config.min_pas_spacing",
            description=(
                "Hard distance floor (bp) for the post-detection PAS merger. "
                "Adjacent PAS within one peak whose gap is smaller than this "
                "are merged unconditionally (Tier 1). "
                "-1 (default) auto-detects median read length from each BAM. "
                "0 disables the distance tier."
            ),
        ),
    )
    min_pas_prominence: float = field(
        default=5.0,
        metadata=_spec(
            cli_flag="--min-pas-prominence", yaml_key="min_pas_prominence",
            legacy_dataclass_attr="variable_config.min_pas_prominence",
            click_type=click.FLOAT,
            description=(
                "Static valley-depth threshold for the Tier-2 post-detection "
                "merger.  Used by NON-lambda strategies (original, "
                "sierra_iterative).  Lambda strategies (lambda_poisson, "
                "lambda_gradient) ignore this and use their own "
                "compute_lambda(heights) as the threshold.  "
                "Negative disables Tier 2 entirely.  Default 5.0."
            ),
        ),
    )
    cleavage_offset: int = field(
        default=0,
        metadata=_spec(
            cli_flag="--cleavage-offset", yaml_key="cleavage_offset",
            legacy_dataclass_attr="variable_config.cleavage_offset",
            description=(
                "3' cleavage-site offset correction (bp; issue #72).  Called "
                "peak 3' ends stop ~90-105 nt short of the true cleavage site "
                "because 10x R2 coverage runs out before the poly(A) junction. "
                "When > 0, the reported PAS 3' end is shifted downstream by "
                "this many bp after peak calling, so tight-cutoff benchmarks "
                "and atlas annotation use the inferred cleavage position. "
                "A sane data-driven constant is ~90-100 (try 95). "
                "0 (default) preserves legacy behaviour (no shift). "
                "Use --auto-cleavage-offset to estimate this from the data."
            ),
        ),
    )
    auto_cleavage_offset: bool = field(
        default=False,
        metadata=_spec(
            cli_flag="--auto-cleavage-offset", yaml_key="auto_cleavage_offset",
            is_flag=True,
            legacy_dataclass_attr="variable_config.auto_cleavage_offset",
            description=(
                "Data-driven 3' cleavage-offset estimation (issue #72).  When "
                "set, the offset is inferred per run from the called peaks and "
                "--genome-fasta (genomic A-fraction crest + AATAAA density "
                "downstream of each peak 3' end) instead of using the fixed "
                "--cleavage-offset constant, then applied the same way.  "
                "Requires --genome-fasta; falls back to ~95 bp if the profiles "
                "are inconclusive.  Off (default) preserves legacy behaviour."
            ),
        ),
    )

    # ─── filters (D6: wired into `ema run`; see ema/main.py::_apply_pas_filters) ──
    ip_filter: bool = field(
        default=False,
        metadata=_spec(
            cli_flag="--ip-filter", yaml_key="ip_filter", is_flag=True,
            legacy_args_attr="ip_filter",
            description=(
                "Enable internal-priming filter: drops PAS near genomic "
                "A-rich stretches (requires --genome-fasta). Applied to the "
                "pos/neg PAS BEDs before gene assignment."
            ),
        ),
    )
    genome_fasta: Optional[str] = field(
        default=None,
        metadata=_spec(
            cli_flag="--genome-fasta", yaml_key="genome_fasta",
            legacy_args_attr="genome_fasta",
            click_type=click.Path(exists=True),
            description="Genome FASTA (indexed with pyfaidx/.fai) required by --ip-filter.",
        ),
    )
    ip_filter_mode: str = field(
        default="annotate",
        metadata=_spec(
            cli_flag="--ip-filter-mode", yaml_key="ip_filter_mode",
            choice=("annotate", "filter"),
            legacy_args_attr="ip_filter_mode",
            description=(
                "How the internal-priming filter handles a flagged PAS "
                "(only relevant when --ip-filter is set). 'annotate' "
                "(default) KEEPS every PAS and records the internal_priming "
                "flag on the PAS ledger + pasbed. 'filter' restores the "
                "pre-D9 behaviour of dropping flagged PAS."
            ),
        ),
    )
    annot_filter: bool = field(
        default=False,
        metadata=_spec(
            cli_flag="--annot-filter", yaml_key="annot_filter", is_flag=True,
            legacy_args_attr="annot_filter",
            description=(
                "Enable annotation filter: keeps only PAS overlapping an "
                "annotated gene region. Uses --annotation-bed if given, "
                "otherwise the GTF-derived gene BED (requires --gtf)."
            ),
        ),
    )
    annotation_bed: Optional[str] = field(
        default=None,
        metadata=_spec(
            cli_flag="--annotation-bed", yaml_key="annotation_bed",
            legacy_args_attr="annotation_bed",
            click_type=click.Path(exists=True, dir_okay=False),
            description=(
                "Annotation BED for --annot-filter. Optional -- defaults to "
                "the GTF-derived gene BED (directory_config.endbed) when unset."
            ),
        ),
    )
    ip_a_stretch: int = field(
        default=6,
        metadata=_spec(
            cli_flag="--ip-a-stretch", yaml_key="ip_a_stretch",
            legacy_args_attr="ip_a_stretch",
            description="Minimum consecutive A's (or T's on - strand) flagged as internal priming.",
        ),
    )
    ip_a_fraction: float = field(
        default=0.7,
        metadata=_spec(
            cli_flag="--ip-a-fraction", yaml_key="ip_a_fraction",
            legacy_args_attr="ip_a_fraction",
            click_type=click.FLOAT,
            description="Alternative internal-priming trigger: A/T fraction within the window.",
        ),
    )
    ip_window_left: int = field(
        default=10,
        metadata=_spec(
            cli_flag="--ip-window-left", yaml_key="ip_window_left",
            legacy_args_attr="ip_window_left",
            description="Internal-priming check window, bp upstream of the PAS.",
        ),
    )
    ip_window_right: int = field(
        default=30,
        metadata=_spec(
            cli_flag="--ip-window-right", yaml_key="ip_window_right",
            legacy_args_attr="ip_window_right",
            description="Internal-priming check window, bp downstream of the PAS.",
        ),
    )

    # ─── read-level poly(A) evidence (Stage 1; ema/countmatrix/polya.py) ──
    polya_evidence: str = field(
        default="on",
        metadata=_spec(
            cli_flag="--polya-evidence", yaml_key="polya_evidence",
            choice=("on", "off"),
            legacy_args_attr="polya_evidence",
            description=(
                "Collect read-level poly(A) soft-clip evidence during peak "
                "calling (default on). Annotate-only by default: each PAS's "
                "clip-read support is written into BED column 5 of "
                "pasbed.bed (previously hardcoded 0); coordinates, counts "
                "and PAS selection are unchanged unless --polya-mode or the "
                "clip_seeded strategy says otherwise. 'off' restores the "
                "pre-Stage-1 byte-identical output."
            ),
        ),
    )
    polya_mode: str = field(
        default="annotate",
        metadata=_spec(
            cli_flag="--polya-mode", yaml_key="polya_mode",
            choice=("annotate", "filter", "require"),
            legacy_args_attr="polya_mode",
            description=(
                "How poly(A) support is enforced (only relevant when "
                "--polya-evidence is on). 'annotate' (default) keeps every "
                "PAS. 'filter' drops PAS with zero clip support (BED score "
                "0 == coverage-only tier) at the same seam as --ip-filter, "
                "before gene assignment. 'require' additionally FAILS the "
                "run when no clip-supported PAS remain — use it when an "
                "unsupported call set must never ship silently."
            ),
        ),
    )
    polya_min_clip: int = field(
        default=6,
        metadata=_spec(
            cli_flag="--polya-min-clip", yaml_key="polya_min_clip",
            legacy_args_attr="polya_min_clip",
            description=(
                "Minimum terminal soft-clip length AND minimum A/T run "
                "adjacent to the alignment boundary for a read to count as "
                "poly(A) evidence (default 6; measured 92x wrong-end "
                "specificity on PBMC)."
            ),
        ),
    )
    polya_min_purity: float = field(
        default=0.8,
        metadata=_spec(
            cli_flag="--polya-min-purity", yaml_key="polya_min_purity",
            legacy_args_attr="polya_min_purity",
            click_type=click.FLOAT,
            description=(
                "Minimum A (forward) / T (reverse) fraction over the "
                "clipped bases (default 0.8)."
            ),
        ),
    )
    polya_window: int = field(
        default=100,
        metadata=_spec(
            cli_flag="--polya-window", yaml_key="polya_window",
            legacy_args_attr="polya_window",
            description=(
                "Half-window in bp around a PAS's strand-aware 3' base "
                "within which clip reads count as support for that PAS "
                "(default 100 — the benchmark's matching cutoff)."
            ),
        ),
    )
    polya_seed_window: int = field(
        default=25,
        metadata=_spec(
            cli_flag="--polya-seed-window", yaml_key="polya_seed_window",
            legacy_args_attr="polya_seed_window",
            description=(
                "Single-linkage gap in bp for clustering clip sites into "
                "PAS candidates (clip_seeded strategy only; default 25)."
            ),
        ),
    )
    polya_min_reads: int = field(
        default=1,
        metadata=_spec(
            cli_flag="--polya-min-reads", yaml_key="polya_min_reads",
            legacy_args_attr="polya_min_reads",
            description=(
                "Minimum distinct molecules (UMI-deduplicated; reads "
                "without a UB tag count as one molecule each) for a clip "
                "cluster to be emitted as a tier-1 PAS (clip_seeded "
                "strategy only; default 1)."
            ),
        ),
    )
    min_pas_per_cell: int = field(
        default=50,
        metadata=_spec(
            cli_flag="--min-pas-per-cell", yaml_key="min_pas_per_cell",
            legacy_alias="min_genes",
            legacy_dataclass_attr="filter_config.min_pas_per_cell",
            description="Minimum PAS per cell (also bridges to filter_config.min_genes).",
        ),
    )
    min_read: int = field(
        default=1500,
        metadata=_spec(
            cli_flag="--min-read", yaml_key="min_read",
            legacy_dataclass_attr="filter_config.min_read",
            description="Minimum reads per cell barcode.",
        ),
    )
    min_cells: int = field(
        default=3,
        metadata=_spec(
            cli_flag="--min-cells", yaml_key="min_cells",
            legacy_dataclass_attr="filter_config.min_cells",
            description="Minimum cells expressing a PAS.",
        ),
    )

    # ─── annotation ──────────────────────────────────────────────────────
    max_gene_distance: int = field(
        default=5000,
        metadata=_spec(
            cli_flag="--max-gene-distance", yaml_key="max_gene_distance",
            legacy_args_attr="max_gene_distance",
            description="Max distance for gene-end annotation (bp).",
        ),
    )
    utr_multiplier: float = field(
        default=2.0,
        metadata=_spec(
            cli_flag="--utr-multiplier", yaml_key="utr_multiplier",
            legacy_args_attr="utr_multiplier",
            description="3'UTR length multiplier for extended-3' annotation.",
        ),
    )
    include_extended: bool = field(
        default=False,
        metadata=_spec(
            cli_flag="--include-extended", yaml_key="include_extended",
            is_flag=True, legacy_args_attr="include_extended",
            description="Include extended-3' annotations.",
        ),
    )

    # ─── clustering ──────────────────────────────────────────────────────
    cluster_method: str = field(
        default="leiden_tfidf",
        metadata=_spec(
            cli_flag="--cluster-method", yaml_key="cluster_method",
            legacy_args_attr="clustering_method",
            description="Clustering strategy (leiden_tfidf / leiden_libsize / external).",
        ),
    )
    resolution: float = field(
        default=1.0,
        metadata=_spec(
            cli_flag="--resolution", yaml_key="resolution",
            legacy_args_attr="resolution",
            description="Leiden resolution.",
        ),
    )
    n_pcs: int = field(
        default=40,
        metadata=_spec(
            cli_flag="--n-pcs", yaml_key="n_pcs",
            legacy_args_attr="n_pcs",
            description="Number of principal components.",
        ),
    )
    external_clusters: Optional[str] = field(
        default=None,
        metadata=_spec(
            cli_flag="--external-clusters", yaml_key="external_clusters",
            legacy_args_attr="external_clusters",
            click_type=click.Path(exists=True),
            description="Path to external cluster labels TSV.",
        ),
    )
    random_seed: int = field(
        default=42,
        metadata=_spec(
            cli_flag="--random-seed", yaml_key="random_seed",
            legacy_args_attr="random_seed",
            description="RNG seed for clustering reproducibility.",
        ),
    )

    # ─── cross-dataset matching ─────────────────────────────────────────
    match_method: str = field(
        default="marker_overlap",
        metadata=_spec(
            cli_flag="--match-method", yaml_key="match_method",
            legacy_alias="cluster_match_method",
            legacy_args_attr="cluster_match_method",
            description="Cross-dataset cluster match strategy.",
        ),
    )
    n_top_markers: int = field(
        default=50,
        metadata=_spec(
            cli_flag="--n-top-markers", yaml_key="n_top_markers",
            legacy_args_attr="n_top_markers",
            description="Number of top marker PAS per cluster.",
        ),
    )

    # ─── switch diff parameters ──────────────────────────────────────────
    # These fields are consumed by `ema switch diff`.  They live in RunConfig
    # so the YAML loader recognises them (they appear in _LIVE_KEYS) and so
    # defaults are schema-derived instead of duplicated in _SUBCOMMAND_DEFAULTS.
    # skip_legacy_bridge=True keeps them out of the `ema run` legacy bridge.
    marker_top_n: int = field(
        default=200,
        metadata=_spec(
            cli_flag="--marker-top-n", yaml_key="marker_top_n",
            skip_legacy_bridge=True,
            description="Top-N markers per cluster for differential APA.",
            applies_to=frozenset({"switch_diff"}),
        ),
    )
    marker_method: str = field(
        default="wilcoxon",
        metadata=_spec(
            cli_flag="--marker-method", yaml_key="marker_method",
            skip_legacy_bridge=True,
            description="Marker ranking method for differential APA (wilcoxon / t-test ...).",
            applies_to=frozenset({"switch_diff"}),
        ),
    )
    fdr: float = field(
        default=0.05,
        metadata=_spec(
            cli_flag="--fdr", yaml_key="fdr",
            skip_legacy_bridge=True,
            description="FDR threshold for differential APA calls.",
            applies_to=frozenset({"switch_diff"}),
        ),
    )
    per_worker_mb: int = field(
        default=300,
        metadata=_spec(
            cli_flag="--per-worker-mb", yaml_key="per_worker_mb",
            skip_legacy_bridge=True,
            description="Memory cap (MB) per parallel worker in switch diff / length.",
            applies_to=frozenset({"switch_diff", "switch_length"}),
        ),
    )

    # ─── switch length parameters ────────────────────────────────────────
    isoform_agg: str = field(
        default="per_gene",
        metadata=_spec(
            cli_flag="--isoform-agg", yaml_key="isoform_agg",
            skip_legacy_bridge=True,
            choice=("per_gene", "within_utr", "between_utr", "per_isoform"),
            description=(
                "Aggregation/scope level. For `switch length` (PDUI): "
                "per_gene / per_isoform. For `switch diff`: per_gene (default; "
                "each PAS tested against the rest of its GENE) / within_utr "
                "(each PAS tested against the other PAS sharing its 3'UTR "
                "isoform -- tandem-UTR APA; `per_isoform` is accepted as a "
                "legacy alias) / between_utr (collapse PAS to 3'UTR-level "
                "counts and test differential 3'UTR PREFERENCE between "
                "groups; genes with >=2 UTRs only). Combine with "
                "--cluster-key to contrast ANY obs column (stage, celltype, "
                "leiden, ...) and --cluster-pairs to select specific pairs."
            ),
            applies_to=frozenset({"switch_diff", "switch_length"}),
        ),
    )
    isoform_collapse: str = field(
        default="none",
        metadata=_spec(
            cli_flag="--isoform-collapse", yaml_key="isoform_collapse",
            skip_legacy_bridge=True,
            choice=("none", "mean", "majority"),
            description="How to collapse isoforms when isoform_agg=per_gene.",
            applies_to=frozenset({"switch_length"}),
        ),
    )
    pdui_pseudocount: float = field(
        default=0.0,
        metadata=_spec(
            cli_flag="--pdui-pseudocount", yaml_key="pdui_pseudocount",
            skip_legacy_bridge=True,
            description=(
                "Pseudocount added to per-cell PAS counts before PDUI / entropy "
                "computation (classic / proportion / shannon). Default 0.0 preserves "
                "original behaviour; raise to e.g. 1.0 to avoid NaN on zero-count cells."
            ),
            applies_to=frozenset({"switch_length"}),
        ),
    )
    utr_unmatched: str = field(
        default="gene",
        metadata=_spec(
            cli_flag="--utr-unmatched", yaml_key="utr_unmatched",
            skip_legacy_bridge=True,
            choice=("drop", "gene"),
            description=(
                "How to handle PAS that overlap no annotated UTR under "
                "isoform_agg=per_isoform (switch length) or "
                "within_utr/between_utr (switch diff). 'gene' (default) KEEPS "
                "them via a gene-level fallback (UTR-agnostic); 'drop' "
                "restores the old drop-if-no-UTR behaviour."
            ),
            applies_to=frozenset({"switch_diff", "switch_length"}),
        ),
    )

    # ─── switch diff volcano parameters ──────────────────────────────────
    log2fc_thresh: float = field(
        default=1.0,
        metadata=_spec(
            cli_flag="--log2fc-thresh", yaml_key="log2fc_thresh",
            skip_legacy_bridge=True,
            description=(
                "log2 fold-change threshold drawn on the volcano plot. "
                "Also used to shade the 'significant' region. Default 1.0."
            ),
            applies_to=frozenset({"diff"}),
        ),
    )

    # ─── switch diff / switch test parameters ────────────────────────────
    min_cells_per_group: int = field(
        default=10,
        metadata=_spec(
            cli_flag="--min-cells-per-group", yaml_key="min_cells_per_group",
            skip_legacy_bridge=True,
            description=(
                "Minimum cells (with non-zero counts for NB strategies) in each "
                "cluster group for a PAS to be included in differential testing. "
                "Default 10."
            ),
            applies_to=frozenset({"diff"}),
        ),
    )

    count_mode: str = field(
        default="cells",
        metadata=_spec(
            cli_flag="--count-mode", yaml_key="count_mode",
            skip_legacy_bridge=True,
            description=(
                "Unit the fisher strategy aggregates for its 2xN contingency "
                "table. 'cells' (default, D4, calibrated): each cell "
                "contributes at most once via per-cell PAS detection among "
                "gene-expressing cells, de-pseudoreplicating the test. 'reads' "
                "(legacy, opt-in): sum read/UMI counts per group -- reads "
                "within a cell are correlated, so read-level fisher q-values "
                "are NOT FDR-calibrated (a permutation null reports q<0.05 "
                "hits in 100% of runs; see issue #74) and should be treated as "
                "a ranking screen only. The default was flipped reads->cells "
                "in issue #74 so the out-of-the-box path is calibrated; pass "
                "'reads' explicitly only for backward comparison. Ignored by "
                "the NB strategies, which model per-cell overdispersion "
                "directly."
            ),
            applies_to=frozenset({"diff"}),
        ),
    )

    # ─── cross-dataset MNN parameters ────────────────────────────────────
    mnn_components: int = field(
        default=30,
        metadata=_spec(
            cli_flag="--mnn-components", yaml_key="mnn_components",
            skip_legacy_bridge=True,
            description=(
                "Number of LSI/SVD components for the MNN shared embedding "
                "(switch match --strategy mnn). Default 30."
            ),
            applies_to=frozenset({"switch_match"}),
        ),
    )
    mnn_k_neighbors: int = field(
        default=10,
        metadata=_spec(
            cli_flag="--mnn-k-neighbors", yaml_key="mnn_k_neighbors",
            skip_legacy_bridge=True,
            description=(
                "Number of nearest neighbours for MNN search "
                "(switch match --strategy mnn). Default 10."
            ),
            applies_to=frozenset({"switch_match"}),
        ),
    )

    # ─── clustering fine-tuning (leiden_tfidf) ───────────────────────────
    n_neighbors: int = field(
        default=30,
        metadata=_spec(
            cli_flag="--n-neighbors", yaml_key="n_neighbors",
            legacy_args_attr="n_neighbors",
            description=(
                "Number of nearest neighbours for the kNN graph used by Leiden "
                "(leiden_tfidf default: 30; leiden_libsize default: 10 — set "
                "--n-neighbors explicitly to override)."
            ),
        ),
    )
    tfidf_scale_factor: float = field(
        default=1e4,
        metadata=_spec(
            cli_flag="--tfidf-scale-factor", yaml_key="tfidf_scale_factor",
            legacy_args_attr="tfidf_scale_factor",
            description=(
                "Scale factor for Signac Method 1 TF-IDF (leiden_tfidf strategy). "
                "Default 10000. Adjust if your counts have very different dynamic range."
            ),
        ),
    )
    depth_corr_threshold: float = field(
        default=0.75,
        metadata=_spec(
            cli_flag="--depth-corr-threshold", yaml_key="depth_corr_threshold",
            legacy_args_attr="depth_corr_threshold",
            description=(
                "Pearson |r| threshold for removing LSI components correlated with "
                "sequencing depth (leiden_tfidf, ArchR-style). Default 0.75. "
                "Set to 1.0 to disable depth-correlation filtering."
            ),
        ),
    )
    n_svd_components: int = field(
        default=50,
        metadata=_spec(
            cli_flag="--n-svd-components", yaml_key="n_svd_components",
            legacy_args_attr="n_svd_components",
            description=(
                "Number of SVD/PCA components computed before filtering/neighbor "
                "graph (leiden_tfidf and leiden_libsize strategies). Default 50."
            ),
        ),
    )
    n_top_hvg: int = field(
        default=2000,
        metadata=_spec(
            cli_flag="--n-top-hvg", yaml_key="n_top_hvg",
            legacy_args_attr="n_top_hvg",
            description=(
                "Number of highly variable genes/PAS selected before PCA "
                "(leiden_libsize strategy only). Default 2000."
            ),
        ),
    )

    # ─── validation (warn-only knobs, no behaviour) ────────────────────
    benchmark: bool = field(
        default=False,
        metadata=_spec(
            cli_flag="--benchmark", yaml_key="benchmark", is_flag=True,
            skip_legacy_bridge=True,
            description="No-op (use scripts/validate_strategies.py).",
        ),
    )
    validate_db: Optional[str] = field(
        default=None,
        metadata=_spec(
            cli_flag="--validate-db", yaml_key="validate_db",
            skip_legacy_bridge=True,
            click_type=click.Path(exists=True),
            description="No-op (use scripts/validate_strategies.py).",
        ),
    )

    # ─── runtime / observability (CLI-only) ────────────────────────────
    log_level: Optional[str] = field(
        default=None,
        metadata=_spec(
            cli_flag="--log-level", yaml_key=None, skip_legacy_bridge=True,
            description="Logger level or `name=LEVEL` (comma-separated).",
        ),
    )
    no_log_file: bool = field(
        default=False,
        metadata=_spec(
            cli_flag="--no-log-file", yaml_key=None, is_flag=True,
            skip_legacy_bridge=True,
            description="Don't write peakatail_<ts>.log next to outputs.",
        ),
    )
    no_progress: bool = field(
        default=False,
        metadata=_spec(
            cli_flag="--no-progress", yaml_key=None, is_flag=True,
            skip_legacy_bridge=True,
            description="Suppress Rich progress bars.",
        ),
    )
    verbose: int = field(
        default=0,
        metadata=_spec(
            cli_flag="--verbose", yaml_key=None, skip_legacy_bridge=True,
            description="Verbosity counter (-v / -vv).",
        ),
    )
    quiet: bool = field(
        default=False,
        metadata=_spec(
            cli_flag="--quiet", yaml_key=None, is_flag=True,
            skip_legacy_bridge=True,
            description="WARNING-and-up only.",
        ),
    )

    # ---------------------------------------------------------------------
    # Constructors
    # ---------------------------------------------------------------------
    @classmethod
    def from_yaml_dict(cls, cfg: dict[str, Any]) -> "RunConfig":
        """Build a RunConfig from a YAML-loaded dict.

        Honours :attr:`FieldSpec.legacy_alias` so old YAMLs (e.g. with
        ``min_genes`` instead of ``min_pas_per_cell``) keep working.
        Unknown keys are ignored here -- the YAML loader is responsible
        for warning on them.
        """
        kwargs: dict[str, Any] = {}
        yaml_to_field = yaml_key_to_field_name(cls)
        legacy_to_field = legacy_alias_to_field_name(cls)
        for key, val in cfg.items():
            if key in yaml_to_field:
                field_name = yaml_to_field[key]
                # ignore_chro: YAML may supply a list; normalise to comma str.
                if field_name == "ignore_chro" and isinstance(val, list):
                    val = ",".join(str(v) for v in val)
                kwargs[field_name] = val
            elif key in legacy_to_field:
                kwargs[legacy_to_field[key]] = val
            # else: silent skip (loader will have warned)
        return cls(**kwargs)

    # ---------------------------------------------------------------------
    # Apply to legacy module-level config dataclasses
    # ---------------------------------------------------------------------
    def apply_to_legacy_globals(
        self,
        user_set: Iterable[str] | None = None,
    ) -> None:
        """Push field values into ``ema.config`` module-level dataclasses.

        This is the *only* bridge between the new schema and the legacy
        pipeline body.  The body still reads ``args.<attr>``,
        ``filter_config.<attr>``, ``directory_config.<attr>`` etc.

        Args:
            user_set: Set of field names the user explicitly supplied
                (from Click ParameterSource).  When ``None`` we apply
                every non-default field; when given, only those fields
                are bridged so values left at the dataclass default do
                not trample YAML-set globals.
        """
        from ema import config as _cfg
        targets = {
            "directory_config": _cfg.directory_config,
            "variable_config": _cfg.variable_config,
            "filter_config": _cfg.filter_config,
            "args": _cfg.args,
        }
        for f in fields(self):
            spec: FieldSpec = f.metadata.get("spec")  # type: ignore[assignment]
            if spec is None or spec.skip_legacy_bridge:
                continue
            value = getattr(self, f.name)
            if user_set is not None and f.name not in user_set:
                # Skip fields the user didn't set so YAML wins.
                continue
            if spec.legacy_dataclass_attr:
                holder, attr = spec.legacy_dataclass_attr.split(".")
                # ignore_chro is stored as a comma-separated string in
                # RunConfig but variable_config.ignore_chro must be a list.
                if attr == "ignore_chro" and isinstance(value, str):
                    value = [c.strip() for c in value.split(",") if c.strip()]
                if holder == "directory_config":
                    # directory_config is a frozen dataclass — go through the
                    # set_directory_config helper, which uses object.__setattr__
                    # internally to mutate the singleton in-place.
                    _cfg.set_directory_config(**{attr: value})
                else:
                    setattr(targets[holder], attr, value)
                # min_pas_per_cell is read by preprocessing() as
                # filter_config.min_genes; bridge that too.
                if spec.legacy_alias == "min_genes":
                    setattr(targets["filter_config"], "min_genes", value)
            elif spec.legacy_args_attr:
                setattr(targets["args"], spec.legacy_args_attr, value)


# ---------------------------------------------------------------------------
# Reflection helpers used by Click / wizard / YAML loader generators
# ---------------------------------------------------------------------------
def field_specs(cls: type) -> dict[str, FieldSpec]:
    """Return ``{field_name: FieldSpec}`` for every field on ``cls``."""
    out: dict[str, FieldSpec] = {}
    for f in fields(cls):
        spec = f.metadata.get("spec")
        if spec is not None:
            out[f.name] = spec
    return out


def yaml_keys_from_schema(cls: type = RunConfig) -> set[str]:
    """All YAML keys the schema accepts.

    Used by :mod:`ema.cli.yaml_loader` to drive the unknown-key warning.
    """
    keys: set[str] = set()
    for spec in field_specs(cls).values():
        if spec.yaml_key is not None:
            keys.add(spec.yaml_key)
        if spec.legacy_alias is not None:
            keys.add(spec.legacy_alias)
    return keys


def yaml_key_to_field_name(cls: type = RunConfig) -> dict[str, str]:
    out: dict[str, str] = {}
    for name, spec in field_specs(cls).items():
        if spec.yaml_key is not None:
            out[spec.yaml_key] = name
    return out


def legacy_alias_to_field_name(cls: type = RunConfig) -> dict[str, str]:
    out: dict[str, str] = {}
    for name, spec in field_specs(cls).items():
        if spec.legacy_alias is not None:
            out[spec.legacy_alias] = name
    return out


def cli_defaults(cls: type = RunConfig) -> dict[str, Any]:
    """Return ``{cli_flag_kebab: default}`` (e.g. ``{"min-read": 1500}``).

    Used by :mod:`ema.cli.defaults` so DEFAULTS is derived from the
    schema instead of being a hand-maintained second copy.
    """
    out: dict[str, Any] = {}
    for f in fields(cls):
        spec: FieldSpec = f.metadata.get("spec")  # type: ignore[assignment]
        if spec is None or spec.cli_flag is None:
            continue
        kebab = spec.cli_flag.lstrip("-")
        out[kebab] = f.default
    return out


def click_options_from_schema(
    cls: type = RunConfig,
    skip: Iterable[str] = (),
) -> Callable:
    """Decorator factory that adds ``@click.option(...)`` for every
    schema field with a ``cli_flag``.

    Use ``skip`` to omit fields handled separately (e.g. ``--config`` /
    ``--output`` / common-options that ``common_options()`` already
    attaches).

    Click runs decorators in reverse of declaration order, so we apply
    in reverse to preserve schema ordering in --help.
    """
    import click

    skip_set = set(skip)

    def decorator(fn: Callable) -> Callable:
        opts = []
        for f in fields(cls):
            spec: FieldSpec = f.metadata.get("spec")  # type: ignore[assignment]
            if spec is None or spec.cli_flag is None or f.name in skip_set:
                continue
            opt_kwargs: dict = {"default": f.default, "help": spec.description}
            if spec.is_flag:
                opt_kwargs["is_flag"] = True
            elif spec.choice is not None:
                opt_kwargs["type"] = click.Choice(list(spec.choice))
            elif spec.click_type is not None:
                opt_kwargs["type"] = spec.click_type
            else:
                # Derive type from the dataclass annotation default
                ann = f.type
                if ann in (int, "int", "Optional[int]", "int | None"):
                    opt_kwargs["type"] = int
                elif ann in (float, "float", "Optional[float]"):
                    opt_kwargs["type"] = float
                elif ann is bool or ann == "bool":
                    opt_kwargs["is_flag"] = True
                # else: leave Click to default to str
            if spec.click_kwargs:
                opt_kwargs.update(spec.click_kwargs)
            # Variable name in callback kwargs is the field name (snake_case)
            opts.append(click.option(spec.cli_flag, f.name, **opt_kwargs))

        for opt in reversed(opts):
            fn = opt(fn)
        return fn

    return decorator


def wizard_prompts_from_schema(
    cls: type = RunConfig,
    fields_to_ask: Iterable[str] = (),
) -> "list[tuple[str, FieldSpec, Any]]":
    """Return ``[(field_name, spec, default), ...]`` for the wizard.

    The wizard iterates this list and builds a ``questionary`` prompt per
    entry: ``confirm`` for is_flag, ``select`` for choice, ``text``
    otherwise.  Centralising the iteration here means new fields show up
    in the wizard's advanced section automatically.
    """
    out: list[tuple[str, FieldSpec, Any]] = []
    field_map = {f.name: f for f in fields(cls)}
    for name in fields_to_ask:
        if name not in field_map:
            raise KeyError(f"unknown RunConfig field: {name}")
        spec = field_specs(cls)[name]
        out.append((name, spec, field_map[name].default))
    return out
