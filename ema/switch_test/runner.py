"""Runner functions for the APA switch test (differential APA + PDUI).

Provides library-level ``run_diff`` and ``run_length`` entry points consumed
by ``ema/cli/switch_diff.py`` and ``ema/cli/switch_length.py``.  The old
argparse ``cli()`` has been removed; use ``ema switch diff`` / ``ema switch
length`` instead.
"""

from __future__ import annotations
import functools
import logging
import multiprocessing
from itertools import combinations
from pathlib import Path

log = logging.getLogger(__name__)

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from joblib import parallel_backend

from ema.switch_test.strategies import get_diff_strategy, list_diff_strategies
from ema.switch_test.pair_runner import run_one_pair
from ema.quantification.strategies import get_pdui_strategy, list_pdui_strategies
from ema.quantification.marker_selector import (
    select_marker_pas,
    restrict_count_matrix,
    save_markers,
)
from ema.utils import (
    COUNTS_LAYER,
    ResourceManager,
    get_resource_manager,
    reset_resource_manager,
)


def parse_cell_combinations(s: str) -> list[tuple[str, str]]:
    """'0,1;0,2;1,2' -> [('0','1'), ('0','2'), ('1','2')]"""
    out = []
    for pair in s.split(";"):
        parts = pair.split(",")
        if len(parts) != 2:
            raise ValueError(f"Each pair must be 'A,B', got: {pair}")
        out.append((parts[0].strip(), parts[1].strip()))
    return out


def _matrix_is_integral(matrix) -> bool:
    """True when every stored value of *matrix* is a whole number.

    Only the stored (non-zero) values of a sparse matrix are inspected --
    the implicit zeros are integral by definition.  ``NaN``/``inf`` count as
    non-integral: a count matrix has neither.

    Args:
        matrix: dense array or ``scipy.sparse`` matrix.

    Returns:
        ``True`` for integer/bool dtypes and for float dtypes whose every
        stored value equals its own rounding.
    """
    data = matrix.data if sp.issparse(matrix) else np.asarray(matrix)
    data = np.asarray(data)
    if data.size == 0:
        return True
    if data.dtype == bool or np.issubdtype(data.dtype, np.integer):
        return True
    if not np.issubdtype(data.dtype, np.floating):
        return False
    if not np.isfinite(data).all():
        return False
    return bool(np.all(data == np.rint(data)))


def _select_count_matrix(
    adata: ad.AnnData,
    layer: str | None = None,
    allow_non_count_matrix: bool = False,
):
    """Pick the matrix the switch tests should quantify, and prove it is counts.

    Resolution order:

    1. *layer* when the caller names one (the literal ``"X"`` forces
       ``adata.X``).  A named-but-absent layer is an error, never a silent
       fallback -- a silent fallback here is how TF-IDF weights reached the
       PDUI table in the first place.
    2. ``adata.layers['counts']`` when present -- the raw-count stash
       :func:`ema.clustering.clustering.clustering` writes before any
       ``normalize()`` runs.
    3. ``adata.X``.

    The chosen matrix is then required to be integral.  ``leiden_tfidf``
    overwrites ``.X`` in place with ``log1p(TF * IDF * scale_factor)`` and no
    raw copy was kept before the ``counts`` layer landed, so a clusters.h5ad
    written by an older build carries weights.  PDUI computed on those
    weights is not an approximation of the count answer -- measured against
    the same matrix re-derived from counts: r = 0.47 on informative cells,
    33.4% of cells cross the 0.5 boundary, and 10.6% of cluster-pair dPDUI
    signs flip.  Refuse rather than emit a plausible-looking number.

    Args:
        adata: the loaded AnnData.
        layer: explicit layer name, ``"X"``, or ``None`` for auto-selection.
        allow_non_count_matrix: downgrade the integrality failure to a
            warning.  Escape hatch for deliberately quantifying normalised
            input; never appropriate for a published number.

    Returns:
        ``(matrix, source_label)``.

    Raises:
        KeyError: *layer* was named but is absent.
        ValueError: the chosen matrix is non-integral and
            *allow_non_count_matrix* is False.
    """
    layers = getattr(adata, "layers", None)
    available = list(layers.keys()) if layers is not None else []

    if layer is not None and layer != "X":
        if layer not in available:
            raise KeyError(
                f"count layer {layer!r} is not present on this AnnData "
                f"(layers: {available or 'none'}). Pass --counts-layer X to "
                "quantify .X explicitly, or backfill the layer from the "
                "sibling 06_preprocessing/<ds>/preprocessed.h5ad."
            )
        matrix, source = layers[layer], f"layers[{layer!r}]"
    elif layer == "X":
        matrix, source = adata.X, ".X"
    elif COUNTS_LAYER in available:
        matrix, source = layers[COUNTS_LAYER], f"layers[{COUNTS_LAYER!r}]"
    else:
        matrix, source = adata.X, ".X"

    if not _matrix_is_integral(matrix):
        msg = (
            f"count matrix taken from {source} is NOT integral -- it holds "
            "normalised/transformed values, not counts. `leiden_tfidf` "
            "overwrites .X in place with log1p(TF * IDF * scale_factor), so "
            "an h5ad written before the raw-count stash landed has no counts "
            f"in .X (layers present: {available or 'none'}). PDUI and the "
            "differential tests are only meaningful on integer counts. Fix: "
            "backfill layers['counts'] from the sibling "
            "06_preprocessing/<ds>/preprocessed.h5ad, or re-run clustering "
            "with a build that stashes it. To proceed anyway (NOT a "
            "publishable number) pass --allow-non-count-matrix."
        )
        if not allow_non_count_matrix:
            raise ValueError(msg)
        log.warning("%s (continuing: --allow-non-count-matrix)", msg)
    else:
        log.info("count matrix source: %s (integral)", source)
    return matrix, source


def build_count_dfs(
    adata: ad.AnnData,
    which: str = "both",
    layer: str | None = None,
    allow_non_count_matrix: bool = False,
):
    """Return ``(pdui_df, diff_df, cell_index, pas_index)``.

    The matrix is chosen by :func:`_select_count_matrix` -- ``layers['counts']``
    in preference to ``.X`` -- and must be integral; see that function for the
    *layer* / *allow_non_count_matrix* contract.

    ``which`` selects which dense count frame(s) to materialize, so a caller
    that needs only one does not pay for a second full dense copy of the
    matrix (the switch tests densify the count matrix via ``.toarray()``; for
    the largest cell types that copy is tens of GB, and holding both frames
    at once is what pushed parallel ``ema switch`` tasks into OOM):

    - ``"both"`` (default): both frames, byte-identical to the original
      behaviour -- ``pdui_df`` (PAS x cells) and ``diff_df`` (cells x PAS),
      both derived from the same densified matrix.
    - ``"pdui"``: only ``pdui_df`` (used by :func:`run_length`); ``diff_df``
      is ``None``.  Densifies the transpose directly, never holding the
      cells x PAS copy.
    - ``"diff"``: only ``diff_df`` (used by :func:`run_diff`); ``pdui_df`` is
      ``None``.

    ``cell_index`` / ``pas_index`` are always returned.  Whichever frame(s)
    are built are numerically identical to the ``"both"`` path.
    """
    if which not in ("both", "pdui", "diff"):
        raise ValueError(f"which must be 'both'|'pdui'|'diff', got {which!r}")
    try:
        pas_index = [int(v) for v in adata.var_names]
    except (TypeError, ValueError):
        pas_index = list(adata.var_names)
    cell_index = list(adata.obs_names)
    counts, _source = _select_count_matrix(
        adata, layer=layer, allow_non_count_matrix=allow_non_count_matrix
    )
    pdui_df = None
    diff_df = None
    if which in ("both", "diff"):
        X = counts.toarray() if sp.issparse(counts) else counts
        diff_df = pd.DataFrame(X, index=cell_index, columns=pas_index)
    if which == "both":
        # Reuse the already-densified X so the default path is unchanged.
        pdui_df = pd.DataFrame(X.T, index=pas_index, columns=cell_index)
    elif which == "pdui":
        Xt = counts.T.toarray() if sp.issparse(counts) else counts.T
        pdui_df = pd.DataFrame(Xt, index=pas_index, columns=cell_index)
    return pdui_df, diff_df, cell_index, pas_index


def _build_gene_id_map(adata: ad.AnnData) -> dict[int, str]:
    """Return ``{pas_id: gene_id}`` from ``adata.var['gene_id']``.

    PAS whose ``var_names`` entry does not parse as ``int``, or whose
    ``gene_id`` is empty/NaN, are skipped.  Insertion order follows
    ``adata.var`` so any downstream tiebreak stays deterministic.

    Args:
        adata: Loaded AnnData with a PAS->gene assignment on ``adata.var``.

    Returns:
        Dict keyed by integer PAS id; empty when ``adata.var`` has no
        ``gene_id`` column.
    """
    out: dict[int, str] = {}
    var = getattr(adata, "var", None)
    if var is None or "gene_id" not in getattr(var, "columns", ()):
        log.warning(
            "adata.var has no 'gene_id' column (older h5ad?); the gene-level "
            "PAS map is empty and gene-scoped output will be too."
        )
        return out
    for pas_id_str, gene_id in var["gene_id"].items():
        if not gene_id or (isinstance(gene_id, float) and pd.isna(gene_id)):
            continue
        try:
            pas_id_int = int(pas_id_str)
        except (TypeError, ValueError):
            continue
        out[pas_id_int] = str(gene_id)
    return out


def _build_gene_fallback_map(
    adata: ad.AnnData,
    pas_coords: pd.DataFrame | None = None,
) -> dict[int, list[tuple[str, str, int, int, int]]]:
    """Synthesize a PAS -> gene map carrying a STRAND-AWARE proximal->distal rank.

    Each PAS maps to a single ``(gene_id, "_gene_", transcript_pos=0, rank,
    total)`` entry, where ``rank`` is 1 for the gene's proximal PAS and
    ``total`` for its distal one, computed from genomic coordinates by
    :func:`~ema.quantification.pas_to_isoform.rank_pas_by_genomic_position`
    (plus strand: distal = HIGHER coordinate; minus strand: distal = LOWER
    coordinate).  Used by:

    - the ``per_gene`` aggregation branch of :func:`run_length`, and
    - the ``per_isoform`` branch when ``utr_unmatched="gene"``, to keep PAS
      that overlap no annotated UTR (UTR-agnostic gene-level fallback)
      instead of silently dropping them.

    This REPLACES the ``(gene_id, "_gene_", 0, 1, 1)`` sentinel.  With every
    PAS carrying ``rank=1``, ``ClassicPDUIStrategy``'s stable sort fell back
    to insertion order -- which is genomic-ascending -- so the LOWEST
    coordinate was reported as proximal for minus-strand genes as well: all
    1667 ``_gene_``-fallback rows of the shipped 3'UTR table were inverted,
    turning every shortening call into a lengthening call.

    Args:
        adata: Loaded AnnData with a PAS->gene assignment on ``adata.var``.
        pas_coords: DataFrame indexed by PAS id (as ``str``) with ``start``
            and ``strand`` columns -- normally read from ``pasbed.bed``.
            When ``None``, the same two columns are taken off ``adata.var``
            if it carries them.

    Returns:
        ``{pas_id: [(gene_id, "_gene_", 0, rank, total)]}``.  Empty dict when
        ``adata.var`` has no ``gene_id`` column.

    Raises:
        ValueError: no coordinate/strand source is available, or none of the
            PAS have a coordinate.  Ranking without strand silently applies a
            plus-strand convention, so guessing is not an option.
    """
    from ema.quantification.pas_to_isoform import rank_pas_by_genomic_position

    gene_of_pas = _build_gene_id_map(adata)
    if not gene_of_pas:
        return {}

    coords = pas_coords
    if coords is None:
        var = getattr(adata, "var", None)
        var_cols = set(getattr(var, "columns", ()))
        if {"start", "strand"} <= var_cols:
            coords = var[["start", "strand"]].copy()
            coords.index = var.index.astype(str)
    if coords is None:
        raise ValueError(
            "cannot build the gene-level PAS map: no PAS coordinates or "
            "strands are available (neither a pasbed.bed nor 'start'/"
            "'strand' columns on adata.var). Proximal vs distal is defined "
            "by strand, so ranking without it silently applies a plus-strand "
            "convention and inverts every minus-strand gene. Pass --pasbed."
        )

    pas_gene_df = pd.DataFrame(
        {
            "pas_id": list(gene_of_pas.keys()),
            "gene_id": list(gene_of_pas.values()),
        }
    )
    lookup = coords.reindex(pas_gene_df["pas_id"].astype(str))
    pas_gene_df["start"] = lookup["start"].to_numpy()
    pas_gene_df["strand"] = lookup["strand"].to_numpy()

    n_total = len(pas_gene_df)
    n_missing = int(pd.isna(pas_gene_df["start"]).sum())
    if n_missing == n_total:
        raise ValueError(
            f"cannot build the gene-level PAS map: none of the {n_total} PAS "
            "in this h5ad were found in the PAS coordinate table. The h5ad "
            "and the pasbed.bed almost certainly come from different runs -- "
            "pass the matching --pasbed rather than ranking by input order, "
            "which inverts every minus-strand gene."
        )
    if n_missing:
        log.warning(
            "%d/%d PAS had no coordinate for rank ordering; they are ranked "
            "after the coordinate-known PAS of their gene (input order as "
            "tiebreak) and their proximal/distal assignment is NOT "
            "strand-verified.",
            n_missing, n_total,
        )

    ranked = rank_pas_by_genomic_position(pas_gene_df)
    return {
        int(row.pas_id): [(row.gene_id, "_gene_", 0, int(row.rank), int(row.total))]
        for row in ranked.itertuples(index=False)
    }


def _assert_pdui_strand_convention(df: pd.DataFrame, source: str) -> None:
    """Fail loudly when a PDUI table violates the proximal/distal convention.

    For every row carrying both coordinates and a single strand::

        strand '+'  =>  proximal_start < distal_start
        strand '-'  =>  proximal_start > distal_start

    This one comparison is what exposed the whole strand-inversion family on
    real output (2299/2299 minus-strand genes in the shipped ``per_gene``
    table; 1667/1667 ``_gene_``-fallback rows).  It is vectorised and costs
    one pass per written table, so it runs on every write instead of waiting
    for the next review.

    Args:
        df: an augmented PDUI table (post ``_augment_pdui_df``).
        source: label used in the error message (usually the output path).

    Raises:
        RuntimeError: at least one row is inverted or has a proximal/distal
            strand disagreement.
    """
    needed = {"proximal_start", "distal_start", "proximal_strand", "distal_strand"}
    if df is None or df.empty or not needed.issubset(df.columns):
        return
    prox = pd.to_numeric(df["proximal_start"], errors="coerce")
    dist = pd.to_numeric(df["distal_start"], errors="coerce")
    p_strand = df["proximal_strand"].astype(str)
    d_strand = df["distal_strand"].astype(str)

    have_coords = prox.notna() & dist.notna()
    stranded = have_coords & p_strand.isin(("+", "-")) & d_strand.isin(("+", "-"))
    mismatch = stranded & (p_strand != d_strand)
    same = stranded & (p_strand == d_strand)
    inverted = (
        (same & (p_strand == "+") & ~(prox < dist))
        | (same & (p_strand == "-") & ~(prox > dist))
    )
    bad = mismatch | inverted
    n_bad = int(bad.sum())
    if not n_bad:
        return

    cols = [c for c in ("gene_id", "transcript_id", "proximal_pas_id",
                        "distal_pas_id", "proximal_start", "distal_start",
                        "proximal_strand", "distal_strand")
            if c in df.columns]
    genes = (
        df.loc[bad, "gene_id"].nunique() if "gene_id" in df.columns else n_bad
    )
    raise RuntimeError(
        f"{source}: strand/proximal-distal convention violated on {n_bad} of "
        f"{int(stranded.sum())} coordinate-bearing rows ({genes} gene(s)); "
        f"{int(mismatch.sum())} of those have proximal and distal on "
        "DIFFERENT strands. Required: strand '+' => proximal_start < "
        "distal_start, strand '-' => proximal_start > distal_start. A "
        "violation means proximal and distal were selected without strand, "
        "which turns every shortening call into a lengthening call. First "
        f"offending rows:\n{df.loc[bad, cols].head(3).to_string(index=False)}"
    )


_LEGACY_DIFF_AGG = {"per_isoform": "within_utr"}


def _build_diff_isoform_groups(
    adata: ad.AnnData,
    diff_df: pd.DataFrame,
    h5ad_path: str,
    gtf: str,
    pasbed_path: Path,
    isoform_agg: str,
    utr_unmatched: str,
    n_jobs: int,
) -> tuple[pd.DataFrame, list[dict]]:
    """Build the UTR-scoped test matrix + group list for ``run_diff``'s
    ``within_utr`` / ``between_utr`` aggregation scopes.

    This is the SINGLE shared grouping switch consumed by :func:`_run_grouped_diff`
    (in turn consumed by all four diff strategies -- fisher, nb_multi,
    nb_pairwise, mwu_percell -- via the normal ``strategy.test()`` interface;
    no per-strategy code is touched). Mirrors ``run_length``'s ``per_isoform``
    PAS -> UTR mapping (same GTF-cache + ``bedtools intersect`` via
    :func:`~ema.quantification.pas_to_isoform.map_pas_to_isoforms`) without
    modifying ``run_length`` itself.

    Scopes:
      - ``within_utr``: the unit stays PAS.  Returned matrix IS ``diff_df``
        unchanged; each group is one 3'UTR isoform, ``bg_cols``/``report_cols``
        = the PAS assigned to it.  A PAS overlapping >=1 UTR gets one group
        per UTR (tested once per UTR -> one output row per UTR).  Under
        ``utr_unmatched="gene"`` a PAS with NO UTR overlap gets a single
        fallback group (``"<gene>::_gene_"``) whose ``bg_cols`` is the PAS's
        WHOLE GENE (so it is genuinely "tested at gene level") but whose
        ``report_cols`` is just itself, so UTR-mapped PAS of that gene are
        not double-reported under the fallback bucket too.  Under
        ``utr_unmatched="drop"`` such PAS get no group at all (omitted).
      - ``between_utr``: the unit changes to the 3'UTR itself.  Returns a
        NEW cells x UTR matrix (each column = the summed counts of its
        member PAS, including the ``utr_unmatched="gene"`` fallback bucket,
        which collapses a gene's orphan PAS into one UTR-agnostic column).
        Groups are keyed by GENE; a gene's ``bg_cols``/``report_cols`` are
        its UTR columns.  Only genes with >= 2 UTR columns produce a group
        (a single UTR has nothing to contrast against).

    Returns:
        ``(test_matrix, groups)``.  ``groups`` is a list of dicts with keys
        ``group_id``, ``gene_id``, ``bg_cols``, ``report_cols``.
    """
    from ema.annotate.gtf2isoform_utr import parse_isoform_utrs
    from ema.quantification.pas_to_isoform import map_pas_to_isoforms

    if utr_unmatched not in ("drop", "gene"):
        raise ValueError(f"utr_unmatched={utr_unmatched!r} invalid; expected drop|gene")
    if isoform_agg not in ("within_utr", "between_utr"):
        raise ValueError(
            f"isoform_agg={isoform_agg!r} not handled by _build_diff_isoform_groups "
            "(expected within_utr or between_utr)"
        )

    # Same GTF cache resolution as run_length's per_isoform branch, kept
    # independent so run_length is never touched by this change.
    _gtf_cache = Path(h5ad_path).resolve().parent
    for _ in range(4):
        if (_gtf_cache / "gtf_cache").exists() or (_gtf_cache / "run_config.json").exists():
            _gtf_cache = _gtf_cache / "gtf_cache"
            break
        _gtf_cache = _gtf_cache.parent
    else:
        _gtf_cache = Path.home() / ".cache" / "peakatail" / "gtf"
    _gtf_cache.mkdir(parents=True, exist_ok=True)

    isoform_utrs = parse_isoform_utrs(Path(gtf), cache_dir=_gtf_cache, n_workers=min(n_jobs, 4))
    pas_isoform_map_raw = map_pas_to_isoforms(pasbed_path, isoform_utrs)

    diff_cols = set(diff_df.columns)
    utr_pas_members: dict[str, list] = {}
    gene_of_utr: dict[str, str] = {}
    for pas_id, entries in pas_isoform_map_raw.items():
        if pas_id not in diff_cols:
            continue
        for gene_id, transcript_id, *_rest in entries:
            utr_id = f"{gene_id}::{transcript_id}"
            utr_pas_members.setdefault(utr_id, []).append(pas_id)
            gene_of_utr[utr_id] = gene_id

    # Gene-level fallback bucketing.  Only the PAS->gene assignment is
    # consumed here (the bucket sums its members; nothing reads a rank), so
    # this uses the rank-free _build_gene_id_map rather than
    # _build_gene_fallback_map -- which would demand pasbed coordinates that
    # `switch diff` does not need for this step.
    orphans_by_gene: dict[str, list] = {}
    for pas_id, gene_id in _build_gene_id_map(adata).items():
        if pas_id not in diff_cols or pas_id in pas_isoform_map_raw:
            continue
        orphans_by_gene.setdefault(gene_id, []).append(pas_id)

    if utr_unmatched == "gene":
        for gene_id, orphan_pas in orphans_by_gene.items():
            utr_id = f"{gene_id}::_gene_"
            utr_pas_members[utr_id] = orphan_pas
            gene_of_utr[utr_id] = gene_id
        if orphans_by_gene:
            log.info(
                "run_diff: isoform_agg=%s, utr_unmatched='gene': %d gene(s) "
                "got a UTR-agnostic fallback bucket for %d PAS with no UTR overlap",
                isoform_agg, len(orphans_by_gene),
                sum(len(v) for v in orphans_by_gene.values()),
            )
    elif orphans_by_gene:
        log.info(
            "run_diff: isoform_agg=%s, utr_unmatched='drop': %d PAS with no "
            "UTR overlap across %d gene(s) are OMITTED",
            isoform_agg, sum(len(v) for v in orphans_by_gene.values()),
            len(orphans_by_gene),
        )

    if isoform_agg == "within_utr":
        gene_to_all_pas: dict[str, list] = {}
        for pas_id, entries in gene_fallback_map.items():
            if pas_id in diff_cols:
                gene_to_all_pas.setdefault(entries[0][0], []).append(pas_id)

        groups: list[dict] = []
        for utr_id, members in utr_pas_members.items():
            gene_id = gene_of_utr[utr_id]
            if utr_id.endswith("::_gene_"):
                # Fallback bucket: background = the WHOLE gene (genuinely
                # "tested at gene level"), but only report the orphan PAS
                # -- UTR-mapped PAS of this gene already get their own row
                # from their own UTR group above.
                bg_cols = gene_to_all_pas.get(gene_id, members)
                report_cols = members
            else:
                bg_cols = members
                report_cols = members
            groups.append({
                "group_id": utr_id, "gene_id": gene_id,
                "bg_cols": bg_cols, "report_cols": report_cols,
            })
        return diff_df, groups

    # between_utr: collapse PAS -> UTR-level counts (sum), then group UTRs
    # by gene (>=2 UTRs required to have anything to contrast).
    utr_ids = sorted(utr_pas_members)
    if utr_ids:
        utr_matrix = pd.DataFrame(
            {utr_id: diff_df[members].sum(axis=1) for utr_id, members in utr_pas_members.items()},
            index=diff_df.index,
        )[utr_ids]
    else:
        utr_matrix = pd.DataFrame(index=diff_df.index)

    genes_to_utrs: dict[str, list[str]] = {}
    for utr_id in utr_ids:
        genes_to_utrs.setdefault(gene_of_utr[utr_id], []).append(utr_id)
    groups = [
        {"group_id": gene_id, "gene_id": gene_id, "bg_cols": utrs, "report_cols": utrs}
        for gene_id, utrs in genes_to_utrs.items()
        if len(utrs) >= 2
    ]
    return utr_matrix, groups


def _run_grouped_diff(
    strategy,
    test_matrix: pd.DataFrame,
    groups: list[dict],
    cluster_labels: pd.Series,
    cluster1: str | None,
    cluster2: str | None,
    min_cells_per_group: int,
    n_jobs: int,
    count_mode: str = "cells",
) -> pd.DataFrame:
    """Run ``strategy.test()`` once per group, scoping each call's count
    matrix to that group's ``bg_cols``.

    This is the SINGLE choke point all four diff strategies (fisher,
    nb_multi, nb_pairwise, mwu_percell) go through for ``within_utr`` /
    ``between_utr`` scoping: restricting the columns a strategy sees
    changes its internal background/denominator (fisher's within-group
    2xK contingency table, the NB strategies' per-cell library-size
    offset/total) from "whole gene" / "whole transcriptome" to "this
    group's members" -- WITHOUT touching any strategy's own code.

    Results are filtered down to each group's ``report_cols`` (so PAS only
    present to widen a fallback group's background aren't spuriously
    reported under that group too), tagged with ``diff_group_id``,
    concatenated across groups, and ``qvalue`` (if present) is RECOMPUTED
    via BH-FDR over the POOLED p-values -- a per-group BH correction over a
    handful of rows would be meaningless.
    """
    collected: list[pd.DataFrame] = []
    for group in groups:
        bg_cols = [c for c in group["bg_cols"] if c in test_matrix.columns]
        if len(bg_cols) < 2:
            continue
        sub = test_matrix[bg_cols]
        pas_group_map = {str(c): group["group_id"] for c in bg_cols}
        kwargs: dict = dict(
            count_matrix=sub,
            cluster_labels=cluster_labels,
            min_cells_per_group=min_cells_per_group,
            n_jobs=n_jobs,
            pas_gene_map=pas_group_map,
            count_mode=count_mode,
        )
        if cluster1 is not None:
            kwargs["cluster1"] = cluster1
            kwargs["cluster2"] = cluster2
        result = strategy.test(**kwargs)
        if result is None or result.empty:
            continue
        report_cols = set(group["report_cols"])
        result = result.loc[result.index.isin(report_cols)]
        if result.empty:
            continue
        result = result.copy()
        result["diff_group_id"] = group["group_id"]
        collected.append(result)

    if not collected:
        return pd.DataFrame()

    out = pd.concat(collected, axis=0)
    if "pvalue" in out.columns and len(out) > 0:
        from scipy.stats import false_discovery_control
        out["qvalue"] = false_discovery_control(out["pvalue"].values, method="bh")
    return out


def _dispatch_pair(
    pair: tuple[str, str],
    *,
    strategy_name: str,
    diff_df: pd.DataFrame,
    cluster_labels: pd.Series,
    n_jobs_inner: int,
    min_cells_per_group: int = 10,
    pas_gene_map: dict[str, str] | None = None,
    count_mode: str = "cells",
) -> tuple[str, str, pd.DataFrame]:
    """Top-level wrapper for ``run_one_pair`` suitable for ``Pool.imap_unordered``.

    ``functools.partial`` cannot carry keyword-only args in all Python versions,
    so this thin wrapper unpacks the positional ``pair`` argument and forwards
    the rest as keyword args to :func:`run_one_pair`.

    This function must remain at module level (not a closure) so that the
    ``spawn`` pool can pickle it reliably.

    Args:
        pair: ``(cluster1_label, cluster2_label)`` tuple from the pair list.
        strategy_name: Registry key for the differential strategy.
        diff_df: Full count matrix (cells × PAS).
        cluster_labels: Cell-to-cluster assignment Series.
        n_jobs_inner: Inner worker budget from ``ResourceManager.split_jobs()``.
        min_cells_per_group: Minimum cells per group for a PAS to enter testing.

    Returns:
        Forwarded ``(c1, c2, result_df)`` from :func:`run_one_pair`.
    """
    c1, c2 = pair
    return run_one_pair(
        strategy_name,
        diff_df,
        cluster_labels,
        c1,
        c2,
        n_jobs_inner=n_jobs_inner,
        min_cells_per_group=min_cells_per_group,
        pas_gene_map=pas_gene_map,
        count_mode=count_mode,
    )


def _resolve_pasbed(h5ad_path: str, explicit: str | None = None) -> Path | None:
    """Resolve the pasbed.bed for a clustered h5ad (B4).

    Resolution order:
      1. ``explicit`` (the user's ``--pasbed``) if it exists.
      2. Walk up from the h5ad's directory (up to 4 levels) looking for a
         ``pasbed.bed`` — the run layout puts it a sibling dir away from the
         clustering h5ad, so a naive ``h5ad.parent/pasbed.bed`` misses it.
    Returns the resolved ``Path`` or ``None`` if nothing is found.
    """
    if explicit:
        p = Path(explicit)
        if p.exists():
            return p
    search_root = Path(h5ad_path).resolve().parent
    for _ in range(4):
        candidate = search_root / "pasbed.bed"
        if candidate.exists():
            return candidate
        search_root = search_root.parent
    return None


def _resolve_pas_gene_map(h5ad_path: str) -> dict[str, str]:
    """Resolve ``{pas_id: gene_id}`` for a clustered h5ad (B7).

    Walks up from the h5ad directory (up to 5 levels) looking for a per-dataset
    ``annotatedpas.bed`` (col4=pas_id, col7=gene_id) or ``pas_gene.tsv``
    (pas_id<TAB>gene_id). Returns an empty dict if nothing is found. This is the
    authoritative gene_id source; ``adata.var['gene_id']`` may be ~59% NaN after
    an across-dataset concat (bug B7), so it must never be trusted alone.
    """
    search_root = Path(h5ad_path).resolve().parent
    for _ in range(5):
        annot = search_root / "annotatedpas.bed"
        if annot.exists():
            mapping: dict[str, str] = {}
            try:
                with open(annot) as fh:
                    for line in fh:
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) >= 7 and parts[6]:
                            mapping[str(parts[3])] = parts[6]
                if mapping:
                    return mapping
            except OSError:
                pass
        pg = search_root / "pas_gene.tsv"
        if pg.exists():
            mapping = {}
            try:
                with open(pg) as fh:
                    for line in fh:
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) >= 2 and parts[0] != "pas_id" and parts[1]:
                            mapping[str(parts[0])] = parts[1]
                if mapping:
                    return mapping
            except OSError:
                pass
        search_root = search_root.parent
    return {}


def _repair_gene_id(gene_id_map, pas_gene_map: dict[str, str]):
    """Fill NaN/empty entries of a ``var['gene_id']`` Series from a pas_gene map.

    Joins on the Series INDEX (var_names = pas_id), which is intact even when the
    concat outer-join NaN'd the gene_id *column* (B7). Returns ``(repaired_series,
    n_repaired)``.
    """
    import pandas as pd  # local import — heavy

    series = gene_id_map.copy()
    if not pas_gene_map:
        return series, 0
    # var['gene_id'] frequently arrives as a pandas Categorical (AnnData stores
    # string var columns that way). Writing a repaired value that is not already
    # a category (e.g. a gene id absent from this cell type, or the "-" strand/
    # unknown placeholder) raises under pandas 2.x. Repair on a plain object
    # series so any replacement can be assigned; downstream only uses the values
    # as grouping/lookup keys, so the dtype change is behaviour-preserving.
    if isinstance(series.dtype, pd.CategoricalDtype):
        series = series.astype(object)
    is_missing = series.isna() | (series.astype(str).str.strip().isin(["", "nan", "None"]))
    n_repaired = 0
    for idx in series.index[is_missing]:
        repl = pas_gene_map.get(str(idx))
        if repl:
            series.at[idx] = repl
            n_repaired += 1
    return series, n_repaired


def run_diff(
    h5ad_paths: list[str],
    pasbed: str | None,
    gtf: str | None,
    output_dir: str,
    cluster_pairs: str | None,
    cluster_key: str,
    marker_top_n: int,
    marker_method: str,
    strategy: str,
    fdr: float,
    threads: int | None,
    per_worker_mb: int,
    min_cells_per_group: int = 10,
    count_mode: str = "cells",
    isoform_agg: str = "per_gene",
    utr_unmatched: str = "gene",
    counts_layer: str | None = None,
    allow_non_count_matrix: bool = False,
    progress_manager=None,
) -> dict[tuple[str, str], pd.DataFrame]:
    """Library-level entry point for differential APA testing.

    Runs differential APA (Fisher / NB regression) across cluster pairs on
    one or more h5ad files.  Body is the post-argparse logic of the legacy
    ``cli()`` diff pathway, promoted so ``ema/cli/switch_diff.py`` can call it
    without going through argparse.

    Args:
        h5ad_paths: One or more paths to clustered AnnData h5ad files.
        pasbed: Optional PAS BED path for isoform-aware PDUI.
        gtf: Optional GTF path for isoform-aware PDUI.  Required when
            ``isoform_agg`` is ``"within_utr"`` or ``"between_utr"``.
        output_dir: Directory to write differential/* and pdui_*.tsv files.
        cluster_pairs: ``'c1,c2;c3,c4'``-style pair string, or None for all.
        cluster_key: adata.obs column with cluster labels.  This can be ANY
            obs column, not just a clustering result -- e.g. ``"stage"`` or
            ``"celltype"`` -- so "differential between stages/cell
            types/anything" is already supported by picking the column here
            and (optionally) narrowing to specific contrasts with
            ``cluster_pairs``.
        marker_top_n: Top-N marker PAS per cluster (0 = disabled).
        marker_method: Marker ranking method (wilcoxon / t-test / logreg).
        strategy: Registered differential APA strategy name.
        fdr: FDR q-value threshold for significance.
        threads: Max parallel workers ceiling (or None for auto).
        per_worker_mb: Estimated peak RAM per parallel worker (MB).
        min_cells_per_group: Minimum cells (with nonzero counts for NB strategies)
            in each cluster for a PAS to enter differential testing. Default 10.
        count_mode: Aggregation unit for the fisher strategy's contingency
            table. ``"cells"`` (default, D4) counts each cell at most once via
            per-cell PAS detection, de-pseudoreplicating the test -- this is the
            FDR-calibrated path (issue #74). ``"reads"`` (legacy, opt-in) sums
            read/UMI counts per group and is anti-conservative because
            within-cell reads are correlated -- treat its q-values as a ranking
            screen only. The default was flipped ``"reads"``->``"cells"`` in
            issue #74 so the out-of-the-box path is calibrated; pass
            ``"reads"`` explicitly only for backward comparison. Ignored by the
            NB strategies (they model per-cell overdispersion directly).
        isoform_agg: Scope of each test's background/denominator --
            ``"per_gene"`` (default; unchanged/byte-identical to legacy
            behaviour): each PAS vs the rest of its gene.  ``"within_utr"``:
            each PAS vs the other PAS sharing its 3'UTR isoform (tandem-UTR
            APA); ``"per_isoform"`` is accepted as a legacy alias.
            ``"between_utr"``: PAS are collapsed to 3'UTR-level counts and
            the test asks whether 3'UTR PREFERENCE differs between groups
            (genes with >=2 UTRs only).  ``"within_utr"``/``"between_utr"``
            require ``gtf``; falls back to ``"per_gene"`` with a warning
            when it is missing/unresolvable.
        utr_unmatched: How to handle a PAS with no annotated UTR overlap
            under ``isoform_agg in {"within_utr", "between_utr"}``.
            ``"gene"`` (default) keeps it via a gene-level fallback bucket
            (still tested/counted at the gene level); ``"drop"`` omits it.
            Ignored when ``isoform_agg == "per_gene"``.
        counts_layer: Which matrix to test.  ``None`` (default) prefers
            ``layers['counts']`` and falls back to ``.X``; ``"X"`` forces
            ``.X``; any other value names a layer that must exist.  See
            :func:`_select_count_matrix`.
        allow_non_count_matrix: Downgrade the "matrix is not integral" error
            to a warning.  The differential strategies cast to int
            (``fisher.py``'s contingency table), so running them on TF-IDF
            weights makes the effective n a sum of IDF-weighted logs.
        progress_manager: Optional :class:`~ema.progress.ProgressManager`.  When
            supplied, a ``"Cluster-pair testing"`` stage is registered and
            advanced once per pair completed (both serial and parallel paths).

    Returns:
        Dict mapping ``(c1, c2)`` pairs to their result DataFrames (all h5ads
        accumulated). For omnibus strategies the key is ``("omnibus", "")`` and
        the value is the omnibus DataFrame. Empty dict if no pairs were run.
    """
    isoform_agg = _LEGACY_DIFF_AGG.get(isoform_agg, isoform_agg)
    if isoform_agg not in ("per_gene", "within_utr", "between_utr"):
        raise ValueError(
            f"isoform_agg={isoform_agg!r} invalid; expected "
            "per_gene|within_utr|between_utr (or legacy alias per_isoform)"
        )
    if utr_unmatched not in ("drop", "gene"):
        raise ValueError(f"utr_unmatched={utr_unmatched!r} invalid; expected drop|gene")

    reset_resource_manager()
    rm = ResourceManager(
        user_max_threads=threads,
        user_per_worker_mb=per_worker_mb if per_worker_mb != 300 else None,
    )
    import ema.utils as _ema_utils
    _ema_utils._RM_INSTANCE = rm

    n_jobs = rm.get_n_jobs(per_worker_mb=per_worker_mb, stage="run_diff")
    log.info("run_diff: resources: %s", rm.report())
    log.info("run_diff: n_jobs=%d", n_jobs)

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    diff_dir = out_dir / "differential"
    diff_dir.mkdir(exist_ok=True)

    if count_mode not in ("reads", "cells"):
        raise ValueError(
            f"count_mode must be 'reads' or 'cells', got {count_mode!r}"
        )

    diff_strat = get_diff_strategy(strategy)

    # Issue #74: read-level fisher is anti-conservative (pseudoreplication over
    # correlated within-cell reads). The calibrated per-cell mode is now the
    # DEFAULT (count_mode="cells"); this warning only fires when a caller has
    # explicitly opted back into the legacy reads path, so it is an informed
    # ranking-screen choice rather than a silent trap.
    if strategy == "fisher" and count_mode == "reads":
        log.warning(
            "run_diff: fisher count_mode='reads' is NOT FDR-calibrated -- its "
            "q-values are anti-conservative (a permutation null reports q<0.05 "
            "hits in 100% of runs; issue #74). Use as a RANKING SCREEN only. "
            "For calibrated inference use the default count_mode='cells' or "
            "strategy='nb_pairwise'."
        )

    # Accumulate all pair results across h5ads (returned to caller for viz).
    all_pair_results: dict[tuple[str, str], pd.DataFrame] = {}

    for h5ad_path in h5ad_paths:
        log.info("run_diff: loading %s", h5ad_path)
        adata = ad.read_h5ad(h5ad_path)

        # Marker selection
        if marker_top_n > 0:
            log.info("run_diff: selecting top %d markers per cluster (%s)", marker_top_n, marker_method)
            markers = select_marker_pas(
                adata,
                cluster_key=cluster_key,
                top_n_per_cluster=marker_top_n,
                method=marker_method,
            )
            save_markers(markers, out_dir / "markers.tsv")
            log.info(
                "run_diff: %d unique marker PAS selected (union across %d clusters)",
                len(markers), adata.obs[cluster_key].nunique(),
            )
        else:
            markers = None
            log.info("run_diff: marker subsetting disabled — using all PAS")

        _, diff_df_full, _, _ = build_count_dfs(
            adata, which="diff", layer=counts_layer,
            allow_non_count_matrix=allow_non_count_matrix,
        )
        if markers is not None:
            marker_set = set(markers)
            diff_df = restrict_count_matrix(diff_df_full, list(marker_set), axis="cols")
            log.info("run_diff: count matrix restricted to %d PAS", diff_df.shape[1])
        else:
            diff_df = diff_df_full

        cluster_labels = pd.Series(adata.obs[cluster_key].values, index=adata.obs_names)
        unique_clusters = sorted(cluster_labels.unique().astype(str).tolist())

        # within_utr / between_utr: build the UTR-scoped test matrix + group
        # list ONCE per h5ad (shared by both the omnibus and per-pair
        # branches below via the single _run_grouped_diff choke point).
        # Falls back to per_gene (with a warning) when --gtf / pasbed.bed
        # cannot be resolved, mirroring run_length's per_isoform fallback.
        _isoform_agg = isoform_agg
        _isoform_test_matrix: pd.DataFrame | None = None
        _isoform_groups: list[dict] | None = None
        if _isoform_agg != "per_gene":
            _resolved_gtf = gtf if gtf and Path(gtf).exists() else None
            _resolved_pasbed = _resolve_pasbed(h5ad_path, pasbed) if _resolved_gtf else None
            if _resolved_gtf is None:
                log.warning(
                    "run_diff: --isoform-agg=%s requires --gtf; falling back to per_gene",
                    _isoform_agg,
                )
                _isoform_agg = "per_gene"
            elif _resolved_pasbed is None:
                log.warning(
                    "run_diff: --isoform-agg=%s requires a resolvable pasbed.bed "
                    "(pass --pasbed explicitly); falling back to per_gene",
                    _isoform_agg,
                )
                _isoform_agg = "per_gene"
            else:
                _isoform_test_matrix, _isoform_groups = _build_diff_isoform_groups(
                    adata=adata,
                    diff_df=diff_df,
                    h5ad_path=h5ad_path,
                    gtf=_resolved_gtf,
                    pasbed_path=_resolved_pasbed,
                    isoform_agg=_isoform_agg,
                    utr_unmatched=utr_unmatched,
                    n_jobs=n_jobs,
                )
                log.info(
                    "run_diff: isoform_agg=%s: %d group(s) built for differential testing",
                    _isoform_agg, len(_isoform_groups),
                )

        if diff_strat.supports_multi_condition:
            log.info("run_diff: running %s omnibus across %d clusters", strategy, len(unique_clusters))
            if _isoform_agg == "per_gene":
                df = diff_strat.test(
                    count_matrix=diff_df,
                    cluster_labels=cluster_labels,
                    n_jobs=n_jobs,
                    min_cells_per_group=min_cells_per_group,
                )
            else:
                df = _run_grouped_diff(
                    diff_strat, _isoform_test_matrix, _isoform_groups, cluster_labels,
                    None, None, min_cells_per_group, n_jobs, count_mode=count_mode,
                )
            out_path = diff_dir / f"{strategy}_omnibus.tsv"
            df.to_csv(out_path, sep="\t")
            sig = (df["qvalue"] < fdr).sum() if "qvalue" in df.columns else 0
            log.info("run_diff: omnibus: %d significant PAS (q<%s) -> %s", sig, fdr, out_path)
            all_pair_results[("omnibus", "")] = df
        else:
            if cluster_pairs:
                pairs = parse_cell_combinations(cluster_pairs)
            else:
                pairs = list(combinations(unique_clusters, 2))
                if strategy != "fisher":
                    log.warning(
                        "run_diff: --strategy=%s on all %d pairs may be slow. "
                        "Use --cluster-pairs to subset.",
                        strategy, len(pairs),
                    )
            log.info("run_diff: running %s on %d cluster pairs", strategy, len(pairs))

            n_outer, n_inner = rm.split_jobs(n_outer=len(pairs), per_inner_mb=300)
            log.info(
                "run_diff: parallelism: %d pairs -> (%d outer, %d inner per pair)",
                len(pairs), n_outer, n_inner,
            )

            # Register the cluster-pair testing progress stage.  Always show
            # the bar so the user gets uniform feedback — even a 1/1 bar is
            # better than no visible activity.
            _pair_client = None
            if progress_manager is not None:
                _pair_stage = progress_manager.add_stage(
                    "Cluster-pair testing", total=len(pairs)
                )
                _pair_client = progress_manager.client(_pair_stage)

            # Build pas_gene_map for the WITHIN-GENE Fisher framing.  The
            # fisher strategy uses this to group PAS by gene; nb_pairwise /
            # nb_multi accept-and-ignore via **_ignored.  When the h5ad has
            # no gene_id column the map is left as None and Fisher logs a
            # warning + falls back to the (less APA-correct) global path.
            pas_gene_map: dict[str, str] | None = None
            if "gene_id" in adata.var.columns:
                pas_gene_map = {
                    str(k): str(v)
                    for k, v in adata.var["gene_id"].dropna().items()
                    if str(v).strip()
                }
                log.info(
                    "run_diff: pas_gene_map: %d PAS -> %d unique genes "
                    "(within-gene Fisher comparison enabled)",
                    len(pas_gene_map), len(set(pas_gene_map.values())),
                )
            else:
                log.warning(
                    "run_diff: adata.var has no 'gene_id' column; Fisher "
                    "will fall back to cross-gene comparison."
                )

            pair_results: dict[tuple[str, str], pd.DataFrame] = {}

            if _isoform_agg != "per_gene":
                # within_utr / between_utr: run each pair through the shared
                # grouped-testing choke point (_run_grouped_diff) instead of
                # run_one_pair/_dispatch_pair.  Serial across pairs for now
                # (groups are already the inner unit of work); per_gene's
                # Pool-based pair parallelism below is untouched.
                for c1, c2 in pairs:
                    df = _run_grouped_diff(
                        diff_strat, _isoform_test_matrix, _isoform_groups,
                        cluster_labels, c1, c2, min_cells_per_group, n_inner,
                        count_mode=count_mode,
                    )
                    pair_results[(c1, c2)] = df
                    if _pair_client is not None:
                        _pair_client.advance(1)
            elif n_outer <= 1 or len(pairs) <= 1:
                for c1, c2 in pairs:
                    _, _, df = run_one_pair(
                        strategy, diff_df, cluster_labels, c1, c2,
                        n_jobs_inner=n_inner,
                        min_cells_per_group=min_cells_per_group,
                        pas_gene_map=pas_gene_map,
                        count_mode=count_mode,
                    )
                    pair_results[(c1, c2)] = df
                    if _pair_client is not None:
                        _pair_client.advance(1)
            else:
                worker_fn = functools.partial(
                    _dispatch_pair,
                    strategy_name=strategy,
                    diff_df=diff_df,
                    cluster_labels=cluster_labels,
                    n_jobs_inner=n_inner,
                    min_cells_per_group=min_cells_per_group,
                    pas_gene_map=pas_gene_map,
                    count_mode=count_mode,
                )
                ctx = multiprocessing.get_context("spawn")
                with ctx.Pool(n_outer) as pool:
                    for c1, c2, df in pool.imap_unordered(worker_fn, pairs, chunksize=1):
                        pair_results[(c1, c2)] = df
                        if _pair_client is not None:
                            _pair_client.advance(1)

            # Build annotation lookup once per h5ad, shared across all pairs.
            # gene_id from adata.var (index = pas_id as str).
            _gene_id_map: pd.Series | None = None
            if "gene_id" in adata.var.columns:
                _gene_id_map = adata.var["gene_id"].copy()
                # B7: adata.var['gene_id'] can be ~59% NaN after an across-dataset
                # concat (outer join drops it for PAS absent in some datasets),
                # which silently drops most genes from the findings. Repair the
                # NaN entries from the authoritative per-dataset pas_gene mapping,
                # joined on var_names (pas_id), which survives the concat.
                _pas_gene_map = _resolve_pas_gene_map(h5ad_path)
                _gene_id_map, _n_fixed = _repair_gene_id(_gene_id_map, _pas_gene_map)
                if _n_fixed:
                    log.info(
                        "run_diff: repaired %d NaN gene_id entries from pas_gene map",
                        _n_fixed,
                    )
            else:
                # No column at all — rebuild it entirely from the pas_gene map.
                _pas_gene_map = _resolve_pas_gene_map(h5ad_path)
                if _pas_gene_map:
                    _gene_id_map = pd.Series(
                        {str(v): _pas_gene_map.get(str(v), "") for v in adata.var_names}
                    )
                    log.info(
                        "run_diff: adata.var had no 'gene_id'; rebuilt %d entries "
                        "from pas_gene map", len(_gene_id_map),
                    )
                else:
                    log.warning(
                        "run_diff: adata.var has no 'gene_id' column and no "
                        "pas_gene map found; gene_id will be blank in TSVs."
                    )

            # chrom/start/end/strand from pasbed.bed.
            # B4: honour the explicit --pasbed argument first, then resolve via
            # the run layout (walk up from the h5ad). The old code hardcoded
            # ``h5ad.parent/pasbed.bed`` and ignored the passed ``pasbed`` arg —
            # under the current layout the pasbed is a sibling *dir* away, so the
            # coordinate columns came back silently blank.
            _pasbed_cols: pd.DataFrame | None = None
            _pasbed_path = _resolve_pasbed(h5ad_path, pasbed)
            if _pasbed_path is not None and _pasbed_path.exists():
                try:
                    _bed = pd.read_csv(
                        _pasbed_path,
                        sep="\t",
                        header=None,
                        usecols=[0, 1, 2, 3, 5],
                        names=["chrom", "start", "end", "pas_id", "strand"],
                        dtype=str,
                    )
                    _bed["pas_id"] = _bed["pas_id"].astype(str)
                    _pasbed_cols = _bed.set_index("pas_id")[["chrom", "start", "end", "strand"]]
                    log.info(
                        "run_diff: loaded pasbed for coordinate annotation: %s (%d rows)",
                        _pasbed_path, len(_pasbed_cols),
                    )
                except Exception as _e:
                    log.warning("run_diff: failed to parse pasbed.bed at %s: %s", _pasbed_path, _e)
            else:
                log.debug("run_diff: no pasbed.bed found at %s; coords will be absent", _pasbed_path)

            def _augment_diff_df(
                df_raw: pd.DataFrame,
                c1_label: str,
                c2_label: str,
            ) -> pd.DataFrame:
                """Left-join gene_id and coordinates onto a per-pair result DataFrame.

                Augmented column order: pas_id, gene_id, chrom, start, end,
                strand, cluster1, cluster2, <original stat columns>.
                If a lookup is unavailable the corresponding columns are filled
                with empty strings so the TSV structure stays consistent.
                """
                aug = df_raw.copy()
                # Normalise index name so joins work regardless of whether the
                # strategy set it or left it as a positional range index.
                if aug.index.name != "pas_id":
                    aug.index.name = "pas_id"
                aug = aug.reset_index()  # pas_id becomes a regular column
                aug["pas_id"] = aug["pas_id"].astype(str)

                # --- gene_id join ---
                if _gene_id_map is not None:
                    gene_series = _gene_id_map.rename_axis("pas_id").reset_index()
                    gene_series["pas_id"] = gene_series["pas_id"].astype(str)
                    aug = aug.merge(gene_series, on="pas_id", how="left")
                else:
                    aug["gene_id"] = ""

                # --- coordinate join ---
                if _pasbed_cols is not None:
                    coord_df = _pasbed_cols.rename_axis("pas_id").reset_index()
                    aug = aug.merge(coord_df, on="pas_id", how="left")
                    aug[["chrom", "start", "end", "strand"]] = (
                        aug[["chrom", "start", "end", "strand"]].fillna("")
                    )
                else:
                    aug["chrom"] = ""
                    aug["start"] = ""
                    aug["end"] = ""
                    aug["strand"] = ""

                # --- self-describing cluster columns ---
                aug.insert(0, "cluster2", c2_label)
                aug.insert(0, "cluster1", c1_label)

                # Reorder: pas_id, gene_id, chrom, start, end, strand, cluster1, cluster2, <stats>
                stat_cols = [
                    c for c in aug.columns
                    if c not in {"pas_id", "gene_id", "chrom", "start", "end",
                                 "strand", "cluster1", "cluster2"}
                ]
                aug = aug[
                    ["pas_id", "gene_id", "chrom", "start", "end", "strand",
                     "cluster1", "cluster2"] + stat_cols
                ]
                return aug

            total_sig = 0
            for c1, c2 in pairs:
                df = pair_results[(c1, c2)]
                df_out = _augment_diff_df(df, c1, c2)
                out_path = diff_dir / f"{strategy}_{c1}_vs_{c2}.tsv"
                df_out.to_csv(out_path, sep="\t", index=False)
                if "qvalue" in df.columns:
                    total_sig += (df["qvalue"] < fdr).sum()
                # Store the augmented df (includes gene_id, chrom, strand) so
                # downstream consumers such as the viz gene_track auto-top-N
                # ranking can access gene annotations without re-parsing disk.
                all_pair_results[(c1, c2)] = df_out

            log.info(
                "run_diff: %d pairs: %d total significant PAS (q<%s) -> %s",
                len(pairs), total_sig, fdr, diff_dir,
            )

            # E5: emit the normalized findings_long table (FindingRow contract),
            # keyed by pas_uid + canonical_cluster with an explicit direction.
            try:
                from ema.switch_test.long_output import findings_long, write_long_table
                _arm = f"switch_diff:{strategy}"
                _long = findings_long(
                    {(c1, c2): all_pair_results[(c1, c2)] for c1, c2 in pairs},
                    strategy=strategy, arm=_arm, fdr=fdr,
                )
                _written = write_long_table(_long, str(diff_dir / "switch_diff_long"))
                log.info("run_diff: findings_long (%d rows) -> %s", len(_long), _written)
            except Exception as _e:
                log.warning("run_diff: findings_long emit failed: %s", _e)

    log.info("run_diff: done.")
    return all_pair_results


def run_length(
    h5ad_paths: list[str],
    gtf: str | None,
    output_dir: str,
    cluster_pairs: str | None,
    cluster_key: str,
    strategy: str,
    isoform_agg: str,
    isoform_collapse: str,
    threads: int | None,
    pseudocount: float = 0.0,
    utr_unmatched: str = "gene",
    pasbed: str | None = None,
    counts_layer: str | None = None,
    allow_non_count_matrix: bool = False,
    progress_client=None,
) -> tuple[pd.DataFrame | None, "ad.AnnData | None"]:
    """Library-level entry point for 3'UTR length / PDUI quantification.

    Computes per-cluster PDUI scores using the requested strategy on one or
    more h5ad files.  Body is the PDUI pathway of the legacy ``cli()`` function,
    promoted so ``ema/cli/switch_length.py`` can call it without argparse.

    Args:
        h5ad_paths: One or more paths to clustered AnnData h5ad files.
        gtf: Optional GTF path required for isoform-level aggregation.
        output_dir: Directory to write pdui_*.tsv files.
        cluster_pairs: Reserved for future filtering; currently unused.
        cluster_key: adata.obs column with cluster labels.
        strategy: PDUI strategy name (classic / proportion / shannon).
        isoform_agg: Aggregation level — ``"per_gene"`` or ``"per_isoform"``.
            Legacy values ``"gene"`` / ``"isoform"`` are translated to the
            ``per_*`` form so existing scripts/YAMLs do not silently invoke
            the wrong branch.
        isoform_collapse: How to summarise isoforms — ``"none"`` / ``"mean"`` /
            ``"majority"``.
        threads: Max parallel workers ceiling (or None for auto).
        pseudocount: Added to each per-cell count before PDUI / entropy
            computation.  Default 0.0 preserves original behaviour.  Set to
            e.g. 1.0 to eliminate NaN on zero-count cells.
        utr_unmatched: How to handle PAS that overlap no annotated UTR when
            ``isoform_agg == "per_isoform"``.  ``"gene"`` (default) keeps
            them via a gene-level fallback (UTR-agnostic, same synthesis
            used by the ``per_gene`` path).  ``"drop"`` restores the old
            behaviour where such PAS are silently absent from the PDUI
            output.  Ignored when ``isoform_agg == "per_gene"`` (those PAS
            are always gene-level already).
        pasbed: Explicit ``pasbed.bed`` path.  When ``None`` the file is
            resolved by walking up from each h5ad (see
            :func:`_resolve_pasbed`).  It is REQUIRED: PAS strand + start is
            what defines proximal vs distal, and there is no safe default.
        counts_layer: Which matrix to quantify.  ``None`` (default) prefers
            ``layers['counts']`` and falls back to ``.X``; ``"X"`` forces
            ``.X``; any other value names a layer that must exist.  See
            :func:`_select_count_matrix`.
        allow_non_count_matrix: Downgrade the "matrix is not integral" error
            to a warning.  PDUI on normalised values is not a count ratio;
            this exists only for deliberate exploration.
        progress_client: Optional :class:`~ema.progress.ProgressClient`.  When
            supplied, ``advance(1)`` is called after each h5ad is processed so
            the CLI progress bar ticks forward.

    Returns:
        Tuple ``(pdui_df, adata)`` from the last h5ad processed, or
        ``(None, None)`` if nothing was computed.  ``pdui_df`` is the PDUI
        result DataFrame; ``adata`` is the loaded AnnData (for cluster obs).
    """
    # Backward-compat: translate legacy {"gene","isoform"} tokens to the
    # canonical {"per_gene","per_isoform"} the strategy code branches on.
    # Without this, isoform_agg="gene" silently runs the per_isoform path.
    _LEGACY_AGG = {"gene": "per_gene", "isoform": "per_isoform"}
    if isoform_agg in _LEGACY_AGG:
        log.warning(
            "isoform_agg=%r is a legacy alias; using %r. Please update callers.",
            isoform_agg, _LEGACY_AGG[isoform_agg],
        )
        isoform_agg = _LEGACY_AGG[isoform_agg]
    if isoform_agg not in ("per_gene", "per_isoform"):
        raise ValueError(
            f"isoform_agg={isoform_agg!r} invalid; expected per_gene or per_isoform"
        )
    if isoform_collapse not in ("none", "mean", "majority"):
        raise ValueError(
            f"isoform_collapse={isoform_collapse!r} invalid; expected none|mean|majority"
        )
    if utr_unmatched not in ("drop", "gene"):
        raise ValueError(
            f"utr_unmatched={utr_unmatched!r} invalid; expected drop|gene"
        )

    reset_resource_manager()
    rm = ResourceManager(user_max_threads=threads)
    import ema.utils as _ema_utils
    _ema_utils._RM_INSTANCE = rm

    n_jobs = rm.get_n_jobs(per_worker_mb=300, stage="run_length")
    log.info("run_length: resources: %s", rm.report())

    # Loud warning for the documented no-op flag rather than silent acceptance.
    if cluster_pairs is not None and cluster_pairs != "":
        log.warning(
            "ema switch length: --cluster-pairs is currently unused for length "
            "analysis (your value %r will NOT filter the output).",
            cluster_pairs,
        )

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    pdui_methods = [strategy] if strategy and strategy != "none" else []

    # Track the last computed PDUI df + adata for viz (returned to caller).
    _last_pdui_df: pd.DataFrame | None = None
    _last_adata: "ad.AnnData | None" = None

    for h5ad_path in h5ad_paths:
        log.info("run_length: loading %s", h5ad_path)
        adata = ad.read_h5ad(h5ad_path)
        _last_adata = adata

        pdui_df_full, _, _, _ = build_count_dfs(
            adata, which="pdui", layer=counts_layer,
            allow_non_count_matrix=allow_non_count_matrix,
        )

        if not pdui_methods:
            log.info("run_length: PDUI skipped (no strategy)")
            continue

        # Resolve pasbed.bed ONCE per h5ad (--pasbed, else walk up from the
        # h5ad).  The same coordinate table drives all three consumers: the
        # strand-aware proximal->distal rank the PDUI strategies consume, the
        # per_isoform PAS->UTR intersection below, and the coordinate columns
        # written onto every output row.
        #
        # This is a HARD requirement, not a nicety.  Without start+strand the
        # rank degrades to adata.var input order, which is genomic-ascending
        # -- a silent plus-strand convention that swaps proximal and distal
        # for every minus-strand gene (2299/2299 measured on the shipped
        # table; minus strand is 50.1% of all PDUI genes).  A wrong answer
        # that looks right is worse than no answer, so refuse to guess.
        if pasbed is not None and not Path(pasbed).exists():
            raise FileNotFoundError(
                f"ema switch length: --pasbed {pasbed!r} does not exist. "
                "Refusing to fall back to the directory walk-up, which could "
                "rank PAS against a different run's coordinates."
            )
        _pb = _resolve_pasbed(h5ad_path, pasbed)
        if _pb is None:
            raise FileNotFoundError(
                "ema switch length: no pasbed.bed could be resolved for "
                f"{h5ad_path} (searched its directory and 3 parents"
                + (f"; --pasbed {pasbed!r} does not exist" if pasbed else "")
                + "). PAS start + strand are REQUIRED to order proximal->"
                "distal; without them every minus-strand gene is silently "
                "inverted. Pass --pasbed <pasbed.bed>."
            )
        try:
            _pasbed_cols = pd.read_csv(
                _pb, sep="\t", header=None,
                names=["chrom", "start", "end", "pas_id", "score", "strand"],
                usecols=["chrom", "start", "end", "pas_id", "strand"],
                dtype={"pas_id": str, "chrom": str, "strand": str,
                       "start": "Int64", "end": "Int64"},
            ).set_index("pas_id")[["chrom", "start", "end", "strand"]]
        except Exception as _e:
            raise ValueError(
                f"ema switch length: could not read PAS coordinates from {_pb} "
                f"({_e}). A BED6 (chrom/start/end/pas_id/score/strand) is "
                "required; ranking proximal->distal without it inverts every "
                "minus-strand gene."
            ) from _e
        if _pasbed_cols.index.has_duplicates:
            n_dup = int(_pasbed_cols.index.duplicated().sum())
            log.warning(
                "run_length: %s has %d duplicate pas_id rows; keeping the "
                "first occurrence of each.", _pb, n_dup,
            )
            _pasbed_cols = _pasbed_cols[~_pasbed_cols.index.duplicated()]
        log.info("run_length: PAS coordinates: %d rows from %s",
                 len(_pasbed_cols), _pb)

        # Build pas_isoform_map. Per_isoform needs the GTF (transcript-level
        # UTRs); per_gene only needs the PAS->gene assignment which now lives
        # directly on adata.var (written by ema.matrixfilter.preprocessing when
        # gene_ids are forwarded from annotate).
        pas_isoform_map: dict[int, list[tuple[str, str, int, int, int]]] = {}
        if gtf and Path(gtf).exists() and isoform_agg == "per_isoform":
            # Isoform-aware path requires both GTF and a PAS BED with strand info.
            # The pasbed path is derived from the run dir (set by the new
            # auto-routing helper) so the deprecated `emaout/pasbed.bed`
            # global lookup is gone.
            from ema.annotate.gtf2isoform_utr import parse_isoform_utrs
            from ema.quantification.pas_to_isoform import map_pas_to_isoforms
            # Cache the parsed isoform UTR map next to the source h5ad's run
            # dir (or fall back to the user-level ~/.cache).  Without this the
            # human GTF is re-parsed (~1 GB on disk, ~3-5 GB in memory) on
            # every `ema switch length --isoform-agg per_isoform` invocation.
            # Cap workers via the same ResourceManager n_jobs we already
            # computed for this run so we don't spawn one-per-chromosome.
            _gtf_cache = Path(h5ad_path).resolve().parent
            for _ in range(4):
                if (_gtf_cache / "gtf_cache").exists() or (_gtf_cache / "run_config.json").exists():
                    _gtf_cache = _gtf_cache / "gtf_cache"
                    break
                _gtf_cache = _gtf_cache.parent
            else:
                _gtf_cache = Path.home() / ".cache" / "peakatail" / "gtf"
            _gtf_cache.mkdir(parents=True, exist_ok=True)
            isoform_utrs = parse_isoform_utrs(
                Path(gtf),
                cache_dir=_gtf_cache,
                n_workers=min(n_jobs, 4),
            )
            # Same pasbed.bed the coordinate table above came from -- it was
            # resolved (or the run aborted) before this branch, so there is
            # no second walk-up and no way for the two to disagree.
            pas_isoform_map = map_pas_to_isoforms(_pb, isoform_utrs)
            log.info(
                "run_length: isoform map: %d PAS mapped from %s",
                len(pas_isoform_map), _pb,
            )

            # --utr-unmatched=gene: PAS that overlap NO annotated UTR are
            # silently absent from pas_isoform_map (bedtools inner join).
            # Keep them via the same gene-level fallback the per_gene path
            # uses, rather than dropping them.  PAS that DID map to >=1 UTR
            # keep their UTR entries unchanged (multi-UTR PAS still get one
            # entry per overlapping transcript).
            if utr_unmatched == "gene":
                # Strand-aware ranks, from the same pasbed coordinates: the
                # old sentinel gave every fallback PAS rank=1, so the classic
                # strategy picked proximal/distal by insertion order and
                # inverted all 1667 `_gene_` rows of the shipped table.
                gene_fallback_map = _build_gene_fallback_map(
                    adata, pas_coords=_pasbed_cols
                )
                n_fallback_added = 0
                for pas_id, fallback_entry in gene_fallback_map.items():
                    if pas_id not in pas_isoform_map:
                        pas_isoform_map[pas_id] = fallback_entry
                        n_fallback_added += 1
                if n_fallback_added:
                    log.info(
                        "run_length: utr_unmatched='gene': added %d PAS with "
                        "no UTR overlap via gene-level fallback "
                        "(isoform map now %d PAS)",
                        n_fallback_added, len(pas_isoform_map),
                    )
        elif isoform_agg == "per_isoform" and (not gtf or not Path(gtf).exists()):
            log.warning(
                "run_length: --isoform-agg=per_isoform requires --gtf; "
                "falling back to per_gene"
            )
            isoform_agg = "per_gene"

        # Per_gene path: synthesise the pas_isoform_map from adata.var's
        # PAS->gene assignment, ranked proximal(1)->distal(N) from the real
        # genomic coordinates loaded above.  The helper is the SAME one the
        # per_isoform utr_unmatched="gene" fallback uses, so the two paths
        # can no longer disagree about which end of a gene is proximal --
        # which is exactly how the `_gene_` sentinel survived the per_gene
        # fix and kept inverting 1667/1667 fallback rows.
        if isoform_agg == "per_gene" and not pas_isoform_map:
            pas_isoform_map = _build_gene_fallback_map(
                adata, pas_coords=_pasbed_cols
            )
            log.info(
                "run_length: per_gene map built from adata.var['gene_id']: "
                "%d PAS mapped (rank = strand-aware genomic "
                "proximal->distal order, 1-based)",
                len(pas_isoform_map),
            )

        # Use the **threading** joblib backend for PDUI strategies.  Two
        # reasons, both observed on the full-BAM regression run:
        #   1. The default loky backend pickles a copy of ``pdui_df_full``
        #      (a 1051 x 22419 dense DataFrame ~ 190 MB) into each worker.
        #      With n_jobs=4 the parent peaked at ~13 GB RSS and crashed the
        #      host; threads share memory so this stays at ~600 MB.
        #   2. loky workers also reinitialise OpenBLAS thread pools and
        #      deadlocked even with the BLAS env caps in ema/__init__.py.
        # Inner pandas/numpy ops release the GIL, so threading still gets
        # real parallelism for the gene-level loops.  threadpool_limits(1)
        # belt-and-braces caps any BLAS calls inside the threads.
        from contextlib import contextmanager
        try:
            from threadpoolctl import threadpool_limits as _threadpool_limits
        except ImportError:
            @contextmanager
            def _threadpool_limits(limits=1):  # type: ignore[misc]
                yield
        # `_pasbed_cols` (genomic coordinates per PAS) was already loaded
        # above, before the pas_isoform_map was built, so both the per_gene
        # rank synthesis AND the row-decoration below (`_augment_pdui_df`)
        # share the same lookup — no second read of pasbed.bed.

        # Cell -> cluster map from the AnnData (typically "leiden", but
        # honour whatever the user passed via --cluster-key).
        _cluster_series: pd.Series | None = None
        if cluster_key and cluster_key in adata.obs.columns:
            _cluster_series = adata.obs[cluster_key].astype(str)

        def _augment_pdui_df(df_raw: pd.DataFrame) -> pd.DataFrame:
            """Left-join cluster + per-PAS coordinates onto a strategy output.

            Adds a `cluster` column (cell -> cluster lookup) and, for every
            PAS-id-like column present (`pas_id`, `proximal_pas_id`,
            `distal_pas_id`), four coordinate columns: ``<col>_chrom``,
            ``<col>_start``, ``<col>_end``, ``<col>_strand``.  Falls back
            to empty strings when a lookup is unavailable so the TSV
            structure stays consistent across runs.
            """
            aug = df_raw.copy()

            # Cluster column
            if _cluster_series is not None and "cell" in aug.columns:
                aug["cluster"] = aug["cell"].map(_cluster_series).fillna("")
            else:
                aug["cluster"] = ""

            # Coordinate joins per PAS-id column.  Use a join on a temporary
            # str-cast column so we don't disturb the original dtype.
            pas_cols = [c for c in ("pas_id", "proximal_pas_id", "distal_pas_id")
                        if c in aug.columns]
            # `_pasbed_cols` is guaranteed non-None here -- run_length aborts
            # above when no pasbed resolves.  Blank coordinate columns would
            # also silently disarm _assert_pdui_strand_convention, which has
            # nothing to compare when start/strand are empty.
            for col in pas_cols:
                prefix = "" if col == "pas_id" else col.rsplit("_pas_id", 1)[0] + "_"
                key = aug[col].astype(str)
                joined = key.map(_pasbed_cols["chrom"]).rename(f"{prefix}chrom")
                aug[f"{prefix}chrom"] = joined.fillna("")
                aug[f"{prefix}start"] = key.map(_pasbed_cols["start"]).astype("Int64")
                aug[f"{prefix}end"] = key.map(_pasbed_cols["end"]).astype("Int64")
                aug[f"{prefix}strand"] = key.map(_pasbed_cols["strand"]).fillna("")
            return aug

        for method in pdui_methods:
            strat = get_pdui_strategy(method)
            with _threadpool_limits(limits=1):
                with parallel_backend("threading", n_jobs=n_jobs):
                    df = strat.compute(
                        count_matrix=pdui_df_full,
                        pas_isoform_map=pas_isoform_map,
                        aggregation=isoform_agg,
                        isoform_collapse=isoform_collapse,
                        pseudocount=pseudocount,
                    )
            # Filename comes from the strategy class (not f"pdui_{method}"):
            # proportion and shannon produce non-PDUI quantities so prefixing
            # them with "pdui_" was misleading.  See PDUIStrategy.output_filename.
            out_path = out_dir / getattr(strat, "output_filename", f"pdui_{method}.tsv")
            df_out = _augment_pdui_df(df)
            # Cheap permanent guard: never WRITE a table that violates the
            # strand convention. One vectorised comparison; it is the check
            # that produced the 2299/2299 evidence in the first place.
            _assert_pdui_strand_convention(df_out, str(out_path))
            df_out.to_csv(out_path, sep="\t", index=False)
            log.info(
                "run_length: %s: %d rows -> %s (augmented with cluster + coords)",
                method, len(df_out), out_path,
            )
            _last_pdui_df = df

        # Advance the progress bar once per h5ad (all methods for this h5ad done).
        if progress_client is not None:
            try:
                progress_client.advance(1)
            except Exception:
                pass

    log.info("run_length: done.")
    return _last_pdui_df, _last_adata


# NOTE: The old argparse cli() has been removed.
# Use `ema switch diff` / `ema switch length` instead.
