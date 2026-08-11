"""`ema switch geneview` — gene-track visualisation for per-cluster PAS usage."""
from __future__ import annotations

import logging

import click

from ema.cli.common import (
    apply_yaml_to_kwargs,
    common_options,
    parse_log_overrides,
    parse_plot_engines,
    resolve_subcommand_output_dir,
)
from ema.progress import ProgressManager

log = logging.getLogger(__name__)


@click.command(name="geneview")
@common_options(output_default="switch_geneview")
@click.option(
    "--diff-tsv",
    "diff_tsv",
    multiple=True,
    type=click.Path(exists=True, dir_okay=False),
    help=(
        "Differential TSV from `ema switch diff` (repeatable). "
        "Used to auto-select top-N genes by volcano score. "
        "Optional when --gene-id is supplied."
    ),
)
@click.option(
    "--pasbed",
    "pasbed",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="pasbed.bed — PAS coordinates (BED6 format).",
)
@click.option(
    "--h5ad",
    "-i",
    "h5ad",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="clusters.h5ad with leiden cluster labels and gene_id in var.",
)
@click.option(
    "--gtf",
    "gtf",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    help="Ensembl GTF for isoform structure (optional).",
)
@click.option(
    "--top-genes",
    "top_genes",
    type=int,
    default=10,
    show_default=True,
    help="Number of top genes to auto-pick from --diff-tsv results.",
)
@click.option(
    "--gene-id",
    "gene_id",
    multiple=True,
    type=str,
    help="Explicit gene ID(s) to render (repeatable). Union-ed with --top-genes list.",
)
@click.option(
    "--cluster-key",
    "cluster_key",
    type=str,
    default="leiden",
    show_default=True,
    help="obs column carrying cluster labels.",
)
@click.option(
    "--color-key",
    "color_key",
    type=str,
    default=None,
    help=(
        "obs column carrying a condition per track (e.g. healthy / primary "
        "tumour / metastasis). Tracks are tinted by it and a legend is drawn, "
        "so the grouping is readable without long y-labels."
    ),
)
@click.option(
    "--subtitle",
    "subtitle",
    type=str,
    default="",
    help=(
        "Line shown under the gene title, e.g. the cell type the panel is "
        "restricted to. Naming it once in the header beats repeating it on "
        "every track label."
    ),
)
@click.option(
    "--pas-distance-table/--no-pas-distance-table",
    "pas_distance_table_opt",
    default=False,
    show_default=True,
    help=(
        "Draw a table of PAS coordinates and the distance between adjacent "
        "PAS beneath each panel, and write the full table next to the figure "
        "as gene_<id>_pas_distances.csv."
    ),
)
@click.option(
    "--isoform-map",
    "isoform_map",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help=(
        "pdui_classic.tsv from `switch length --isoform-agg per_isoform`. Used "
        "to show which transcript's 3'UTR each PAS was assigned to in the "
        "distance table. The assignment is read from that file, not recomputed."
    ),
)
@click.pass_context
def geneview(ctx: click.Context, **kwargs) -> None:
    """Gene-track visualisation: per-cluster PAS coverage and proportions.

    Renders one panel per gene showing:
    - Per-cluster PAS read coverage (depth-normalised)
    - Per-cluster within-gene PAS proportions
    - Gene isoform structure (when --gtf is supplied)

    Genes are selected from the union of --gene-id and the top-N genes
    ranked by volcano score from --diff-tsv results.
    """
    # Merge YAML config first; CLI flags (explicitly set) always win.
    apply_yaml_to_kwargs(ctx, kwargs)

    # Validate: need at least one source of gene identifiers.
    if not kwargs["diff_tsv"] and not kwargs["gene_id"]:
        raise click.UsageError(
            "Provide at least one of --diff-tsv or --gene-id so genes can be selected."
        )

    # Resolve output dir — auto-route inside the run dir when possible.
    user_explicit_output = ctx.get_parameter_source("output").name == "COMMANDLINE"
    out_dir = resolve_subcommand_output_dir(
        kwargs["output"],
        user_explicit=user_explicit_output,
        source_paths=[kwargs["h5ad"]],
        subdir="switch_geneview",
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    from ema.logging_config import setup_logging, teardown_logging

    setup_logging(
        level=(kwargs["log_level"] or "INFO").split(",")[0],
        output_dir=out_dir,
        quiet=kwargs["quiet"],
        verbose_count=kwargs["verbose"],
        no_log_file=kwargs["no_log_file"],
        per_logger_overrides=parse_log_overrides(kwargs["log_level"]),
    )

    try:
        log.info(
            "ema switch geneview: h5ad=%s top_genes=%d cluster_key=%s",
            kwargs["h5ad"],
            kwargs["top_genes"],
            kwargs["cluster_key"],
        )

        import pandas as pd

        # ------------------------------------------------------------------ #
        # Load AnnData                                                         #
        # ------------------------------------------------------------------ #
        import anndata as ad

        adata = ad.read_h5ad(kwargs["h5ad"])
        log.info("Loaded h5ad: %d cells x %d PAS", adata.n_obs, adata.n_vars)

        # ------------------------------------------------------------------ #
        # Load PAS BED                                                         #
        # ------------------------------------------------------------------ #
        pasbed = pd.read_csv(
            kwargs["pasbed"],
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "pas_id", "score", "strand"],
            dtype={"pas_id": str, "chrom": str, "strand": str},
        )
        log.info("Loaded pasbed: %d PAS entries", len(pasbed))

        # ------------------------------------------------------------------ #
        # Build gene list                                                       #
        # ------------------------------------------------------------------ #
        # Start with explicit gene IDs (stripped, case preserved to match var column).
        explicit_genes: list[str] = [g.strip() for g in kwargs["gene_id"] if g.strip()]

        auto_genes: list[str] = []
        if kwargs["diff_tsv"]:
            # Parse each TSV and group into pair_results dict.
            pair_results: dict[tuple[str, str | None], pd.DataFrame] = {}
            for tsv_path in kwargs["diff_tsv"]:
                try:
                    df = pd.read_csv(tsv_path, sep="\t")
                except Exception as exc:
                    log.warning("Could not read diff TSV %s: %s", tsv_path, exc)
                    continue
                # Derive pair key from cluster1/cluster2 columns when present.
                if "cluster1" in df.columns and "cluster2" in df.columns:
                    c1 = str(df["cluster1"].iloc[0]) if not df.empty else "?"
                    c2 = str(df["cluster2"].iloc[0]) if not df.empty else "?"
                    pair_key: tuple[str, str | None] = (c1, c2)
                else:
                    # Fall back to filename as an opaque key.
                    import pathlib
                    stem = pathlib.Path(tsv_path).stem
                    pair_key = (stem, None)

                # Accumulate — if multiple TSVs map to the same pair, concat.
                if pair_key in pair_results:
                    pair_results[pair_key] = pd.concat(
                        [pair_results[pair_key], df], ignore_index=True
                    )
                else:
                    pair_results[pair_key] = df

            if pair_results:
                from ema.viz._gene_track_helpers import rank_top_genes

                auto_genes = rank_top_genes(pair_results, n=kwargs["top_genes"])
                log.info(
                    "Ranked top-%d genes from %d pair(s): %s",
                    kwargs["top_genes"],
                    len(pair_results),
                    auto_genes,
                )

        # Union explicit + auto while preserving order and deduping.
        seen: set[str] = set()
        gene_list: list[str] = []
        for g in explicit_genes + auto_genes:
            if g not in seen:
                seen.add(g)
                gene_list.append(g)

        if not gene_list:
            log.warning(
                "No genes selected (diff TSVs produced no ranked genes and "
                "--gene-id was not supplied). Nothing to render."
            )
            return

        log.info("Will render %d gene(s): %s", len(gene_list), gene_list)

        # ------------------------------------------------------------------ #
        # PAS -> 3'UTR (transcript) assignment, as recorded by per-isoform     #
        # length. Only the id columns are read; the file has millions of       #
        # (cell, gene, transcript) rows and we need the unique triples only.   #
        # ------------------------------------------------------------------ #
        pas_isoforms: dict[int, list[str]] = {}
        if kwargs.get("isoform_map"):
            try:
                im = pd.read_csv(
                    kwargs["isoform_map"], sep="\t",
                    usecols=["transcript_id", "proximal_pas_id", "distal_pas_id"],
                ).drop_duplicates()
                acc: dict[int, set[str]] = {}
                for col in ("proximal_pas_id", "distal_pas_id"):
                    for pas_id, tid in zip(im[col], im["transcript_id"]):
                        if pd.isna(pas_id) or str(tid) == "_gene_":
                            continue
                        acc.setdefault(int(pas_id), set()).add(str(tid))
                pas_isoforms = {k: sorted(v) for k, v in acc.items()}
                log.info(
                    "Loaded 3'UTR assignment for %d PAS from %s",
                    len(pas_isoforms), kwargs["isoform_map"],
                )
                if not pas_isoforms:
                    log.warning(
                        "--isoform-map %s yielded no PAS->transcript pairs "
                        "(is it a per_gene run? transcript_id would be '_gene_')",
                        kwargs["isoform_map"],
                    )
            except Exception as exc:
                log.warning("Could not read --isoform-map %s: %s",
                            kwargs["isoform_map"], exc)

        # ------------------------------------------------------------------ #
        # Resolve plot engines                                                  #
        # ------------------------------------------------------------------ #
        engines = parse_plot_engines(
            kwargs.get("plot_engine", "matplotlib"),
            kwargs.get("no_plots", False),
        )

        # ------------------------------------------------------------------ #
        # Render one panel per gene                                             #
        # ------------------------------------------------------------------ #
        from pathlib import Path

        from ema.viz._gene_track_helpers import (
            build_gene_panel, load_isoforms_for_gene, load_gene_name_from_gtf,
        )
        from ema.viz import render_all

        figs_dir = out_dir / "figures"
        figs_dir.mkdir(parents=True, exist_ok=True)

        rendered_count = 0
        with ProgressManager(disable=kwargs.get("no_progress", False)) as pm:
            if len(gene_list) >= 5:
                _gene_stage = pm.add_stage("Gene rendering", total=len(gene_list))
                _gene_client = pm.client(_gene_stage)
            else:
                _gene_client = None

            for gene_id in gene_list:
                # Optionally load isoform structure + human-readable name.
                isoforms = None
                gene_name = ""
                if kwargs.get("gtf"):
                    try:
                        isoforms = load_isoforms_for_gene(Path(kwargs["gtf"]), gene_id)
                        log.debug(
                            "Loaded %d isoform(s) for %s from GTF", len(isoforms), gene_id
                        )
                    except Exception as exc:
                        log.warning(
                            "load_isoforms_for_gene failed for %s: %s", gene_id, exc
                        )
                    try:
                        gene_name = load_gene_name_from_gtf(
                            Path(kwargs["gtf"]), gene_id
                        )
                    except Exception as exc:
                        log.warning(
                            "load_gene_name_from_gtf failed for %s: %s", gene_id, exc
                        )

                panel = build_gene_panel(
                    gene_id=gene_id,
                    adata=adata,
                    pasbed=pasbed,
                    cluster_key=kwargs["cluster_key"],
                    isoforms=isoforms,
                    gene_name=gene_name,
                    color_key=kwargs.get("color_key"),
                    subtitle=kwargs.get("subtitle", ""),
                    show_distance_table=kwargs.get("pas_distance_table_opt", False),
                    pas_isoforms=pas_isoforms,
                )
                if panel is None:
                    log.warning(
                        "Gene %s has no PAS in the AnnData or pasbed — skipping.", gene_id
                    )
                    if _gene_client is not None:
                        _gene_client.advance(1)
                    continue

                basepath = figs_dir / f"gene_{gene_id}"
                written = render_all("gene_track", panel, basepath, engines=engines)

                # The in-figure table is capped for legibility; the CSV is not.
                if kwargs.get("pas_distance_table_opt", False):
                    from ema.viz._gene_track_helpers import pas_distance_table

                    dist_csv = figs_dir / f"gene_{gene_id}_pas_distances.csv"
                    pas_distance_table(panel).to_csv(dist_csv, index=False)
                    log.info("Gene %s: wrote %s", gene_id, dist_csv.name)
                log.info(
                    "Gene %s: %d figure file(s) written", gene_id, len(written)
                )
                rendered_count += 1
                if _gene_client is not None:
                    _gene_client.advance(1)

        log.info(
            "ema switch geneview complete: %d/%d gene(s) rendered into %s",
            rendered_count,
            len(gene_list),
            figs_dir,
        )

        # ------------------------------------------------------------------ #
        # Write figures manifest                                                #
        # ------------------------------------------------------------------ #
        try:
            from ema.viz._meta import write_figures_index

            write_figures_index(figs_dir, command="ema switch geneview")
        except Exception as exc:
            log.warning("write_figures_index failed: %s", exc)

    finally:
        teardown_logging()
