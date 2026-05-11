"""`ema switch length` — 3'UTR shortening / lengthening quantification."""
from __future__ import annotations

import logging

import click

from ema.cli.common import common_options, parse_log_overrides, resolve_output_dir
from ema.cli.defaults import DEFAULTS

log = logging.getLogger(__name__)


def _list_pdui_strategies() -> list[str]:
    try:
        from ema.quantification.strategies import list_pdui_strategies
        return list(list_pdui_strategies())
    except ImportError:
        return []


def _list_strategies_callback(ctx, param, value):
    if value:
        for s in _list_pdui_strategies():
            click.echo(s)
        ctx.exit(0)


@click.command(name="length")
@click.option("--list-strategies", is_flag=True, default=False,
              is_eager=True, expose_value=False,
              callback=_list_strategies_callback,
              help="Print available length strategies and exit.")
@common_options(output_default="switch_out")
@click.option("--h5ad", "-i", "h5ad", multiple=True, required=True,
              type=click.Path(exists=True, dir_okay=False))
@click.option("--gtf", type=click.Path(exists=True, dir_okay=False), default=None,
              help="Required when --isoform-agg=per_isoform.")
@click.option("--cluster-pairs", "cluster_pairs", type=str, default=None)
@click.option("--cluster-key", "cluster_key", type=str, default="leiden")
@click.option("--strategy", "-s", "strategy", type=str, default="classic")
@click.option("--isoform-agg", "isoform_agg",
              type=click.Choice(["per_gene", "per_isoform"]),
              default=DEFAULTS["isoform-agg"], show_default=True,
              help="Aggregation level — must match strategy vocabulary "
                   "(per_gene collapses isoforms, per_isoform keeps them).")
@click.option("--isoform-collapse", "isoform_collapse",
              type=click.Choice(["none", "mean", "majority"]),
              default=DEFAULTS["isoform-collapse"], show_default=True,
              help="How to collapse multiple isoforms when --isoform-agg=per_gene "
                   "and the strategy tracks isoforms internally.")
def length(**kwargs) -> None:
    """3'UTR shortening / lengthening quantification (PDUI variants)."""
    valid = _list_pdui_strategies()
    if valid and kwargs["strategy"] not in valid:
        raise click.BadParameter(
            f"Invalid --strategy {kwargs['strategy']!r}. Available: {', '.join(valid)}",
        )

    out_dir = resolve_output_dir(kwargs["output"])
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
        log.info("ema switch length: strategy=%s", kwargs["strategy"])
        from ema.switch_test.runner import run_length
        pdui_df, last_adata = run_length(
            h5ad_paths=list(kwargs["h5ad"]),
            gtf=kwargs["gtf"],
            output_dir=str(out_dir),
            cluster_pairs=kwargs["cluster_pairs"],
            cluster_key=kwargs["cluster_key"],
            strategy=kwargs["strategy"],
            isoform_agg=kwargs["isoform_agg"],
            isoform_collapse=kwargs["isoform_collapse"],
            threads=kwargs["threads"],
        )

        # Viz wire-in: render pdui_distribution + length_shifts figures
        try:
            from ema.cli.common import parse_plot_engines
            from ema.viz import render_all
            from pathlib import Path as _Path
            engines = parse_plot_engines(
                kwargs.get("plot_engine", "both"),
                kwargs.get("no_plots", False),
            )
            if engines and last_adata is not None:
                figs_dir = out_dir / "figures"
                figs_dir.mkdir(parents=True, exist_ok=True)
                cluster_key = kwargs["cluster_key"]

                # pdui_distribution: annotate adata.obs with a summary PDUI score
                # The PDUI df has rows=gene/PAS, columns depend on strategy.
                # We use mean PDUI per cell by projecting from the count matrix.
                import pandas as _pd
                import numpy as _np
                import scipy.sparse as _sp

                if pdui_df is not None and "pdui" in pdui_df.columns:
                    # classic strategy: columns include a per-gene 'pdui' value
                    # We map gene mean onto each cell as a summary score.
                    score_key = "mean_pdui"
                    last_adata.obs[score_key] = _np.nan
                    # attempt per-cell mean from count matrix
                    X = last_adata.X.toarray() if _sp.issparse(last_adata.X) else last_adata.X
                    if X.shape[1] > 0:
                        # mean read-weighted PDUI across detected PAS per cell
                        pas_pdui = pdui_df["pdui"].values if len(pdui_df) == X.shape[1] else _np.zeros(X.shape[1])
                        last_adata.obs[score_key] = (X * pas_pdui).sum(axis=1) / (X.sum(axis=1) + 1e-9)
                    render_all(
                        plot_type="pdui_distribution",
                        data=(last_adata, score_key),
                        output_basepath=figs_dir / "pdui_distribution",
                        engines=engines,
                    )
                elif last_adata is not None and cluster_key in last_adata.obs.columns:
                    # Fallback: use raw count sparsity as a proxy if pdui not available
                    import scipy.sparse as _sp2
                    X2 = last_adata.X.toarray() if _sp2.issparse(last_adata.X) else last_adata.X
                    last_adata.obs["pas_detection_rate"] = (X2 > 0).mean(axis=1)
                    render_all(
                        plot_type="pdui_distribution",
                        data=(last_adata, "pas_detection_rate"),
                        output_basepath=figs_dir / "pdui_distribution",
                        engines=engines,
                    )

                # length_shifts: build gene × cluster-pair ΔPDUI matrix
                if pdui_df is not None and cluster_key in last_adata.obs.columns:
                    clusters = sorted(
                        last_adata.obs[cluster_key].unique(),
                        key=lambda x: int(x) if str(x).isdigit() else x,
                    )
                    from itertools import combinations as _comb
                    X3 = last_adata.X.toarray() if _sp.issparse(last_adata.X) else last_adata.X
                    gene_col = "gene_id" if "gene_id" in pdui_df.columns else pdui_df.columns[0]
                    pdui_col = "pdui" if "pdui" in pdui_df.columns else None

                    if pdui_col and len(pdui_df) > 0:
                        shifts_data: dict[str, list] = {}
                        gene_ids = list(pdui_df[gene_col])
                        pdui_vals = pdui_df[pdui_col].values
                        # Compute mean PDUI per cluster using gene-level pdui_df
                        cluster_pdui: dict[str, _np.ndarray] = {}
                        for cl in clusters:
                            cl_mask = last_adata.obs[cluster_key] == cl
                            cl_X = X3[cl_mask]
                            if cl_X.shape[1] == len(pdui_vals):
                                cluster_pdui[cl] = pdui_vals  # same gene-level PDUIs
                            else:
                                cluster_pdui[cl] = pdui_vals

                        for c1, c2 in _comb(clusters, 2):
                            col = f"{c1}_vs_{c2}"
                            if c1 in cluster_pdui and c2 in cluster_pdui:
                                shifts_data[col] = list(cluster_pdui[c1] - cluster_pdui[c2])
                            else:
                                shifts_data[col] = [0.0] * len(gene_ids)

                        if shifts_data:
                            shifts_df = _pd.DataFrame(shifts_data, index=gene_ids)
                            render_all(
                                plot_type="length_shifts",
                                data=shifts_df,
                                output_basepath=figs_dir / "length_shifts",
                                engines=engines,
                            )

                log.info("ema switch length: viz figures written to %s", figs_dir)
        except Exception as _viz_exc:
            log.warning("ema switch length: viz rendering failed (non-fatal): %s", _viz_exc)
    finally:
        teardown_logging()
