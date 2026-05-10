"""`ema parse-gtf` — pre-warm the global GTF cache."""
from __future__ import annotations

import logging
from pathlib import Path

import click

from ema.cli.common import common_options, parse_log_overrides

log = logging.getLogger(__name__)


@click.command(name="parse-gtf")
@common_options(include_output=False)
@click.option("--gtf", "-g", "gtf",
              type=click.Path(exists=True, dir_okay=False),
              default=None,
              help="GTF file to pre-parse (required unless --show-cache).")
@click.option("--cache-dir", "cache_dir", type=click.Path(file_okay=False),
              default=None, help="Override the global ~/.cache/peakatail/gtf/ location.")
@click.option("--force", is_flag=True, default=False,
              help="Re-parse even on cache hit.")
@click.option("--show-cache", is_flag=True, default=False,
              help="List cache entries with sizes/age and exit.")
def parse_gtf(**kwargs) -> None:
    """Pre-warm the GTF cache so subsequent runs hit instantly."""
    from ema.logging_config import setup_logging, teardown_logging
    setup_logging(
        level=(kwargs["log_level"] or "INFO").split(",")[0],
        output_dir=None,
        quiet=kwargs["quiet"],
        verbose_count=kwargs["verbose"],
        no_log_file=True,
        per_logger_overrides=parse_log_overrides(kwargs["log_level"]),
    )

    try:
        if kwargs["show_cache"]:
            from ema.annotate.gtf_cache import global_cache_dir
            cdir = global_cache_dir()
            if not cdir.exists():
                click.echo("(no cache yet)")
                return
            for entry in sorted(cdir.iterdir()):
                if entry.is_dir():
                    size_mb = sum(p.stat().st_size for p in entry.rglob("*")) / 1024 / 1024
                    click.echo(f"  {entry.name}  ({size_mb:.1f} MB)")
            return

        # --gtf is required when not showing the cache.  Surfaced here (rather than
        # via Click's required=True) so --show-cache can be used standalone.
        if not kwargs["gtf"]:
            raise click.UsageError(
                "Missing option '--gtf' / '-g' (required unless --show-cache is used)."
            )

        from ema.annotate.gtf_cache import process_gtf_cached, lookup_global_cache
        if not kwargs["force"]:
            hit = lookup_global_cache(kwargs["gtf"])
            if hit is not None:
                log.info("cache HIT — %s", hit.get("entry_dir", "(global)"))
                return

        log.info("cache MISS — parsing %s", kwargs["gtf"])
        process_gtf_cached(
            gtf_path=kwargs["gtf"],
            output_dir=kwargs["cache_dir"] or str(Path.home() / ".cache" / "peakatail" / "gtf_run"),
            endbed_path=None,
            features_path=None,
        )
        log.info("done.")
    finally:
        teardown_logging()
