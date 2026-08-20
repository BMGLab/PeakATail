from __future__ import annotations

from ema.utils.resource_manager import ResourceManager  # noqa: F401

#: Name of the ``AnnData.layers`` entry that holds RAW integer counts.
#:
#: ``ema.clustering.clustering.clustering()`` stashes the input matrix here
#: BEFORE any strategy's ``normalize()`` runs, because ``leiden_tfidf``
#: overwrites ``.X`` in place with ``log1p(TF * IDF * scale_factor)``.  Every
#: count consumer (``ema switch length`` PDUI, ``ema switch diff``) must read
#: this layer in preference to ``.X``; see
#: :func:`ema.switch_test.runner.build_count_dfs`.
COUNTS_LAYER = "counts"

_RM_INSTANCE: ResourceManager | None = None


def get_resource_manager() -> ResourceManager:
    """Return the process-wide ResourceManager singleton.

    Lazily initialised on first call.  Reads ``--threads`` from parsed CLI
    args if available (via ``ema.config.args``), otherwise falls back to
    auto-detection with no user ceiling.

    Returns:
        The singleton ``ResourceManager`` instance.  Subsequent calls return
        the same object.
    """
    global _RM_INSTANCE
    if _RM_INSTANCE is None:
        user_threads: int | None = None
        try:
            # ema.config is only present after the CLI entry point has run
            # parse_args().  Guard against ImportError, SystemExit (argparse
            # complains when called with unexpected argv such as pytest flags),
            # and any other exception so we always fall back gracefully.
            from ema.config import args as _args  # type: ignore[import]
            user_threads = getattr(_args, "threads", None)
        except (ImportError, SystemExit, Exception):
            user_threads = None
        _RM_INSTANCE = ResourceManager(user_max_threads=user_threads)
    return _RM_INSTANCE


def reset_resource_manager() -> None:
    """Force the singleton to be re-created on the next ``get_resource_manager()`` call.

    Intended for tests that need to inject a fresh instance with different
    configuration (e.g. a mocked ``--threads`` value).
    """
    global _RM_INSTANCE
    _RM_INSTANCE = None
