"""clip_seeded — poly(A)-clip-seeded PAS calling (Stage 1 Phase 2).

The measured motivation (chr19+chr21, PBMC 10k v3, score_tool metrics):

=============================  =======  ======  ======  =======
arm                                  n   P@100   R@100   F1@100
=============================  =======  ======  ======  =======
PeakATail as shipped            12,889  0.1614  0.1962  0.1771
clip clusters, >=1 read         14,889  0.3218  0.3366  0.3291
=============================  =======  ======  ======  =======

~60% of the poly(A) clip evidence lies OUTSIDE every coverage-peak window,
so no post-hoc filter can recover it — PAS candidates must be *seeded* from
the clip sites themselves.

This strategy is a thin hybrid coordinator:

* ``seeds_from_clips = True`` tells the peak-calling loop (all three paths:
  monolithic, 3-stage pipeline, tile workers) to accumulate every
  qualifying clip read into a chromosome-level
  :class:`~ema.countmatrix.polya.ClipSeeder` and to route the emission
  through its two-tier flush:

  - TIER 1 (primary): 25 bp single-linkage clip-site clusters clearing
    ``--polya-min-reads`` distinct molecules, emitted at the read-weighted
    cluster mode ``[mode, mode+1)`` with the raw clip-read count in BED
    column 5 and the cluster's per-CB clip-read counts as the matrix row.
  - TIER 2 (coverage_only): the inner coverage strategy's PAS that overlap
    no tier-1 cluster, emitted with their clip support in column 5 —
    which is 0 by construction, so ``score == 0`` doubles as the tier tag
    and ``--polya-mode filter`` drops exactly this tier.

* the coverage side simply delegates to an existing registered strategy
  (default: ``lambda_gradient``, the production default), so nothing about
  coverage peak detection is re-implemented or silently changed here.
"""

from __future__ import annotations

import inspect

from ema.strategies import register, get_strategy
from ema.strategies.base import PeakFinderStrategy


@register("clip_seeded")
class ClipSeededStrategy(PeakFinderStrategy):
    """Hybrid clip-seeded strategy: clip clusters primary, coverage second tier.

    ``find_pas`` / ``get_cb_dict_for_pas`` / ``valley_threshold`` delegate to
    the inner coverage strategy — those calls produce the TIER-2 candidates.
    The tier-1 clip clusters are produced by the caller-level
    :class:`~ema.countmatrix.polya.ClipSeeder` (peak calling accumulates
    clip sites chromosome-wide; a strategy only ever sees one peak window,
    which is exactly why seeding cannot live inside ``find_pas``).
    """

    #: Consumed by the peak-calling paths to activate chromosome-level clip
    #: accumulation + two-tier emission.
    seeds_from_clips = True

    def __init__(
        self,
        coverage_strategy: str = "lambda_gradient",
        lambda_method: str | None = None,
        max_pas: int | None = None,
        smoothing_window: int | None = None,
        min_prominence: float | None = None,
    ) -> None:
        """Build the inner coverage strategy, forwarding any tunables it accepts.

        Args:
            coverage_strategy: Registered name of the coverage strategy that
                produces the second-tier candidates.
            lambda_method / max_pas / smoothing_window / min_prominence:
                Standard strategy tunables (same names the CLI forwards);
                ``None`` means "leave the inner strategy's default".
        """
        if coverage_strategy == "clip_seeded":
            raise ValueError("clip_seeded cannot wrap itself")
        supplied = {
            "lambda_method": lambda_method,
            "max_pas": max_pas,
            "smoothing_window": smoothing_window,
            "min_prominence": min_prominence,
        }
        from ema.strategies import _REGISTRY
        inner_cls = _REGISTRY.get(coverage_strategy)
        if inner_cls is None:
            raise ValueError(
                f"Unknown coverage strategy {coverage_strategy!r}. "
                f"Available: {sorted(_REGISTRY)}"
            )
        accepted = set(inspect.signature(inner_cls.__init__).parameters)
        kwargs = {
            k: v for k, v in supplied.items()
            if v is not None and k in accepted
        }
        self.coverage_strategy_name = coverage_strategy
        self.coverage = get_strategy(coverage_strategy, **kwargs)

    # ------------------------------------------------------------------
    # PeakFinderStrategy interface — pure delegation (tier 2)
    # ------------------------------------------------------------------

    def find_pas(self, peak):
        return self.coverage.find_pas(peak)

    def get_cb_dict_for_pas(self, peak, pas_1, pas_2):
        return self.coverage.get_cb_dict_for_pas(peak, pas_1, pas_2)

    def valley_threshold(self, peak, user_min_pas_prominence: float) -> float:
        return self.coverage.valley_threshold(peak, user_min_pas_prominence)

    def get_params(self) -> dict:
        inner = {}
        try:
            inner = self.coverage.get_params()
        except Exception:
            pass
        return {
            "strategy": "clip_seeded",
            "coverage_strategy": self.coverage_strategy_name,
            "coverage_params": inner,
            "description": (
                "clip-site clusters are primary PAS candidates; coverage "
                "peaks without clip support are a second tier (score=0)"
            ),
        }
