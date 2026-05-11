from ema.strategies.base import PeakFinderStrategy
from ema.strategies import register


@register("original")
class OriginalStrategy(PeakFinderStrategy):
    """Wraps the existing Peak.pasfind() method unchanged.

    This strategy produces IDENTICAL output to the current PeakATail code.
    It exists to allow A/B comparison with new strategies.
    """

    def find_pas(self, peak):
        pas_1, pas_2 = peak.pasfind()
        if pas_1 == 0:
            return []
        return [(pas_1, pas_2)]

    def get_cb_dict_for_pas(self, peak, pas_1, pas_2):
        # Original behavior: return the full cb_dict for the entire peak
        return peak.cb_dict

    def get_params(self):
        return {"strategy": "original", "description": "wraps pasfind() unchanged"}
