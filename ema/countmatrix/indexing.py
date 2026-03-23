class BarcodeIndex:
    """Encapsulates barcode-to-column mapping state for one run."""

    def __init__(self):
        self._cb_total = {}
        self._cb_index = 0

    def get_index(self, cb: str) -> int:
        try:
            return self._cb_total[cb]
        except KeyError:
            self._cb_index += 1
            self._cb_total[cb] = self._cb_index
            return self._cb_index

    def reset(self):
        self._cb_total.clear()
        self._cb_index = 0

    @property
    def mapping(self) -> dict:
        return dict(self._cb_total)


# Module-level singleton — maintains backward compatibility
_index = BarcodeIndex()


def indexing(cb: str) -> int:
    """Return column index for barcode cb. Assigns new index on first encounter."""
    return _index.get_index(cb)


def reset_index():
    """Reset the module-level barcode index. Call between independent runs."""
    _index.reset()


def get_mapping() -> dict:
    """Return a copy of the current barcode-to-index mapping."""
    return _index.mapping
