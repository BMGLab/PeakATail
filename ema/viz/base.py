"""VizStrategy abstract base class."""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class VizStrategy(ABC):
    """Render a single plot kind in a single engine.

    Subclasses MUST set:
      - name (e.g. "umap_matplotlib") — unique across registry
      - plot_type (e.g. "umap") — groups same-purpose strategies
      - engine ("matplotlib" | "plotly" | "scanpy")

    And implement render(data, output_basepath) -> list[Path].
    """
    name: str = field(init=False, default="")
    plot_type: str = field(init=False, default="")
    engine: str = field(init=False, default="")

    @abstractmethod
    def render(self, data: Any, output_basepath: Path) -> list[Path]:
        """Render the plot(s) and return paths to all written files."""
        ...
