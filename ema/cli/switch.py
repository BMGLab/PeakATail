"""`ema switch` — group for per-cluster differential APA analyses."""
import click

from ema.cli.switch_diff import diff
from ema.cli.switch_geneview import geneview
from ema.cli.switch_length import length
from ema.cli.switch_match import match
from ema.cli.switch_trend import trend


@click.group(name="switch")
def switch() -> None:
    """Per-cluster differential APA analyses."""
    pass


switch.add_command(diff)
switch.add_command(geneview)
switch.add_command(length)
switch.add_command(match)
switch.add_command(trend)
