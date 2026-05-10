import click


@click.command(name="wizard")
def wizard() -> None:
    """Interactive setup wizard (TODO Phase 11)."""
    run()


def run() -> int:
    """Module-level entry called by bare `ema`. Returns exit code."""
    click.echo("ema wizard — not yet implemented (Phase 11)")
    return 0
