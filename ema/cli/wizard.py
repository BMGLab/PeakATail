"""Interactive wizard for PeakATail.

Bare `ema` (no subcommand) routes here. Builds a config via questionary
prompts then dispatches into the matching Click subcommand's logic
(reusing the same code path as the explicit CLI invocations).
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import click
import questionary

log = logging.getLogger(__name__)


_MODES = {
    "Run full pipeline (peak calling -> cluster)": "run",
    "Differential APA test (ema switch diff)": "switch_diff",
    "3'UTR shortening / lengthening (ema switch length)": "switch_length",
    "Cross-dataset cluster matching (ema switch match)": "switch_match",
    "Merge BAMs (ema merge)": "merge",
    "Pre-warm GTF cache (ema parse-gtf)": "parse_gtf",
}


@click.command(name="wizard")
def wizard() -> None:
    """Launch the interactive wizard explicitly."""
    import sys
    sys.exit(run())


def run() -> int:
    """Module-level wizard entry. Returns exit code."""
    mode = _ask_mode()
    if mode is None:
        click.echo("Cancelled — no run started.")
        return 1

    handlers = {
        "run": (_ask_run_config, _dispatch_run),
        "switch_diff": (_ask_switch_diff_config, _dispatch_switch_diff),
        "switch_length": (_ask_switch_length_config, _dispatch_switch_length),
        "switch_match": (_ask_switch_match_config, _dispatch_switch_match),
        "merge": (_ask_merge_config, _dispatch_merge),
        "parse_gtf": (_ask_parse_gtf_config, _dispatch_parse_gtf),
    }
    ask_fn, dispatch_fn = handlers[mode]
    cfg = ask_fn()
    if cfg is None:
        click.echo("Cancelled — no run started.")
        return 1

    if _confirm_save_yaml():
        _save_yaml(cfg, Path(".peakatail/last_run.yaml"))

    if not _confirm_proceed(cfg):
        click.echo("Cancelled — no run started.")
        return 1

    return dispatch_fn(cfg)


# ─── prompts ─────────────────────────────────────────────────────────

def _ask_mode() -> Optional[str]:
    answer = questionary.select("What do you want to do?", choices=list(_MODES.keys())).ask()
    if answer is None:
        return None
    # _MODES.get(answer, answer) handles both: real questionary choices (display
    # string → internal key) and test mocks that return the internal key directly.
    return _MODES.get(answer, answer)


def _ask_run_config() -> Optional[dict]:
    n_str = questionary.select(
        "How many datasets?",
        choices=["1", "2", "3", "4", "5+", "(load YAML)"],
    ).ask()
    if n_str is None:
        return None
    if n_str == "(load YAML)":
        path = questionary.path("YAML config path:").ask()
        if not path:
            return None
        from ema.cli.yaml_loader import load_run_yaml
        return load_run_yaml(path)
    n = 5 if n_str == "5+" else int(n_str)
    if n_str == "5+":
        more = questionary.text("How many?", default="5").ask()
        if not more:
            return None
        n = int(more)

    datasets = []
    for i in range(n):
        ds = _ask_one_dataset(i)
        if ds is None:
            return None
        datasets.append(ds)

    gtf = questionary.path("GTF file:", validate=lambda p: Path(p).exists() or "file not found").ask()
    if gtf is None:
        return None
    atlas = questionary.path(
        "Atlas BED (Enter to skip):",
        validate=lambda p: not p or Path(p).exists() or "file not found",
    ).ask()
    threads = _ask_threads()
    output = questionary.text("Output dir prefix:", default="emaout").ask()

    advanced = questionary.confirm("Configure advanced options?", default=False).ask()
    extras: dict = {}
    if advanced:
        extras = _ask_run_advanced()

    cfg = {
        "datasets": datasets,
        "gtf": gtf,
        "output_dir": output,
        **({"atlas": atlas} if atlas else {}),
        **({"threads": threads} if threads else {}),
        **extras,
    }
    return cfg


def _ask_one_dataset(i: int) -> Optional[dict]:
    ds_id = questionary.text(f"Dataset {i+1} id:", default=f"sample{chr(65+i)}").ask()
    if ds_id is None:
        return None
    merge = questionary.select(
        f"merge_strategy for {ds_id}:",
        choices=["none — independent peak calling per BAM",
                 "before — merge BAMs before peak calling (samtools addreplacerg)",
                 "after — merge per-dataset PAS BEDs after peak calling"],
    ).ask()
    if merge is None:
        return None
    merge_short = merge.split(" ", 1)[0]

    bams = []
    while True:
        bam = questionary.path(
            f"BAM #{len(bams)+1} for {ds_id}:",
            validate=lambda p: Path(p).exists() or "file not found",
        ).ask()
        if bam is None:
            return None
        bams.append(bam)
        more = questionary.confirm("Add another BAM for this dataset?", default=False).ask()
        if not more:
            break
    return {"id": ds_id, "merge_strategy": merge_short, "bams": bams}


def _ask_threads() -> Optional[int]:
    answer = questionary.select(
        "Threads:",
        choices=["auto", "4", "8", "16", "custom"],
    ).ask()
    if answer is None or answer == "auto":
        return None
    if answer == "custom":
        custom = questionary.text("How many?", default="8").ask()
        return int(custom) if custom else None
    return int(answer)


def _ask_run_advanced() -> dict:
    """Tiered advanced options (one prompt per group)."""
    extras: dict = {}
    if questionary.confirm("Tile-based peak calling?", default=False).ask():
        extras["tiles"] = True
    # ... (full grouping per spec §5; abbreviated here for plan length)
    return extras


# Stubs for other modes — pattern is the same shape:
def _ask_switch_diff_config(): return _ask_switch_common("diff")
def _ask_switch_length_config(): return _ask_switch_common("length")
def _ask_switch_match_config(): return _ask_switch_common("match")


def _ask_switch_common(kind: str) -> Optional[dict]:
    h5ads = []
    while True:
        p = questionary.path(
            f"clusters.h5ad #{len(h5ads)+1} (Enter to finish):",
            validate=lambda x: not x or Path(x).exists() or "file not found",
        ).ask()
        if not p:
            break
        h5ads.append(p)
    if not h5ads:
        return None
    if kind == "diff":
        from ema.switch_test.strategies import list_diff_strategies
        choices = list(list_diff_strategies())
    elif kind == "length":
        from ema.quantification.strategies import list_pdui_strategies
        choices = list(list_pdui_strategies())
    else:
        from ema.clustering.cross_dataset import list_match_strategies
        choices = list(list_match_strategies())
    strategy = questionary.select(f"{kind} strategy:", choices=choices).ask()
    if strategy is None:
        return None
    output = questionary.text("Output dir:", default=f"switch_{kind}_out").ask()
    return {"_mode": kind, "h5ad": h5ads, "strategy": strategy, "output": output}


def _ask_merge_config() -> Optional[dict]:
    bams = []
    while True:
        p = questionary.path(f"BAM #{len(bams)+1} (Enter to finish):").ask()
        if not p:
            break
        bams.append(p)
    if not bams:
        return None
    output = questionary.text("Output BAM:", default="merged.bam").ask()
    return {"_mode": "merge", "bam_files": bams, "output": output}


def _ask_parse_gtf_config() -> Optional[dict]:
    gtf = questionary.path("GTF file:").ask()
    if not gtf:
        return None
    return {"_mode": "parse_gtf", "gtf": gtf}


def _confirm_save_yaml() -> bool:
    return questionary.confirm("Save this config to .peakatail/last_run.yaml?", default=True).ask() or False


def _save_yaml(cfg: dict, path: Path) -> None:
    import yaml
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(cfg, sort_keys=False))


def _confirm_proceed(cfg: dict) -> bool:
    click.echo("\n--- Review ---")
    for k, v in cfg.items():
        click.echo(f"  {k}: {v}")
    click.echo("")
    return questionary.confirm("Proceed?", default=True).ask() or False


# ─── dispatch helpers (call into Click commands' callbacks) ──────────

def _dispatch_run(cfg: dict) -> int:
    from ema.cli.run import run as run_cmd
    # Re-construct argv-style invocation. Cleanest path: write cfg to a temp
    # YAML and pass --config to the click command.
    import tempfile, yaml
    tf = tempfile.NamedTemporaryFile("w", suffix=".yaml", delete=False)
    yaml.safe_dump(cfg, tf)
    tf.close()
    try:
        from click.testing import CliRunner
        runner = CliRunner()
        result = runner.invoke(run_cmd, ["--config", tf.name], standalone_mode=False)
        return result.exit_code
    finally:
        Path(tf.name).unlink(missing_ok=True)


def _dispatch_switch_diff(cfg: dict) -> int:
    from ema.cli.switch_diff import diff
    args = ["--strategy", cfg["strategy"], "--output", cfg["output"]]
    for h in cfg["h5ad"]:
        args += ["--h5ad", h]
    from click.testing import CliRunner
    return CliRunner().invoke(diff, args, standalone_mode=False).exit_code


def _dispatch_switch_length(cfg: dict) -> int:
    from ema.cli.switch_length import length
    args = ["--strategy", cfg["strategy"], "--output", cfg["output"]]
    for h in cfg["h5ad"]:
        args += ["--h5ad", h]
    from click.testing import CliRunner
    return CliRunner().invoke(length, args, standalone_mode=False).exit_code


def _dispatch_switch_match(cfg: dict) -> int:
    from ema.cli.switch_match import match
    args = ["--strategy", cfg["strategy"], "--output", cfg["output"]]
    for h in cfg["h5ad"]:
        args += ["--h5ad", h]
    from click.testing import CliRunner
    return CliRunner().invoke(match, args, standalone_mode=False).exit_code


def _dispatch_merge(cfg: dict) -> int:
    from ema.cli.merge import merge
    args = ["--output", cfg["output"]]
    for b in cfg["bam_files"]:
        args += ["--bam-files", b]
    from click.testing import CliRunner
    return CliRunner().invoke(merge, args, standalone_mode=False).exit_code


def _dispatch_parse_gtf(cfg: dict) -> int:
    from ema.cli.parse_gtf import parse_gtf
    from click.testing import CliRunner
    return CliRunner().invoke(parse_gtf, ["--gtf", cfg["gtf"]], standalone_mode=False).exit_code
