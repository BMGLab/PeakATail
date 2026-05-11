"""Allow `python -m ema ...` to dispatch to the Click CLI.

Mirrors the behaviour of the installed ``ema`` console-script entry point
(see pyproject.toml / setup.py). Used by the wizard's ``_spawn_ema``
helper so it can run pipeline commands as subprocesses with live
stdout/stderr and honest exit codes.
"""
from __future__ import annotations

import sys

from ema.cli import main


if __name__ == "__main__":
    sys.exit(main())
