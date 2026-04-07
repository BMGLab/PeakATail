"""CLI entry point for the APA switch test."""

import argparse
from ema.switch_test.switchs import apa_switch


def parse_cell_combinations(cell_combinations_str) -> list:
    """Parse cell combination string into list of pairs.

    Format: '0,1;0,2;1,2' -> [['0','1'], ['0','2'], ['1','2']]
    """
    return [
        list(map(str, pair.split(',')))
        for pair in cell_combinations_str.split(';')
    ]


def cli():
    parser = argparse.ArgumentParser(prog="ema_switch")
    parser.add_argument(
        "--cell_combinations", type=str, required=True,
        help="Cluster pairs to compare. Format: '0,1;0,2;1,2'",
    )
    parser.add_argument(
        "--fdr-threshold", dest="fdr_threshold", type=float, default=0.05,
        help="FDR q-value threshold for significance (default: 0.05)",
    )
    parser.add_argument(
        "--diff-method", dest="diff_method", type=str, default="fisher",
        choices=["fisher"],
        help="Differential APA test method (default: fisher)",
    )
    parser.add_argument(
        "--no-pdui", dest="compute_pdui", action="store_false", default=True,
        help="Skip PDUI calculation",
    )
    args = parser.parse_args()

    cell_combinations = parse_cell_combinations(args.cell_combinations)
    apa_switch(
        cell_combinations,
        fdr_threshold=args.fdr_threshold,
        compute_pdui=args.compute_pdui,
    )
