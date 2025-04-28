import argparse
from ema.switch_test.switchs import apa_switch

def parse_cell_combinations(cell_combinations_str) -> list:
    return [list(map(str, pair.split(','))) for pair in cell_combinations_str.split(';')]

def cli():
    parser = argparse.ArgumentParser(prog="ema_switch")
    parser.add_argument("--cell_combinations", type=str, required=True,
                        help="List of cell combinations. The list must be of size 2*N, where N is the number of pairs. "
                             "Each pair should be a list of two elements representing the cell indices. "
                             "Example: '1,2;3,4;5,6'")
    args = parser.parse_args()
    cell_combinations = parse_cell_combinations(args.cell_combinations)
    print(cell_combinations)
    apa_switch(cell_combinations)