#!/usr/bin/env python
"""Capture current PeakATail output as golden reference for regression testing.

Usage:
    python scripts/capture_golden.py --bed emaout/posbed.bed --mtx emaout/posmatrix.mtx --output tests/golden/
"""
import argparse
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ema.validation.regression import capture_golden


def main():
    parser = argparse.ArgumentParser(
        description="Capture current PeakATail output as golden reference"
    )
    parser.add_argument("--bed", required=True, help="Path to BED output file")
    parser.add_argument("--mtx", required=True, help="Path to MTX output file")
    parser.add_argument("--output", required=True, help="Directory for golden reference")

    args = parser.parse_args()

    if not os.path.exists(args.bed):
        print(f"Error: BED file not found: {args.bed}")
        sys.exit(1)
    if not os.path.exists(args.mtx):
        print(f"Error: MTX file not found: {args.mtx}")
        sys.exit(1)

    stats_path = capture_golden(args.bed, args.mtx, args.output)
    print(f"\nGolden reference saved to: {args.output}")
    print(f"Stats: {stats_path}")


if __name__ == "__main__":
    main()
