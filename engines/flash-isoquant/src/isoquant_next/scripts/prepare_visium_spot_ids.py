"""Generate barcode-to-spot-ID mapping for Visium HD composite barcodes."""

from __future__ import annotations

import argparse
import csv
import sys
from contextlib import ExitStack
from itertools import product
from typing import Optional


def load_coordinate_file(filepath: str) -> list[tuple[str, list[str]]]:
    entries: list[tuple[str, list[str]]] = []
    with open(filepath, newline="", encoding="utf-8") as f:
        for row in csv.reader(f, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 2:
                sys.exit(f"Error: line has fewer than 2 columns in {filepath}: {row}")
            entries.append((row[0], row[1:]))
    return entries


def generate_spot_ids(
    file1: str,
    file2: str,
    output: str,
    prefixes: list[str],
    suffixes: list[str],
    delimiter: str,
    barcode_delimiter: str,
) -> None:
    entries1 = load_coordinate_file(file1)
    entries2 = load_coordinate_file(file2)
    if not entries1:
        sys.exit(f"Error: no entries loaded from {file1}")
    if not entries2:
        sys.exit(f"Error: no entries loaded from {file2}")

    n_cols = len(entries1[0][1])
    n_cols2 = len(entries2[0][1])
    if n_cols != n_cols2:
        sys.exit(f"Error: column count mismatch: {n_cols} in {file1} vs {n_cols2} in {file2}")

    if len(prefixes) == 1 and n_cols > 1:
        prefixes = prefixes * n_cols
    if len(suffixes) == 1 and n_cols > 1:
        suffixes = suffixes * n_cols
    if len(prefixes) != n_cols or len(suffixes) != n_cols:
        sys.exit("Error: prefix/suffix count does not match coordinate columns")

    with ExitStack() as stack:
        out_handle = (
            sys.stdout
            if output == "-"
            else stack.enter_context(open(output, "w", newline="", encoding="utf-8"))
        )
        writer = csv.writer(out_handle, delimiter="\t")
        for (bc1, coords1), (bc2, coords2) in product(entries1, entries2):
            composite_barcode = f"{bc1}{barcode_delimiter}{bc2}"
            spot_ids = [f"{prefixes[i]}{coords1[i]}{delimiter}{coords2[i]}{suffixes[i]}" for i in range(n_cols)]
            writer.writerow([composite_barcode] + spot_ids)


def main(args: Optional[list[str]] = None) -> None:
    parser = argparse.ArgumentParser(description="Generate barcode-to-spot-ID mapping for Visium HD.")
    parser.add_argument("part1")
    parser.add_argument("part2")
    parser.add_argument("-o", "--output", default="-")
    parser.add_argument("--prefix", nargs="+", default=["s_002um_", "s_008um_", "s_016um_"])
    parser.add_argument("--suffix", nargs="+", default=[""])
    parser.add_argument("--delimiter", default="_")
    parser.add_argument("--barcode_delimiter", default="|")
    parsed = parser.parse_args(args)
    generate_spot_ids(
        parsed.part1, parsed.part2, parsed.output, parsed.prefix, parsed.suffix, parsed.delimiter, parsed.barcode_delimiter
    )
