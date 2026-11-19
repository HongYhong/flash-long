"""Replace barcodes with 2um, 8um, and 16um spot IDs in separate TSV outputs."""

from __future__ import annotations

import argparse
import csv


def load_mapping(file_path: str) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with open(file_path, newline="", encoding="utf-8") as f:
        for row in csv.reader(f, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            mapping[row[0]] = row[1]
    return mapping


def load_spot_map(file_path: str) -> dict[str, dict[str, str]]:
    mapping: dict[str, dict[str, str]] = {}
    with open(file_path, newline="", encoding="utf-8") as f:
        for row in csv.reader(f, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            spot2um, spot8um, spot16um = row
            mapping[spot2um[:-2]] = {"8um": spot8um[:-2], "16um": spot16um[:-2]}
    return mapping


def convert_barcode(
    barcode: str, part1_to_y: dict[str, str], part2_to_x: dict[str, str], spot_map: dict[str, dict[str, str]]
) -> tuple[str | None, str | None, str | None]:
    try:
        p1, p2 = barcode.split("|")
    except ValueError:
        return None, None, None
    if p1 not in part1_to_y or p2 not in part2_to_x:
        return None, None, None
    spot2um = f"s_002um_{part2_to_x[p2]}_{part1_to_y[p1]}"
    if spot2um not in spot_map:
        return spot2um, None, None
    return spot2um, spot_map[spot2um]["8um"], spot_map[spot2um]["16um"]


def transform_reads(
    reads_file: str,
    barcode_col: int,
    part1_to_y: dict[str, str],
    part2_to_x: dict[str, str],
    spot_map: dict[str, dict[str, str]],
    output_prefix: str,
) -> None:
    with (
        open(reads_file, newline="", encoding="utf-8") as infile,
        open(f"{output_prefix}_2um.tsv", "w", newline="", encoding="utf-8") as f2um,
        open(f"{output_prefix}_8um.tsv", "w", newline="", encoding="utf-8") as f8um,
        open(f"{output_prefix}_16um.tsv", "w", newline="", encoding="utf-8") as f16um,
    ):
        reader = csv.reader(infile, delimiter="\t")
        w2 = csv.writer(f2um, delimiter="\t")
        w8 = csv.writer(f8um, delimiter="\t")
        w16 = csv.writer(f16um, delimiter="\t")

        header = next(reader)
        w2.writerow(header)
        w8.writerow(header)
        w16.writerow(header)

        for row in reader:
            spot2um, spot8um, spot16um = convert_barcode(row[barcode_col], part1_to_y, part2_to_x, spot_map)
            for writer, val in ((w2, spot2um), (w8, spot8um), (w16, spot16um)):
                out_row = row[:]
                out_row[barcode_col] = val if val else "*"
                writer.writerow(out_row)


def main() -> None:
    parser = argparse.ArgumentParser(description="Replace barcodes with 2um, 8um, 16um spot IDs.")
    parser.add_argument("reads_tsv")
    parser.add_argument("part1_to_y")
    parser.add_argument("part2_to_x")
    parser.add_argument("spot_map")
    parser.add_argument("-c", "--barcode-col", type=int, default=2)
    parser.add_argument("-o", "--output-prefix", default="reads_converted")
    args = parser.parse_args()
    transform_reads(
        args.reads_tsv,
        args.barcode_col - 1,
        load_mapping(args.part1_to_y),
        load_mapping(args.part2_to_x),
        load_spot_map(args.spot_map),
        args.output_prefix,
    )
