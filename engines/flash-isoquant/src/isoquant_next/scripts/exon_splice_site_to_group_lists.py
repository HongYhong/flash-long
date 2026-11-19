"""Convert exon splice-site per-group counts to per-molecule group-list format."""

from __future__ import annotations

import argparse
import gzip
from typing import TextIO


def _open(path: str, mode: str = "r") -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t", encoding="utf-8")
    return open(path, mode, encoding="utf-8")


def _split(s: str) -> list[str]:
    return [] if s == "NA" else s.split(";")


def _join(items: list[str]) -> str:
    return ";".join(items) if items else "NA"


def _write_feature(out: TextIO, feature_id: str, acc: list[list[str]], emit_read_ids: bool) -> None:
    gf, gl, gr, rf, rl, rr = acc
    fields = [feature_id, str(len(gf)), str(len(gl)), str(len(gr)), _join(gf), _join(gl), _join(gr)]
    if emit_read_ids:
        fields += [_join(rf), _join(rl), _join(rr)]
    out.write("\t".join(fields) + "\n")


def convert(input_path: str, output_path: str) -> None:
    with _open(input_path) as inf, _open(output_path, "w") as out:
        header = inf.readline().rstrip("\n").split("\t")
        emit_read_ids = "read_ids_full" in header
        out_header = ["region_gene_candidate", "n_full", "n_left", "n_right", "groups_full", "groups_left", "groups_right"]
        if emit_read_ids:
            out_header += ["read_ids_full", "read_ids_left", "read_ids_right"]
        out.write("\t".join(out_header) + "\n")

        cur_feature = None
        acc: list[list[str]] = []
        for line in inf:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            feature_id, gid = cols[0], cols[4]
            nf, nl, nr = int(cols[1]), int(cols[2]), int(cols[3])
            if feature_id != cur_feature:
                if cur_feature is not None:
                    _write_feature(out, cur_feature, acc, emit_read_ids)
                cur_feature = feature_id
                acc = [[], [], [], [], [], []]
            acc[0].extend([gid] * nf)
            acc[1].extend([gid] * nl)
            acc[2].extend([gid] * nr)
            if emit_read_ids:
                acc[3].extend(_split(cols[5]))
                acc[4].extend(_split(cols[6]))
                acc[5].extend(_split(cols[7]))
        if cur_feature is not None:
            _write_feature(out, cur_feature, acc, emit_read_ids)


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert exon splice-site per-group counts to per-molecule group lists")
    parser.add_argument("--input", "-i", required=True)
    parser.add_argument("--output", "-o", required=True)
    args = parser.parse_args()
    convert(args.input, args.output)
