"""Enum frequency aggregation helpers with context-managed file I/O."""

from __future__ import annotations

from collections import defaultdict
from enum import Enum
from typing import Optional, Type


class EnumStats:
    """Tracks counts of enum members and supports lossless dump/reload."""

    def __init__(self, inf: Optional[str] = None) -> None:
        self.stats_dict: dict = defaultdict(int)
        if inf is not None:
            self.load(inf)

    def add(self, val, count: int = 1) -> None:
        self.stats_dict[val] += count

    def merge(self, other: "EnumStats") -> None:
        for k, v in other.stats_dict.items():
            self.stats_dict[k] += v

    def dump(self, out_file: str) -> None:
        with open(out_file, "w", encoding="utf-8") as outf:
            for k, v in self.stats_dict.items():
                key_name = f"{k.__class__.__name__}.{k.name}" if isinstance(k, Enum) else str(k)
                outf.write(f"{key_name}\t{v}\n")

    def load(self, in_file: str, enum_cls: Optional[Type[Enum]] = None) -> None:
        with open(in_file, "r", encoding="utf-8") as inf:
            for line in inf:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) != 2:
                    continue
                raw_key, count_str = parts
                if enum_cls is not None:
                    member_name = raw_key.split(".")[-1]
                    key = enum_cls[member_name]
                else:
                    key = raw_key
                self.stats_dict[key] += int(count_str)


def combine_counts(input_data, output: str) -> None:
    """Lazy wrapper delegating multi-sample combination to Polars-backed quantification engine."""
    from isoquant_next.quantification.combine import combine_counts as _combine_counts

    _combine_counts(input_data, output)
