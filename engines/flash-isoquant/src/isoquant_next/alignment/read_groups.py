"""Read grouping strategies and multi-column table loaders with strict file handle closure."""

from __future__ import annotations

import gzip
import logging
import os
import sys
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pysam

from isoquant_next.barcodes.table_splitter import split_read_table_by_chromosomes
from isoquant_next.core.error_codes import IsoQuantExitCode

logger = logging.getLogger("IsoQuant")


def _open_table_text(path: str):
    return gzip.open(path, "rt", encoding="utf-8") if path.lower().endswith((".gz", ".gzip")) else open(path, "rt", encoding="utf-8")


def parse_barcode2spot_spec(spec: str) -> Tuple[str, int, List[int]]:
    parts = spec.split(":")
    filename = parts[0]
    barcode_col = int(parts[1]) if len(parts) > 1 and parts[1] else 0
    spot_spec = parts[2] if len(parts) > 2 and parts[2] else "1"
    spot_cols = [int(x) for x in spot_spec.split(",")]
    return filename, barcode_col, spot_cols


class AbstractReadGrouper:
    default_group_id = "NA"

    def get_group_id(self, alignment, read_assignment=None, filename: Optional[str] = None):
        raise NotImplementedError()


class DefaultReadGrouper(AbstractReadGrouper):

    def __init__(self) -> None:
        self.read_groups: set[str] = {self.default_group_id}

    def get_group_id(self, alignment, read_assignment=None, filename: Optional[str] = None) -> str:
        return self.default_group_id


class AlignmentTagReadGrouper(AbstractReadGrouper):
    def __init__(self, tag: str = "RG") -> None:
        self.tag = tag
        self.read_groups: set[Optional[str]] = set()

    def get_group_id(self, alignment, read_assignment=None, filename: Optional[str] = None) -> Optional[str]:
        if alignment is None or not alignment.has_tag(self.tag):
            return None
        val = alignment.get_tag(self.tag)
        self.read_groups.add(val)
        return val


class ReadIdSplitReadGrouper(AbstractReadGrouper):
    default_group_id = "NA"

    def __init__(self, delim: str = "_") -> None:
        self.delim = delim
        self.read_groups: set[str] = set()

    def get_group_id(self, alignment, read_assignment=None, filename: Optional[str] = None) -> str:
        rid = getattr(alignment, "query_name", None) if alignment is not None else None
        if not rid or self.delim not in rid:
            return self.default_group_id
        gid = rid.rsplit(self.delim, 1)[-1]
        self.read_groups.add(gid)
        return gid


class FileNameGrouper(AbstractReadGrouper):
    def __init__(self, args=None, sample=None) -> None:
        self.sample_map: Dict[str, str] = {}
        self.read_groups: set[str] = set()
        if sample is not None and hasattr(sample, "file_list") and hasattr(sample, "readable_names_dict"):
            for fl in sample.file_list:
                f = fl[0] if isinstance(fl, (list, tuple)) else fl
                self.sample_map[f] = sample.readable_names_dict.get(f, os.path.basename(f))

    def get_group_id(self, alignment, read_assignment=None, filename: Optional[str] = None) -> str:
        if filename and filename in self.sample_map:
            gid = self.sample_map[filename]
        elif filename:
            gid = os.path.basename(filename)
        else:
            gid = "NA"
        self.read_groups.add(gid)
        return gid


class BarcodeGrouper(AbstractReadGrouper):
    default_group_id = "NA"

    def __init__(self, barcode_tag: str = "CB") -> None:
        self.barcode_tag = barcode_tag
        self.read_groups: set[str] = set()

    def get_group_id(self, alignment, read_assignment=None, filename: Optional[str] = None) -> str:
        if read_assignment is not None and getattr(read_assignment, "barcode", None):
            bc = read_assignment.barcode
        elif alignment is not None and alignment.has_tag(self.barcode_tag):
            bc = alignment.get_tag(self.barcode_tag)
        else:
            bc = self.default_group_id
        self.read_groups.add(bc)
        return bc


class SharedTableData:
    """Single in-memory read/barcode -> multi-column tuple table shared across column groupers."""

    def __init__(
        self,
        table_filename: str,
        read_id_column_index: int = 0,
        group_id_column_indices: Sequence[int] = (1,),
        delim: str = "\t",
    ) -> None:
        self.table_filename = table_filename
        self.read_id_column_index = int(read_id_column_index)
        self.group_id_column_indices = list(group_id_column_indices)
        self.delim = delim
        self.read_map: Dict[str, List[str]] = (
            load_multicolumn_table(table_filename, self.read_id_column_index, self.group_id_column_indices, delim)
            if os.path.exists(table_filename)
            else {}
        )

    @property
    def num_columns(self) -> int:
        return len(self.group_id_column_indices)


class ReadTableGrouper(AbstractReadGrouper):
    default_group_id = "NA"

    def __init__(self, shared_data: SharedTableData, column_index: int = 0) -> None:
        self.shared_data = shared_data
        self.column_index = column_index
        self.read_groups: set[str] = set()

    def get_group_id(self, alignment, read_assignment=None, filename: Optional[str] = None) -> str:
        rid = getattr(alignment, "query_name", None) if alignment is not None else None
        if not rid and read_assignment is not None:
            rid = getattr(read_assignment, "read_id", None)
        vals = self.shared_data.read_map.get(rid or "")
        if not vals or self.column_index >= len(vals):
            return self.default_group_id
        gid = vals[self.column_index]
        self.read_groups.add(gid)
        return gid


class BarcodeSpotGrouper(AbstractReadGrouper):
    default_group_id = "NA"

    def __init__(self, shared_data: SharedTableData, column_index: int = 0) -> None:
        self.shared_data = shared_data
        self.column_index = column_index
        self.read_groups: set[str] = set()

    def get_group_id(self, alignment, read_assignment=None, filename: Optional[str] = None) -> str:
        bc = None
        if read_assignment is not None and getattr(read_assignment, "barcode", None):
            bc = read_assignment.barcode
        elif alignment is not None and alignment.has_tag("CB"):
            bc = alignment.get_tag("CB")
        vals = self.shared_data.read_map.get(bc or "")
        if not vals or self.column_index >= len(vals):
            return self.default_group_id
        gid = vals[self.column_index]
        self.read_groups.add(gid)
        return gid


class MultiReadGrouper:
    def __init__(self, groupers: Sequence[AbstractReadGrouper] | AbstractReadGrouper) -> None:
        if not isinstance(groupers, list):
            groupers = list(groupers) if isinstance(groupers, Iterable) else [groupers]
        self.groupers = groupers
        self.read_groups: List[set] = [set() for _ in groupers]

    def get_group_id(self, alignment, read_assignment=None, filename: Optional[str] = None) -> List:
        group_ids = []
        for i, grouper in enumerate(self.groupers):
            gid = grouper.get_group_id(alignment, read_assignment, filename)
            group_ids.append(gid)
            self.read_groups[i].add(gid)
        return group_ids

    def get_all_groups(self) -> List[List]:
        return [list(g) for g in self.read_groups]


def _spec_value(values: List[str], index: int, default: str) -> str:
    return values[index] if len(values) > index else default


def _require_barcoded_reads(args, sample, grouping_name: str) -> None:
    if (not hasattr(sample, "barcodes_split_reads") or not sample.barcodes_split_reads) and not getattr(
        args, "barcoded_bam", False
    ):
        logger.critical("%s grouping requires barcoded reads (use --barcoded_reads or --barcoded_bam)", grouping_name)
        sys.exit(int(IsoQuantExitCode.MISSING_REQUIRED_OPTION))


def _barcode_spot_groupers(shared_data: SharedTableData):
    if shared_data.num_columns == 1:
        return BarcodeSpotGrouper(shared_data, column_index=0)
    return [BarcodeSpotGrouper(shared_data, column_index=c) for c in range(shared_data.num_columns)]


def _make_barcode_mapping_groupers(mapping_spec, mapping_arg_name: str, grouping_name: str, args, sample):
    if not mapping_spec:
        logger.critical("%s grouping requires --%s", grouping_name, mapping_arg_name)
        sys.exit(int(IsoQuantExitCode.MISSING_REQUIRED_OPTION))
    _require_barcoded_reads(args, sample, grouping_name)
    filename, barcode_col, spot_cols = parse_barcode2spot_spec(mapping_spec)
    shared_data = SharedTableData(filename, read_id_column_index=barcode_col, group_id_column_indices=spot_cols, delim="\t")
    return _barcode_spot_groupers(shared_data)


def _make_file_groupers(values: List[str], spec_string: str, sample, chr_id: str, spec_index: int):
    if len(values) < 2:
        logger.critical("Group specification %s is too short, specify at least a file name", spec_string)
        sys.exit(int(IsoQuantExitCode.INVALID_PARAMETER))
    read_col = int(values[2]) if len(values) > 2 else 0
    delim = values[4] if len(values) > 4 else "\t"
    group_col_spec = values[3] if len(values) > 3 else "1"
    chr_file = f"{sample.read_group_file}_spec{spec_index}_{chr_id}"
    if "," in group_col_spec:
        cols = [int(x) for x in group_col_spec.split(",")]
        shared_data = SharedTableData(chr_file, read_col, cols, delim)
        return [ReadTableGrouper(shared_data, i) for i in range(len(cols))]
    col_idx = int(group_col_spec)
    shared_data = SharedTableData(chr_file, read_col, [col_idx], delim)
    return ReadTableGrouper(shared_data, 0)


def parse_grouping_spec(spec_string: str, args, sample, chr_id: str, spec_index: int = 0):
    values = spec_string.split(":")
    key = values[0]
    if key == "file_name":
        return FileNameGrouper(args, sample)
    if key == "no_auto":
        return None
    if key == "barcode_spot":
        return _make_barcode_mapping_groupers(
            getattr(args, "barcode2spot", None), "barcode2spot", "barcode_spot", args, sample
        )
    if key == "barcode_barcode":
        return _make_barcode_mapping_groupers(
            getattr(args, "barcode2barcode", None), "barcode2barcode", "barcode_barcode", args, sample
        )
    if key == "barcode":
        _require_barcoded_reads(args, sample, "barcode")
        return BarcodeGrouper()
    if key == "tag":
        return AlignmentTagReadGrouper(tag=_spec_value(values, 1, "RG"))
    if key == "read_id":
        return ReadIdSplitReadGrouper(delim=_spec_value(values, 1, "_"))
    if key == "file":
        return _make_file_groupers(values, spec_string, sample, chr_id, spec_index)
    logger.critical("Unsupported read grouping option: %s", key)
    sys.exit(int(IsoQuantExitCode.INVALID_PARAMETER))


def create_read_grouper(args, sample, chr_id: str):
    if not hasattr(args, "read_group") or args.read_group is None:
        return DefaultReadGrouper()
    specs = args.read_group.split(";") if isinstance(args.read_group, str) else args.read_group
    groupers = []
    for spec_index, spec in enumerate(specs):
        s = spec.strip()
        if not s:
            continue
        g = parse_grouping_spec(s, args, sample, chr_id, spec_index)
        if g:
            if isinstance(g, list):
                groupers.extend(g)
            else:
                groupers.append(g)
    if not groupers:
        return DefaultReadGrouper()
    if len(groupers) == 1:
        return groupers[0]
    return MultiReadGrouper(groupers)


def prepare_read_groups(args, sample) -> None:
    if not hasattr(args, "read_group") or args.read_group is None:
        return
    specs = args.read_group.split(";") if isinstance(args.read_group, str) else args.read_group
    bam_files = [x[0] if isinstance(x, (list, tuple)) else x for x in sample.file_list]
    if not bam_files:
        return
    chromosomes: List[str] = []
    seen: set[str] = set()
    for bf in bam_files:
        with pysam.AlignmentFile(bf, "rb") as af:
            for c in af.references:
                if c not in seen:
                    seen.add(c)
                    chromosomes.append(c)

    for spec_index, spec in enumerate(specs):
        s = spec.strip()
        if not s:
            continue
        values = s.split(":")
        if values[0] != "file":
            continue
        table_filename = values[1]
        prefix = f"{sample.read_group_file}_spec{spec_index}"
        out_map = split_read_table_by_chromosomes(table_filename, bam_files[0], prefix, chromosomes, has_header=False)
        # Alias output paths to expected `<prefix>_<chr>` filenames cleanly
        for c, written_path in out_map.items():
            expected_path = f"{prefix}_{c}"
            if written_path != expected_path and os.path.exists(written_path):
                os.replace(written_path, expected_path)


def get_grouping_pool_types(args) -> Dict[int, str]:
    if not hasattr(args, "read_group") or args.read_group is None:
        return {}
    specs = args.read_group.split(";") if isinstance(args.read_group, str) else args.read_group
    pool_types: Dict[int, str] = {}
    gi = 0
    for spec_index, spec in enumerate(specs):
        s = spec.strip()
        if not s:
            continue
        values = s.split(":")
        stype = values[0]
        if stype == "file_name":
            pool_types[gi] = "file_name"
            gi += 1
        elif stype in ("barcode_spot", "barcode_barcode"):
            attr = getattr(args, "barcode2spot" if stype == "barcode_spot" else "barcode2barcode", None)
            if attr:
                _, _, cols = parse_barcode2spot_spec(attr)
                for ci in range(len(cols)):
                    pool_types[gi] = f"{stype}:{ci}"
                    gi += 1
            else:
                pool_types[gi] = stype
                gi += 1
        elif stype == "barcode":
            pool_types[gi] = "barcode"
            gi += 1
        elif stype in ("tag", "read_id"):
            pool_types[gi] = "dynamic"
            gi += 1
        elif stype == "file":
            delim = values[4] if len(values) >= 5 else "\t"
            if len(values) >= 4 and "," in values[3]:
                for cidx in values[3].split(","):
                    pool_types[gi] = f"tsv:{spec_index}:{cidx}:{delim}"
                    gi += 1
            else:
                cidx = values[3] if len(values) > 3 else "1"
                pool_types[gi] = f"tsv:{spec_index}:{cidx}:{delim}"
                gi += 1
    return pool_types


def get_grouping_strategy_names(args) -> List[str]:
    if not hasattr(args, "read_group") or args.read_group is None:
        return ["default"]
    specs = args.read_group.split(";") if isinstance(args.read_group, str) else args.read_group
    names: List[str] = []
    for spec_index, spec in enumerate(specs):
        s = spec.strip()
        if not s:
            continue
        values = s.split(":")
        stype = values[0]
        if stype == "file_name":
            names.append("file_name")
        elif stype in ("barcode_spot", "barcode_barcode"):
            attr = getattr(args, "barcode2spot" if stype == "barcode_spot" else "barcode2barcode", None)
            if attr:
                _, _, cols = parse_barcode2spot_spec(attr)
                if len(cols) == 1:
                    names.append(stype)
                else:
                    names.extend(f"{stype}_col{ci}" for ci in range(len(cols)))
            else:
                names.append(stype)
        elif stype == "barcode":
            names.append("barcode")
        elif stype == "tag":
            names.append(f"tag_{values[1] if len(values) > 1 else 'RG'}")
        elif stype == "read_id":
            delim = values[1] if len(values) > 1 else "_"
            names.append(f"read_id_{delim.replace('/', '_').replace(chr(92), '_')}")
        elif stype == "file":
            if len(values) >= 4 and "," in values[3]:
                names.extend(f"file{spec_index}_col{c}" for c in values[3].split(","))
            else:
                names.append(f"file{spec_index}_col{values[3] if len(values) > 3 else '1'}")
    return names if names else ["default"]


def load_table(table_tsv_file: str, read_id_column_index: int, group_id_column_index: int, delim: str) -> Dict[str, str]:
    min_cols = max(read_id_column_index, group_id_column_index)
    read_map: Dict[str, str] = {}
    with _open_table_text(table_tsv_file) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            cols = s.split(delim)
            if len(cols) <= min_cols:
                continue
            read_map[cols[read_id_column_index]] = cols[group_id_column_index]
    return read_map


def load_multicolumn_table(
    table_tsv_file: str,
    read_id_column_index: int,
    group_id_column_indices: Sequence[int],
    delim: str,
) -> Dict[str, List[str]]:
    min_cols = max(read_id_column_index, max(group_id_column_indices))
    read_map: Dict[str, List[str]] = {}
    with _open_table_text(table_tsv_file) as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            cols = s.split(delim)
            if len(cols) <= min_cols:
                continue
            read_map[cols[read_id_column_index]] = [cols[i] for i in group_id_column_indices]
    return read_map


def load_barcode2barcode_mapping(filename: str, barcode_col: int, spot_cols: Sequence[int]) -> Dict[str, List[str]]:
    return load_multicolumn_table(filename, barcode_col, spot_cols, "\t")
