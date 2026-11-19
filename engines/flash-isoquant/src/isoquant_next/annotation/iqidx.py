"""Single-file zero-copy mmap CSR annotation index (.iqidx) eliminating SQLite N+1 queries."""

from __future__ import annotations

import json
import mmap
import os
import sqlite3
import struct
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from isoquant_next.annotation.gtf_validator import (
    detect_annotation_format,
    open_text_maybe_gz,
    parse_gff3_attributes,
    parse_gtf_attributes,
)
from isoquant_next.core.utils import TranscriptNaming, overlaps

IQIDX_MAGIC = b"IQIDX001"
_HEADER_PREFIX = struct.Struct("<8sQ")

DEFAULT_RNA_TYPES = frozenset({
    "transcript", "RNA", "mRNA", "miRNA", "ncRNA", "tRNA", "rRNA",
    "snRNA", "lnc_RNA", "snoRNA", "telomerase_RNA", "Y_RNA", "scRNA",
})
OTHER_FEATURE_CODES: Dict[str, int] = {
    "CDS": 0,
    "UTR": 1,
    "start_codon": 2,
    "stop_codon": 3,
}
OTHER_FEATURE_NAMES: Tuple[str, ...] = ("CDS", "UTR", "start_codon", "stop_codon")

_SKIP_GENE_ATTRS = frozenset({"gene_id", "ID", "level", "Parent"})
_SKIP_TX_ATTRS = frozenset({"transcript_id", "gene_id", "ID", "level", "exons", "Parent"})
_STRAND_TO_I8 = {"+": 1, "-": -1, ".": 0}
_I8_TO_STRAND = {1: "+", -1: "-", 0: "."}


def _format_filtered_attrs(attrs: Dict[str, str], skip_keys: frozenset[str]) -> str:
    parts = [f'{k} "{v}";' for k, v in attrs.items() if k not in skip_keys and v]
    return (" ".join(parts) + " ") if parts else ""


def _resolve_tx_exons_with_cds_fallback(
    exons: List[Tuple[int, int, int]], other_feats: List[Tuple[int, int, int]]
) -> List[Tuple[int, int, int]]:
    """Apply exact GeneInfo.set_introns_and_exons rule: keep exons + non-overlapping CDS if mixed/CDS-only."""
    if not other_feats:
        exons.sort(key=lambda x: (x[0], x[1]))
        return exons
    has_cds = any(t == 0 for _, _, t in other_feats)
    if not has_cds:
        exons.sort(key=lambda x: (x[0], x[1]))
        return exons

    # Merge (is_exon, start, end, exon_id_idx) sorted by start matching gffutils order_by='start'
    combined: List[Tuple[int, int, int, int]] = [(start, end, 0, ex_id) for start, end, ex_id in exons]
    for start, end, fcode in other_feats:
        if fcode == 0:
            combined.append((start, end, 1, -1))
    combined.sort(key=lambda item: (item[0], item[2], item[1]))

    resolved: List[Tuple[int, int, int]] = []
    n = len(combined)
    for idx, (s, e, kind, ex_id) in enumerate(combined):
        if kind == 0:
            resolved.append((s, e, ex_id))
            continue
        prev_c = None if idx == 0 else (combined[idx - 1][0], combined[idx - 1][1])
        next_c = None if idx == n - 1 else (combined[idx + 1][0], combined[idx + 1][1])
        if (prev_c and overlaps((s, e), prev_c)) or (next_c and overlaps((s, e), next_c)):
            continue
        resolved.append((s, e, -1))
    resolved.sort(key=lambda x: (x[0], x[1]))
    return resolved


class _CompactStringTable:
    """Packs deduplicated UTF-8 strings into contiguous bytes + uint64 offsets."""

    def __init__(self) -> None:
        self._map: Dict[str, int] = {"": 0}
        self._buf = bytearray()
        self._offsets: List[int] = [0, 0]

    def intern(self, s: Optional[str]) -> int:
        if not s:
            return 0
        idx = self._map.get(s)
        if idx is not None:
            return idx
        idx = len(self._offsets) - 1
        self._map[s] = idx
        self._buf.extend(s.encode("utf-8"))
        self._offsets.append(len(self._buf))
        return idx

    def finalize(self) -> Tuple[np.ndarray, np.ndarray]:
        offsets = np.asarray(self._offsets, dtype=np.uint64)
        data = np.frombuffer(self._buf, dtype=np.uint8).copy()
        return offsets, data


def _write_iqidx_bundle(
    out_path: str | Path,
    chrom_order: List[str],
    genes_by_chrom: Dict[str, Dict[str, dict]],
    tx_by_id: Dict[str, dict],
) -> Path:
    """Build CSR NumPy arrays per sorted chromosome and atomically write `.iqidx`."""
    str_tab = _CompactStringTable()

    chrom_names = list(chrom_order)
    chrom_gene_ptr: List[int] = [0]

    g_start: List[int] = []
    g_end: List[int] = []
    g_strand: List[int] = []
    g_chrom: List[int] = []
    g_tx_ptr: List[int] = [0]
    g_id_str: List[int] = []
    g_src_str: List[int] = []
    g_attr_str: List[int] = []

    t_start: List[int] = []
    t_end: List[int] = []
    t_strand: List[int] = []
    t_gene: List[int] = []
    t_exon_ptr: List[int] = [0]
    t_other_ptr: List[int] = [0]
    t_id_str: List[int] = []
    t_src_str: List[int] = []
    t_attr_str: List[int] = []

    e_start: List[int] = []
    e_end: List[int] = []
    e_id_str: List[int] = []

    o_start: List[int] = []
    o_end: List[int] = []
    o_type: List[int] = []

    for cid_idx, chrom in enumerate(chrom_names):
        gmap = genes_by_chrom.get(chrom, {})
        # Sort genes strictly by (start, end, gene_id) matching gffutils order_by='start'
        sorted_genes = sorted(gmap.values(), key=lambda g: (g["start"], g["end"], g["id"]))
        for g in sorted_genes:
            tx_ids = g["tx_ids"]
            # Resolve spans and exons for child transcripts sorted by (start, end, id)
            tx_recs = []
            for tid in tx_ids:
                tx = tx_by_id[tid]
                raw_exons = [
                    (s, e, str_tab.intern(ex_id) if ex_id else -1)
                    for (s, e, ex_id) in tx["exons"]
                ]
                raw_others = sorted(tx["others"], key=lambda x: (x[0], x[1], x[2]))
                resolved_exons = _resolve_tx_exons_with_cds_fallback(raw_exons, raw_others)
                if tx["start"] <= 0 and resolved_exons:
                    tx["start"] = min(x[0] for x in resolved_exons)
                    tx["end"] = max(x[1] for x in resolved_exons)
                tx_recs.append((tx, resolved_exons, raw_others))

            tx_recs.sort(key=lambda item: (item[0]["start"], item[0]["end"], item[0]["id"]))
            if g["start"] <= 0 and tx_recs:
                g["start"] = min(item[0]["start"] for item in tx_recs)
                g["end"] = max(item[0]["end"] for item in tx_recs)

            gene_idx = len(g_start)
            g_start.append(g["start"])
            g_end.append(g["end"])
            g_strand.append(_STRAND_TO_I8.get(g["strand"], 0))
            g_chrom.append(cid_idx)
            g_id_str.append(str_tab.intern(g["id"]))
            g_src_str.append(str_tab.intern(g["source"]))
            g_attr_str.append(str_tab.intern(g["attr_str"]))

            for tx, resolved_exons, raw_others in tx_recs:
                t_start.append(tx["start"])
                t_end.append(tx["end"])
                t_strand.append(_STRAND_TO_I8.get(tx["strand"], 0))
                t_gene.append(gene_idx)
                t_id_str.append(str_tab.intern(tx["id"]))
                t_src_str.append(str_tab.intern(tx["source"]))
                t_attr_str.append(str_tab.intern(tx["attr_str"]))

                for es, ee, eid_idx in resolved_exons:
                    e_start.append(es)
                    e_end.append(ee)
                    e_id_str.append(eid_idx)
                t_exon_ptr.append(len(e_start))

                for os_, oe_, ocode in raw_others:
                    o_start.append(os_)
                    o_end.append(oe_)
                    o_type.append(ocode)
                t_other_ptr.append(len(o_start))

            g_tx_ptr.append(len(t_start))
        chrom_gene_ptr.append(len(g_start))

    str_offsets, str_bytes = str_tab.finalize()

    arrays: dict[str, np.ndarray] = {
        "chrom_gene_ptr": np.asarray(chrom_gene_ptr, dtype=np.int64),
        "gene_start": np.asarray(g_start, dtype=np.int32),
        "gene_end": np.asarray(g_end, dtype=np.int32),
        "gene_strand": np.asarray(g_strand, dtype=np.int8),
        "gene_chrom": np.asarray(g_chrom, dtype=np.int32),
        "gene_tx_ptr": np.asarray(g_tx_ptr, dtype=np.int64),
        "gene_id_str": np.asarray(g_id_str, dtype=np.int32),
        "gene_source_str": np.asarray(g_src_str, dtype=np.int32),
        "gene_attr_str": np.asarray(g_attr_str, dtype=np.int32),
        "tx_start": np.asarray(t_start, dtype=np.int32),
        "tx_end": np.asarray(t_end, dtype=np.int32),
        "tx_strand": np.asarray(t_strand, dtype=np.int8),
        "tx_gene": np.asarray(t_gene, dtype=np.int32),
        "tx_exon_ptr": np.asarray(t_exon_ptr, dtype=np.int64),
        "tx_other_ptr": np.asarray(t_other_ptr, dtype=np.int64),
        "tx_id_str": np.asarray(t_id_str, dtype=np.int32),
        "tx_source_str": np.asarray(t_src_str, dtype=np.int32),
        "tx_attr_str": np.asarray(t_attr_str, dtype=np.int32),
        "exon_start": np.asarray(e_start, dtype=np.int32),
        "exon_end": np.asarray(e_end, dtype=np.int32),
        "exon_id_str": np.asarray(e_id_str, dtype=np.int32),
        "other_start": np.asarray(o_start, dtype=np.int32),
        "other_end": np.asarray(o_end, dtype=np.int32),
        "other_type": np.asarray(o_type, dtype=np.uint8),
        "str_offsets": str_offsets,
        "str_bytes": str_bytes,
    }

    # Pre-calculate exact aligned offsets for all sections after metadata JSON
    array_specs: dict[str, list] = {}
    # First reserve metadata size estimate using placeholder offset=0
    skeleton_meta = {
        "chrom_names": chrom_names,
        "arrays": {k: [0, int(arr.size), str(arr.dtype)] for k, arr in arrays.items()},
    }
    raw_meta = json.dumps(skeleton_meta, separators=(",", ":")).encode("utf-8")
    # Leave generous alignment room; exact offset derived from final json length with fixed-width 12-digit numbers?
    # Simpler & exact: two-pass fixed length padding in meta block (aligned to 64B, +256B headroom)
    header_end = _HEADER_PREFIX.size + len(raw_meta) + 256
    data_start = (header_end + 63) & ~63

    cursor = data_start
    for k, arr in arrays.items():
        nbytes = int(arr.nbytes)
        array_specs[k] = [cursor, int(arr.size), str(arr.dtype)]
        cursor = ((cursor + nbytes) + 63) & ~63

    final_meta = json.dumps({"chrom_names": chrom_names, "arrays": array_specs}, separators=(",", ":")).encode("utf-8")
    if _HEADER_PREFIX.size + len(final_meta) > data_start:
        raise RuntimeError("Metadata overflowed aligned header window")

    target = Path(out_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = target.with_suffix(target.suffix + f".tmp.{os.getpid()}")
    with open(tmp_path, "wb") as fh:
        fh.write(_HEADER_PREFIX.pack(IQIDX_MAGIC, len(final_meta)))
        fh.write(final_meta)
        pad_hdr = data_start - (_HEADER_PREFIX.size + len(final_meta))
        if pad_hdr:
            fh.write(b"\x00" * pad_hdr)

        pos = data_start
        for k, arr in arrays.items():
            raw = np.ascontiguousarray(arr).tobytes()
            fh.write(raw)
            pos += len(raw)
            aligned_next = (pos + 63) & ~63
            if aligned_next > pos:
                fh.write(b"\x00" * (aligned_next - pos))
                pos = aligned_next

    os.replace(tmp_path, target)
    return target


def compile_gtf_to_iqidx(gtf_path: str | Path, out_iqidx: str | Path) -> Path:
    """Stream GTF/GFF3 once in linear memory and write an immutable `.iqidx` mmap index."""
    fmt, _ = detect_annotation_format(gtf_path)
    is_gff3 = fmt == "gff3"

    chrom_order: List[str] = []
    seen_chroms: set[str] = set()
    genes_by_chrom: Dict[str, Dict[str, dict]] = defaultdict(dict)
    tx_by_id: Dict[str, dict] = {}
    tx_to_gene: Dict[str, str] = {}

    relevant_types = {"gene", "exon", *DEFAULT_RNA_TYPES, *OTHER_FEATURE_CODES}

    def _ensure_gene(chrom: str, gid: str, start: int, end: int, strand: str, source: str, attr_str: str) -> dict:
        gmap = genes_by_chrom[chrom]
        g = gmap.get(gid)
        if g is None:
            g = {
                "id": gid,
                "start": start,
                "end": end,
                "strand": strand,
                "source": source,
                "attr_str": attr_str,
                "tx_ids": [],
                "_tx_set": set(),
            }
            gmap[gid] = g
        else:
            if start > 0 and (g["start"] <= 0 or start < g["start"]):
                g["start"] = start
            if end > g["end"]:
                g["end"] = end
            if attr_str and not g["attr_str"]:
                g["attr_str"] = attr_str
            if strand != "." and g["strand"] == ".":
                g["strand"] = strand
            if source != "." and g["source"] == ".":
                g["source"] = source
        return g

    def _ensure_tx(
        chrom: str, gid: str, tid: str, start: int, end: int, strand: str, source: str, attr_str: str
    ) -> dict:
        g = _ensure_gene(chrom, gid, 0 if start <= 0 else start, end, strand, source, "")
        tx = tx_by_id.get(tid)
        if tx is None:
            tx = {
                "id": tid,
                "gene_id": gid,
                "start": start,
                "end": end,
                "strand": strand,
                "source": source,
                "attr_str": attr_str,
                "exons": [],
                "others": [],
            }
            tx_by_id[tid] = tx
            tx_to_gene[tid] = gid
        else:
            if start > 0 and (tx["start"] <= 0 or start < tx["start"]):
                tx["start"] = start
            if end > tx["end"]:
                tx["end"] = end
            if attr_str and not tx["attr_str"]:
                tx["attr_str"] = attr_str
            if strand != "." and tx["strand"] == ".":
                tx["strand"] = strand
            if source != "." and tx["source"] == ".":
                tx["source"] = source
        if tid not in g["_tx_set"]:
            g["_tx_set"].add(tid)
            g["tx_ids"].append(tid)
        return tx

    with open_text_maybe_gz(gtf_path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            cols = line.rstrip("\r\n").split("\t")
            if len(cols) < 9:
                continue
            ftype = cols[2]
            if ftype not in relevant_types:
                continue

            chrom = cols[0]
            if chrom not in seen_chroms:
                seen_chroms.add(chrom)
                chrom_order.append(chrom)

            source = cols[1]
            start = int(cols[3])
            end = int(cols[4])
            strand = cols[6]
            attrs = parse_gff3_attributes(cols[8]) if is_gff3 else parse_gtf_attributes(cols[8])

            if ftype == "gene":
                gid = attrs.get("gene_id") or attrs.get("ID")
                if gid:
                    _ensure_gene(
                        chrom, gid, start, end, strand, source,
                        _format_filtered_attrs(attrs, _SKIP_GENE_ATTRS),
                    )
            elif ftype in DEFAULT_RNA_TYPES:
                tid = attrs.get("transcript_id") or attrs.get("ID")
                parent_raw = attrs.get("gene_id") or attrs.get("Parent")
                if not tid:
                    continue
                gid = parent_raw.split(",")[0] if (is_gff3 and parent_raw) else (parent_raw or f"gene_{tid}")
                _ensure_tx(
                    chrom, gid, tid, start, end, strand, source,
                    _format_filtered_attrs(attrs, _SKIP_TX_ATTRS),
                )
            elif ftype == "exon":
                tid_raw = attrs.get("transcript_id") or attrs.get("Parent")
                if not tid_raw:
                    continue
                tids = tid_raw.split(",") if is_gff3 else [tid_raw]
                exon_id = attrs.get("exon_id") or attrs.get("ID") or ""
                for tid in tids:
                    gid = attrs.get("gene_id") or tx_to_gene.get(tid) or f"gene_{tid}"
                    tx = _ensure_tx(chrom, gid, tid, 0, 0, strand, source, "")
                    tx["exons"].append((start, end, exon_id))
            elif ftype in OTHER_FEATURE_CODES:
                tid_raw = attrs.get("transcript_id") or attrs.get("Parent")
                if not tid_raw:
                    continue
                ocode = OTHER_FEATURE_CODES[ftype]
                tids = tid_raw.split(",") if is_gff3 else [tid_raw]
                for tid in tids:
                    gid = attrs.get("gene_id") or tx_to_gene.get(tid) or f"gene_{tid}"
                    tx = _ensure_tx(chrom, gid, tid, 0, 0, strand, source, "")
                    tx["others"].append((start, end, ocode))

    return _write_iqidx_bundle(out_iqidx, chrom_order, genes_by_chrom, tx_by_id)


def compile_gffutils_db_to_iqidx(db_path: str | Path, out_iqidx: str | Path) -> Path:
    """Convert an existing gffutils SQLite database to `.iqidx` via 2 bulk SQL scans (no N+1 queries)."""
    conn = sqlite3.connect(f"file:{Path(db_path).resolve()}?mode=ro", uri=True)
    try:
        cur = conn.cursor()
        chrom_order: List[str] = []
        seen_chroms: set[str] = set()
        genes_by_chrom: Dict[str, Dict[str, dict]] = defaultdict(dict)
        tx_by_id: Dict[str, dict] = {}
        feat_parent: Dict[str, List[str]] = defaultdict(list)

        for parent_id, child_id in cur.execute(
            "SELECT parent, child FROM relations WHERE level = 1 AND parent != child"
        ):
            feat_parent[child_id].append(parent_id)

        raw_children: List[Tuple[str, str, str, int, int, str, str, dict]] = []
        for fid, seqid, source, ftype, start, end, strand, raw_attrs in cur.execute(
            "SELECT id, seqid, source, featuretype, start, end, strand, attributes FROM features"
        ):
            if ftype != "gene" and ftype not in DEFAULT_RNA_TYPES and ftype != "exon" and ftype not in OTHER_FEATURE_CODES:
                continue
            if seqid not in seen_chroms:
                seen_chroms.add(seqid)
                chrom_order.append(seqid)
            attrs_json = json.loads(raw_attrs) if raw_attrs else {}
            flat_attrs = {k: str(v[0]) for k, v in attrs_json.items() if v}

            if ftype == "gene":
                genes_by_chrom[seqid][fid] = {
                    "id": fid,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "source": source,
                    "attr_str": _format_filtered_attrs(flat_attrs, _SKIP_GENE_ATTRS),
                    "tx_ids": [],
                    "_tx_set": set(),
                }
            elif ftype in DEFAULT_RNA_TYPES:
                parents = [p for p in feat_parent.get(fid, []) if p != fid]
                gid = flat_attrs.get("gene_id") or (parents[0] if parents else f"gene_{fid}")
                tx_by_id[fid] = {
                    "id": fid,
                    "gene_id": gid,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "source": source,
                    "attr_str": _format_filtered_attrs(flat_attrs, _SKIP_TX_ATTRS),
                    "exons": [],
                    "others": [],
                    "_chrom": seqid,
                }
            else:
                raw_children.append((fid, seqid, ftype, int(start), int(end), strand, source, flat_attrs))

        for tid, tx in tx_by_id.items():
            chrom = tx["_chrom"]
            gid = tx["gene_id"]
            g = genes_by_chrom[chrom].get(gid)
            if g is None:
                g = {
                    "id": gid,
                    "start": tx["start"],
                    "end": tx["end"],
                    "strand": tx["strand"],
                    "source": tx["source"],
                    "attr_str": "",
                    "tx_ids": [],
                    "_tx_set": set(),
                }
                genes_by_chrom[chrom][gid] = g
            if tid not in g["_tx_set"]:
                g["_tx_set"].add(tid)
                g["tx_ids"].append(tid)

        for fid, seqid, ftype, start, end, strand, source, flat_attrs in raw_children:
            parents = feat_parent.get(fid, [])
            if not parents and "transcript_id" in flat_attrs:
                parents = [flat_attrs["transcript_id"]]
            for pid in parents:
                tx = tx_by_id.get(pid)
                if tx is None:
                    continue
                if ftype == "exon":
                    tx["exons"].append((start, end, flat_attrs.get("exon_id", "")))
                else:
                    tx["others"].append((start, end, OTHER_FEATURE_CODES[ftype]))
    finally:
        conn.close()

    return _write_iqidx_bundle(out_iqidx, chrom_order, genes_by_chrom, tx_by_id)


def load_or_build_iqidx(source_path: str | Path, cache_dir: Optional[str | Path] = None) -> Path:
    """Resolve `.iqidx` path, compiling GTF/GFF3/SQLite DB if necessary."""
    fmt, _ = detect_annotation_format(source_path)
    p = Path(source_path)
    if fmt == "iqidx":
        return p
    target_dir = Path(cache_dir) if cache_dir else p.parent
    target_iqidx = target_dir / (p.name + ".iqidx")
    if target_iqidx.exists() and target_iqidx.stat().st_mtime >= p.stat().st_mtime:
        return target_iqidx
    if fmt == "db":
        return compile_gffutils_db_to_iqidx(p, target_iqidx)
    return compile_gtf_to_iqidx(p, target_iqidx)


@dataclass(frozen=True, slots=True)
class GeneClusterSlice:
    chr_id: str
    gene_indices: Tuple[int, ...]
    start: int
    end: int


class MmapAnnotationIndex:
    """Read-only mmap CSR index shared zero-copy across worker processes (PSS shared cleanly)."""

    supports_fast_id_lookup: bool = True

    def __init__(self, iqidx_path: str | Path) -> None:
        self.path = Path(iqidx_path)
        self._fd = os.open(self.path, os.O_RDONLY)
        try:
            self._mmap = mmap.mmap(self._fd, 0, access=mmap.ACCESS_READ)
            magic, meta_len = _HEADER_PREFIX.unpack_from(self._mmap, 0)
            if magic != IQIDX_MAGIC:
                raise ValueError(f"Invalid .iqidx magic header in {self.path}")
            meta = json.loads(bytes(self._mmap[_HEADER_PREFIX.size : _HEADER_PREFIX.size + meta_len]).decode("utf-8"))
            self.chrom_names: Tuple[str, ...] = tuple(meta["chrom_names"])
            self.chrom_to_id: Dict[str, int] = {name: idx for idx, name in enumerate(self.chrom_names)}
            self._arrays: Dict[str, np.ndarray] = {}
            for name, (offset, count, dtype_str) in meta["arrays"].items():
                self._arrays[name] = np.frombuffer(
                    self._mmap, dtype=np.dtype(dtype_str), count=count, offset=offset
                )
        except Exception:
            self.close()
            raise

    def __enter__(self) -> "MmapAnnotationIndex":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    def close(self) -> None:
        mm = getattr(self, "_mmap", None)
        if mm is not None:
            self._arrays.clear()
            try:
                mm.close()
            except BufferError:
                import gc

                gc.collect()
                try:
                    mm.close()
                except BufferError:
                    pass
            self._mmap = None
        fd = getattr(self, "_fd", -1)
        if fd >= 0:
            os.close(fd)
            self._fd = -1

    def __getattr__(self, name: str) -> np.ndarray:
        arrays = object.__getattribute__(self, "_arrays")
        if name in arrays:
            return arrays[name]
        raise AttributeError(name)

    def get_str(self, str_idx: int) -> str:
        if str_idx <= 0:
            return ""
        offsets = self._arrays["str_offsets"]
        raw = self._arrays["str_bytes"]
        s0 = int(offsets[str_idx])
        s1 = int(offsets[str_idx + 1])
        return bytes(raw[s0:s1]).decode("utf-8") if s1 > s0 else ""

    def seqids(self) -> List[str]:
        return list(self.chrom_names)

    def chrom_gene_slice(self, chr_id: str) -> Tuple[int, int]:
        cid = self.chrom_to_id.get(chr_id)
        if cid is None:
            return 0, 0
        ptr = self._arrays["chrom_gene_ptr"]
        return int(ptr[cid]), int(ptr[cid + 1])

    def iter_gene_ids(self) -> Iterable[str]:
        for idx in self._arrays["gene_id_str"]:
            yield self.get_str(int(idx))

    def iter_transcript_ids(self) -> Iterable[str]:
        for idx in self._arrays["tx_id_str"]:
            yield self.get_str(int(idx))

    def forbidden_numeric_ids(self, chr_id: str) -> set[int]:
        forbidden: set[int] = set()
        g0, g1 = self.chrom_gene_slice(chr_id)
        if g0 >= g1:
            return forbidden
        g_ids = self._arrays["gene_id_str"]
        g_tx = self._arrays["gene_tx_ptr"]
        t_ids = self._arrays["tx_id_str"]
        tx_prefix_len = len(TranscriptNaming.transcript_prefix)

        for gi in range(g0, g1):
            gid = self.get_str(int(g_ids[gi]))
            if gid.startswith(TranscriptNaming.novel_gene_prefix):
                try:
                    forbidden.add(int(gid.split("_")[-1]))
                except (IndexError, ValueError):
                    pass
            t0, t1 = int(g_tx[gi]), int(g_tx[gi + 1])
            for ti in range(t0, t1):
                tid = self.get_str(int(t_ids[ti]))
                if tid.startswith(TranscriptNaming.transcript_prefix):
                    try:
                        forbidden.add(int(tid.split(".")[0][tx_prefix_len:]))
                    except (IndexError, ValueError):
                        pass
        return forbidden

    def feature_id_map(self, chr_id: str, feature: str = "exon") -> Dict[Tuple[str, int, int, str], str]:
        if feature != "exon":
            return {}
        out: Dict[Tuple[str, int, int, str], str] = {}
        g0, g1 = self.chrom_gene_slice(chr_id)
        if g0 >= g1:
            return out
        g_tx = self._arrays["gene_tx_ptr"]
        t_strand = self._arrays["tx_strand"]
        t_ex = self._arrays["tx_exon_ptr"]
        e_s = self._arrays["exon_start"]
        e_e = self._arrays["exon_end"]
        e_id = self._arrays["exon_id_str"]
        for gi in range(g0, g1):
            for ti in range(int(g_tx[gi]), int(g_tx[gi + 1])):
                strand = _I8_TO_STRAND[int(t_strand[ti])]
                for ei in range(int(t_ex[ti]), int(t_ex[ti + 1])):
                    sid = int(e_id[ei])
                    if sid > 0:
                        out[(chr_id, int(e_s[ei]), int(e_e[ei]), strand)] = self.get_str(sid)
        return out

    def iter_gene_clusters(self, chr_id: str) -> Iterable[GeneClusterSlice]:
        """Yield overlap-connected gene clusters along sorted chromosome coordinates."""
        g0, g1 = self.chrom_gene_slice(chr_id)
        if g0 >= g1:
            return
        starts = self._arrays["gene_start"]
        ends = self._arrays["gene_end"]

        cur_indices: List[int] = [g0]
        cur_start = int(starts[g0])
        cur_end = int(ends[g0])

        for gi in range(g0 + 1, g1):
            gs = int(starts[gi])
            ge = int(ends[gi])
            if gs <= cur_end:
                cur_indices.append(gi)
                if ge > cur_end:
                    cur_end = ge
            else:
                yield GeneClusterSlice(chr_id, tuple(cur_indices), cur_start, cur_end)
                cur_indices = [gi]
                cur_start = gs
                cur_end = ge
        yield GeneClusterSlice(chr_id, tuple(cur_indices), cur_start, cur_end)
