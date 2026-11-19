"""Static sorted interval index supporting zero-copy `.iqidx` and duck-typed DBs."""

from __future__ import annotations

from bisect import bisect_left
from collections import defaultdict
from dataclasses import dataclass, field
import logging
from typing import Any, Dict, List, Optional, Sequence, Tuple

from isoquant_next.annotation.gtf_validator import parse_gff3_attributes, parse_gtf_attributes
from isoquant_next.annotation.iqidx import MmapAnnotationIndex

logger = logging.getLogger("IsoQuant")

_BLOCK_SIZE = 64


@dataclass(slots=True)
class IqidxFeatureView:
    """Lightweight feature duck-type backed by `.iqidx` CSR offsets."""

    id: str
    chrom: str
    start: int
    end: int
    strand: str
    featuretype: str
    _attr_str: str
    _gene_idx: int = -1
    _attrs_cache: Optional[Dict[str, List[str]]] = field(default=None, repr=False)

    @property
    def seqid(self) -> str:
        return self.chrom

    @property
    def attributes(self) -> Dict[str, List[str]]:
        if self._attrs_cache is None:
            parsed: Dict[str, List[str]] = {}
            if self._attr_str:
                raw = (
                    parse_gff3_attributes(self._attr_str)
                    if "=" in self._attr_str and '"' not in self._attr_str
                    else parse_gtf_attributes(self._attr_str)
                )
                for k, v in raw.items():
                    parsed[k] = [v]
            if self.featuretype == "gene":
                parsed.setdefault("gene_id", [self.id])
            self._attrs_cache = parsed
        return self._attrs_cache


class _SortedIntervalIndex:
    """Static sorted interval container with block-max pruning matching `[s, e+1)` half-open queries."""

    __slots__ = ("starts", "ends", "values", "block_max_ends")

    def __init__(self, items: List[Tuple[int, int, Any]]) -> None:
        if items:
            items.sort(key=lambda it: (it[0], it[1]))
            self.starts = [it[0] for it in items]
            self.ends = [it[1] for it in items]
            self.values = [it[2] for it in items]
            self.block_max_ends = [
                max(self.ends[b : b + _BLOCK_SIZE])
                for b in range(0, len(self.ends), _BLOCK_SIZE)
            ]
        else:
            self.starts = []
            self.ends = []
            self.values = []
            self.block_max_ends = []

    def __len__(self) -> int:
        return len(self.starts)

    def __iter__(self):
        return zip(self.starts, self.ends, self.values)

    def query(self, start_q: int, end_q: int) -> List[Any]:
        """Return payloads for all `[s, e]` overlapping half-open `[start_q, end_q)`."""
        if end_q <= start_q or not self.starts:
            return []
        limit = bisect_left(self.starts, end_q)
        if limit <= 0:
            return []
        hits: List[Any] = []
        ends = self.ends
        vals = self.values
        bmax = self.block_max_ends
        num_blocks = (limit + _BLOCK_SIZE - 1) // _BLOCK_SIZE
        for blk in range(num_blocks):
            if bmax[blk] < start_q:
                continue
            lo = blk * _BLOCK_SIZE
            hi = min(lo + _BLOCK_SIZE, limit)
            for i in range(lo, hi):
                if ends[i] >= start_q:
                    hits.append(vals[i])
        return hits


class GenomicIntervalIndex:
    """Static sorted interval index mapping genomic coordinates to genes and exons."""

    def __init__(self, genedb: Any, chromosomes: Optional[List[str]] = None) -> None:
        self.db = genedb
        self.chromosomes = set(chromosomes) if chromosomes else None
        self.gene_indices: Dict[str, _SortedIntervalIndex] = {}
        self.exon_indices: Dict[str, _SortedIntervalIndex] = {}
        self._features_by_id: Dict[str, Any] = {}
        self._genes_by_name: Dict[str, List[Any]] = defaultdict(list)
        self._is_iqidx = isinstance(genedb, MmapAnnotationIndex)
        self._build_indices()

    def _build_indices(self) -> None:
        if self._is_iqidx:
            self._build_from_iqidx(self.db)
        else:
            self._build_from_featuredb(self.db)

    def _build_from_iqidx(self, idx: MmapAnnotationIndex) -> None:
        strand_map = {1: "+", -1: "-", 0: "."}
        g_start = idx.gene_start
        g_end = idx.gene_end
        g_strand = idx.gene_strand
        g_tx_ptr = idx.gene_tx_ptr
        g_id_str = idx.gene_id_str
        g_attr_str = idx.gene_attr_str
        t_ex_ptr = idx.tx_exon_ptr
        e_start = idx.exon_start
        e_end = idx.exon_end
        e_id_str = idx.exon_id_str
        t_strand = idx.tx_strand

        gene_count = 0
        exon_count = 0
        for chrom in idx.seqids():
            if self.chromosomes and chrom not in self.chromosomes:
                continue
            g0, g1 = idx.chrom_gene_slice(chrom)
            if g0 >= g1:
                continue
            gene_items: List[Tuple[int, int, str]] = []
            exon_items: List[Tuple[int, int, str]] = []
            seen_exons: set[Tuple[int, int, str]] = set()
            for gi in range(g0, g1):
                gid = idx.get_str(int(g_id_str[gi]))
                gs = int(g_start[gi])
                ge = int(g_end[gi])
                gstr = strand_map.get(int(g_strand[gi]), ".")
                gattrs = idx.get_str(int(g_attr_str[gi]))
                gf = IqidxFeatureView(
                    id=gid,
                    chrom=chrom,
                    start=gs,
                    end=ge,
                    strand=gstr,
                    featuretype="gene",
                    _attr_str=gattrs,
                    _gene_idx=gi,
                )
                self._features_by_id[gid] = gf
                gname = gf.attributes.get("gene_name", [None])[0]
                if gname:
                    self._genes_by_name[gname].append(gf)
                gene_items.append((gs, ge, gid))
                gene_count += 1

                for ti in range(int(g_tx_ptr[gi]), int(g_tx_ptr[gi + 1])):
                    tstr = strand_map.get(int(t_strand[ti]), ".")
                    for ei in range(int(t_ex_ptr[ti]), int(t_ex_ptr[ti + 1])):
                        es = int(e_start[ei])
                        ee = int(e_end[ei])
                        ex_id = idx.get_str(int(e_id_str[ei])) or f"{chrom}:{es}-{ee}:{tstr}"
                        key = (es, ee, ex_id)
                        if key in seen_exons:
                            continue
                        seen_exons.add(key)
                        if ex_id not in self._features_by_id:
                            self._features_by_id[ex_id] = IqidxFeatureView(
                                id=ex_id,
                                chrom=chrom,
                                start=es,
                                end=ee,
                                strand=tstr,
                                featuretype="exon",
                                _attr_str=f'gene_id "{gid}";',
                                _gene_idx=gi,
                            )
                        exon_items.append((es, ee, ex_id))
                        exon_count += 1

            if gene_items:
                self.gene_indices[chrom] = _SortedIntervalIndex(gene_items)
            if exon_items:
                self.exon_indices[chrom] = _SortedIntervalIndex(exon_items)

        logger.info("Indexed %d genes across %d chromosomes and %d exons (.iqidx)", gene_count, len(self.gene_indices), exon_count)

    def _build_from_featuredb(self, db: Any) -> None:
        raw_genes: Dict[str, List[Tuple[int, int, str]]] = defaultdict(list)
        raw_exons: Dict[str, List[Tuple[int, int, str]]] = defaultdict(list)
        try:
            for gene in db.features_of_type("gene"):
                chrom = str(gene.chrom)
                if self.chromosomes and chrom not in self.chromosomes:
                    continue
                raw_genes[chrom].append((int(gene.start), int(gene.end), gene.id))
                attrs = getattr(gene, "attributes", {}) or {}
                gname = attrs.get("gene_name", [None])[0]
                if gname:
                    self._genes_by_name[gname].append(gene)
        except Exception as exc:
            logger.warning("Error building gene index: %s", exc)

        try:
            for exon in db.features_of_type("exon"):
                chrom = str(exon.chrom)
                if self.chromosomes and chrom not in self.chromosomes:
                    continue
                raw_exons[chrom].append((int(exon.start), int(exon.end), exon.id))
        except Exception as exc:
            logger.warning("Error building exon index: %s", exc)

        for chrom, items in raw_genes.items():
            self.gene_indices[chrom] = _SortedIntervalIndex(items)
        for chrom, items in raw_exons.items():
            self.exon_indices[chrom] = _SortedIntervalIndex(items)

    def lookup_feature(self, feature_id: str) -> Any:
        if self._is_iqidx:
            feat = self._features_by_id.get(feature_id)
            if feat is None:
                raise KeyError(feature_id)
            return feat
        return self.db[feature_id]

    def _resolve_features(self, chrom: str, pos: int, ids: Sequence[str], kind: str) -> list:
        out = []
        for fid in ids:
            try:
                out.append(self.lookup_feature(fid))
            except Exception:
                logger.debug("Failed to resolve %s %r at %s:%s", kind, fid, chrom, pos, exc_info=True)
        return out

    def get_genes_at(self, chrom: str, pos: int, window: Optional[int] = None) -> list:
        """Return all gene features overlapping `pos` or `[max(1, pos-window), pos+window)`."""
        idx = self.gene_indices.get(chrom)
        if idx is None:
            return []
        if window is None:
            ids = idx.query(pos, pos + 1)
        else:
            ids = idx.query(max(1, pos - window), pos + window)
        return self._resolve_features(chrom, pos, ids, "gene")

    def get_exons_at(self, chrom: str, pos: int, window: Optional[int] = None) -> list:
        """Return all exon features overlapping `pos` or `[max(1, pos-window), pos+window)`."""
        idx = self.exon_indices.get(chrom)
        if idx is None:
            return []
        if window is None:
            ids = idx.query(pos, pos + 1)
        else:
            ids = idx.query(max(1, pos - window), pos + window)
        return self._resolve_features(chrom, pos, ids, "exon")

    def get_exons_of_gene(self, gene_feature: Any) -> list:
        """Return all exon children for `gene_feature` ordered by start coordinate."""
        if self._is_iqidx and isinstance(gene_feature, IqidxFeatureView) and gene_feature._gene_idx >= 0:
            gi = gene_feature._gene_idx
            idx: MmapAnnotationIndex = self.db
            g_tx_ptr = idx.gene_tx_ptr
            t_ex_ptr = idx.tx_exon_ptr
            t_strand = idx.tx_strand
            e_start = idx.exon_start
            e_end = idx.exon_end
            e_id_str = idx.exon_id_str
            strand_map = {1: "+", -1: "-", 0: "."}
            exons: List[IqidxFeatureView] = []
            seen: set[Tuple[int, int, str]] = set()
            for ti in range(int(g_tx_ptr[gi]), int(g_tx_ptr[gi + 1])):
                tstr = strand_map.get(int(t_strand[ti]), ".")
                for ei in range(int(t_ex_ptr[ti]), int(t_ex_ptr[ti + 1])):
                    es = int(e_start[ei])
                    ee = int(e_end[ei])
                    ex_id = idx.get_str(int(e_id_str[ei])) or f"{gene_feature.chrom}:{es}-{ee}:{tstr}"
                    key = (es, ee, ex_id)
                    if key in seen:
                        continue
                    seen.add(key)
                    exons.append(
                        IqidxFeatureView(
                            id=ex_id,
                            chrom=gene_feature.chrom,
                            start=es,
                            end=ee,
                            strand=tstr,
                            featuretype="exon",
                            _attr_str=f'gene_id "{gene_feature.id}";',
                            _gene_idx=gi,
                        )
                    )
            exons.sort(key=lambda e: (e.start, e.end))
            return exons
        try:
            return list(self.db.children(gene_feature, featuretype="exon", order_by="start"))
        except Exception:
            return []

    def get_all_genes(self) -> list:
        """Return all gene features stored in the index."""
        genes = []
        for chrom, idx in self.gene_indices.items():
            for gene_id in idx.values:
                try:
                    genes.append(self.lookup_feature(gene_id))
                except (KeyError, AttributeError, TypeError):
                    logger.debug("Failed to resolve gene %r on %s", gene_id, chrom, exc_info=True)
        return genes

    def find_genes_by_name(self, gene_name: str) -> list:
        """Return all gene features whose `gene_name` attribute matches `gene_name`."""
        if gene_name in self._genes_by_name:
            return list(self._genes_by_name[gene_name])
        matching = []
        for gene in self.get_all_genes():
            attrs = getattr(gene, "attributes", {}) or {}
            if attrs.get("gene_name", [None])[0] == gene_name:
                matching.append(gene)
        return matching
