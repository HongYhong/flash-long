"""GeneInfo cluster representations backed by mmap `.iqidx` or legacy gffutils FeatureDB."""

from __future__ import annotations

import logging
from collections import defaultdict
from enum import Enum
from functools import partial
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from isoquant_next.annotation.iqidx import (
    MmapAnnotationIndex,
    OTHER_FEATURE_NAMES,
    _I8_TO_STRAND,
)
from isoquant_next.core.serialization import read_int, read_string, write_int, write_string
from isoquant_next.core.utils import (
    contains,
    equal_ranges,
    interval_len,
    junctions_from_blocks,
    overlaps,
)

logger = logging.getLogger("IsoQuant")


class TranscriptModelType(Enum):
    known = 1
    novel_in_catalog = 2
    novel_not_in_catalog = 10


class ExonRegion:
    def __init__(self, chr_id: str, start: int, end: int, strand: str) -> None:
        self.chr_id: str = chr_id
        self.start: int = start
        self.end: int = end
        self.strand: str = strand
        self.member_exon_indices: list[int] = []
        self.gene_ids: set[str] = set()

    @property
    def id(self) -> Tuple[str, int, int, str]:
        return (self.chr_id, self.start, self.end, self.strand)


class ExonSpliceSiteRegion:
    def __init__(self, chr_id: str, start: int, end: int, strand: str, gene_id: str) -> None:
        self.chr_id: str = chr_id
        self.start: int = start
        self.end: int = end
        self.strand: str = strand
        self.gene_id: str = gene_id
        self.candidates: list[Tuple[int, int]] = []
        self.left_unique: list[bool] = []
        self.right_unique: list[bool] = []


class FeatureInfo:
    def __init__(self, chr_id: str, start: int, end: int, strand: str, ftype: str, gene_ids: list[str]) -> None:
        self.id = (chr_id, start, end, strand)
        self.chr_id = chr_id
        self.start = start
        self.end = end
        self.strand = strand
        self.type = ftype
        self.gene_ids = gene_ids

    @staticmethod
    def header() -> str:
        return "chr\tstart\tend\tstrand\tflags\tgene_ids"

    def to_str(self) -> str:
        return f"{self.chr_id}\t{self.start}\t{self.end}\t{self.strand}\t{self.type}\t{','.join(self.gene_ids)}"


class GeneList:
    def __init__(self, gene_id_list, delta: int, chr_id: str, start: int, end: int) -> None:
        self.gene_id_set = set(gene_id_list)
        self.delta = delta
        self.chr_id = chr_id
        self.start = start
        self.end = end

    def merge(self, other: "GeneList") -> None:
        assert self.chr_id == other.chr_id
        self.gene_id_set.update(other.gene_id_set)
        self.start = min(self.start, other.start)
        self.end = max(self.end, other.end)

    def overlaps(self, other: "GeneList") -> bool:
        return self.chr_id == other.chr_id and bool(self.gene_id_set & other.gene_id_set)

    @classmethod
    def deserialize(cls, infile) -> "GeneList":
        obj = cls.__new__(cls)
        obj.delta = read_int(infile)
        gene_count = read_int(infile)
        obj.gene_id_set = {read_string(infile) for _ in range(gene_count)}
        obj.chr_id = read_string(infile)
        obj.start = read_int(infile)
        obj.end = read_int(infile)
        return obj

    def serialize(self, outfile) -> None:
        write_int(self.delta, outfile)
        write_int(len(self.gene_id_set), outfile)
        for g in self.gene_id_set:
            write_string(g.id if hasattr(g, "id") else str(g), outfile)
        write_string(self.chr_id, outfile)
        write_int(self.start, outfile)
        write_int(self.end, outfile)


class FeatureProfiles:
    """Stores sorted reference features and {-2, -1, 1} isoform profiles."""

    def __init__(self) -> None:
        self.profiles: Dict[str, List[int]] = {}
        self.profile_ranges: Dict[str, Tuple[int, int]] = {}
        self.features: List[Tuple[int, int]] = []

    def set_features(self, features: List[Tuple[int, int]]) -> None:
        self.features = features

    def set_profiles(
        self,
        transcript_id: str,
        transcript_features: Sequence[Tuple[int, int]],
        transcript_region: Tuple[int, int],
        comparator,
    ) -> None:
        n = len(self.features)
        t_start, t_end = transcript_region
        profile = [-1 if not (f[1] < t_start or f[0] > t_end) else -2 for f in self.features]
        self.profiles[transcript_id] = profile

        pos = 0
        for feat in transcript_features:
            while pos < n and not comparator(feat, self.features[pos]):
                pos += 1
            while pos < n and comparator(feat, self.features[pos]):
                profile[pos] = 1
                pos += 1

        start_pos = 0
        while start_pos < n and profile[start_pos] < 1:
            start_pos += 1
        end_pos = n - 1
        while end_pos >= 0 and profile[end_pos] < 1:
            end_pos -= 1
        self.profile_ranges[transcript_id] = (start_pos, end_pos + 1)


class TranscriptModel:
    def __init__(
        self,
        chr_id: str,
        strand: str,
        transcript_id: str,
        gene_id: str,
        exon_blocks: List[Tuple[int, int]],
        transcript_type: TranscriptModelType,
        other_features=None,
        source: str = "IsoQuant",
    ) -> None:
        self.chr_id = chr_id
        self.strand = strand
        self.transcript_id = transcript_id
        self.gene_id = gene_id
        self.exon_blocks = exon_blocks
        self.transcript_type = transcript_type
        self.source = source
        self.other_features = [] if not other_features else list(other_features)
        self.additional_info: Dict[str, str] = {}
        self.intron_path: Tuple[Tuple[int, int], ...] = ()
        self.polya_confirmed = False
        self.tss_confirmed = False

    @classmethod
    def from_reference_transcript(cls, gene_info: "GeneInfo", isoform_id: str) -> "TranscriptModel":
        tm = cls.__new__(cls)
        tm.chr_id = gene_info.chr_id
        tm.strand = gene_info.isoform_strands[isoform_id]
        tm.transcript_id = isoform_id
        tm.gene_id = gene_info.gene_id_map[isoform_id]
        tm.exon_blocks = gene_info.all_isoforms_exons[isoform_id]
        tm.transcript_type = TranscriptModelType.known
        tm.source = gene_info.sources[isoform_id]
        tm.other_features = gene_info.other_features[isoform_id]
        tm.additional_info = {}
        tm.intron_path = ()
        tm.polya_confirmed = False
        tm.tss_confirmed = False
        return tm

    def get_start(self) -> int:
        return self.exon_blocks[0][0]

    def get_end(self) -> int:
        return self.exon_blocks[-1][1]

    def check_additional(self, attribute: str) -> bool:
        return attribute in self.additional_info

    def add_additional_attribute(self, attribute: str, value: str) -> None:
        self.additional_info[attribute] = value

    def additional_attributes_str(self) -> str:
        return " ".join(f'{k} "{v}";' for k, v in self.additional_info.items())


class StrandDetector:
    def __init__(self, chr_record) -> None:
        self.strand_dict: Dict[Tuple[int, int], str] = {}
        self.chr_record = chr_record

    def set_strand(self, intron: Tuple[int, int], strand: Optional[str] = None) -> None:
        from isoquant_next.core.utils import get_intron_strand

        if strand:
            self.strand_dict[intron] = strand
        elif self.chr_record is not None:
            self.strand_dict[intron] = get_intron_strand(intron, self.chr_record)

    def count_canonical_sites(self, introns: Iterable[Tuple[int, int]]) -> Tuple[int, int]:
        from isoquant_next.core.utils import get_intron_strand

        count_fwd = 0
        count_rev = 0
        for intron in introns:
            strand = self.strand_dict.get(intron)
            if strand is None:
                strand = get_intron_strand(intron, self.chr_record) if self.chr_record is not None else "."
                self.strand_dict[intron] = strand
            if strand == "+":
                count_fwd += 1
            elif strand == "-":
                count_rev += 1
        return count_fwd, count_rev

    def get_clean_strand(self, introns: Iterable[Tuple[int, int]]) -> str:
        count_fwd, count_rev = self.count_canonical_sites(introns)
        if count_fwd == 0 and count_rev > 0:
            return "-"
        if count_fwd > 0 and count_rev == 0:
            return "+"
        return "."

    def get_strand(self, introns: Iterable[Tuple[int, int]], has_polya: bool = False, has_polyt: bool = False) -> str:
        count_fwd, count_rev = self.count_canonical_sites(introns)
        if count_fwd == count_rev:
            if has_polya and not has_polyt:
                return "+"
            if has_polyt and not has_polya:
                return "-"
            return "."
        return "+" if count_rev < count_fwd else "-"


class GeneInfo:
    """Cluster of overlapping reference genes with sweep-line feature properties and .iqidx support."""

    EXTRA_BASES_FOR_SEQ = 20
    OTHER_FEATURES = {"UTR", "stop_codon", "start_codon", "CDS"}
    DEFAULT_RNA_TYPES = (
        "transcript", "RNA", "mRNA", "miRNA", "ncRNA", "tRNA", "rRNA",
        "snRNA", "lnc_RNA", "snoRNA", "telomerase_RNA", "Y_RNA", "scRNA",
    )

    def __init__(self, gene_db_list, db, delta: int = 0, prepare_profiles: bool = True) -> None:
        if db is None:
            return
        assert gene_db_list
        self.db = db
        self.gene_db_list = gene_db_list
        self.chr_id, self.start, self.end = self.get_gene_region()
        self.delta = delta
        self.all_read_region_start = self.start
        self.all_read_region_end = self.end
        self.canonical_sites: Dict = {}
        self.gene_regions: Dict[str, Tuple[int, int]] = {}
        self.reference_region: Optional[str] = None
        self.intron_profiles = FeatureProfiles()
        self.exon_profiles = FeatureProfiles()
        self.split_exon_profiles = FeatureProfiles()
        self.isoform_strands: Dict[str, str] = {}
        self.gene_strands: Dict[str, str] = {}
        self.other_features: Dict[str, list] = defaultdict(list)
        self.sources: Dict[str, str] = {}
        self.gene_id_map: Dict[str, str] = {}
        self.feature_attributes: Dict[str, str] = defaultdict(str)
        self.all_isoforms_introns: Dict[str, List[Tuple[int, int]]] = {}
        self.all_isoforms_exons: Dict[str, List[Tuple[int, int]]] = {}
        self.exon_property_map: Optional[List[FeatureInfo]] = None
        self.intron_property_map: Optional[List[FeatureInfo]] = None

        self._populate_from_gffutils_single_pass()
        self._finalize_profiles_and_properties(prepare_profiles)

    def _populate_from_gffutils_single_pass(self) -> None:
        """Eliminates 5+ redundant db.children() queries & O(T*E) attribute loop per GeneInfo."""
        skip_g = {"gene_id", "ID", "level", "Parent"}
        skip_t = {"transcript_id", "gene_id", "ID", "level", "exons", "Parent"}
        skip_e = {"transcript_id", "gene_id", "ID", "Parent", "level", "exon_id", "exon", "exon_number"}
        used_cds = False

        for gene_db in self.gene_db_list:
            self.gene_strands[gene_db.id] = gene_db.strand
            self.sources[gene_db.id] = gene_db.source
            g_attr_parts = [
                f'{k} "{v[0]}"; ' for k, v in gene_db.attributes.items() if k not in skip_g and v
            ]
            if g_attr_parts:
                self.feature_attributes[gene_db.id] = "".join(g_attr_parts)

            gene_exons: List[Tuple[int, int, str]] = []
            tx_attrs_list: List[Tuple[str, str]] = []
            for t in self.db.children(gene_db, featuretype=self.DEFAULT_RNA_TYPES, order_by="start"):
                tid = t.id
                self.isoform_strands[tid] = t.strand
                self.sources[tid] = t.source
                self.gene_id_map[tid] = gene_db.id

                t_attr_str = "".join(
                    f'{k} "{v[0]}"; ' for k, v in t.attributes.items() if k not in skip_t and v
                )
                if t_attr_str:
                    self.feature_attributes[tid] = t_attr_str
                tx_exon_attr_str = "".join(
                    f'{k} "{v[0]}"; ' for k, v in t.attributes.items() if k not in skip_e and v
                )
                if tx_exon_attr_str:
                    tx_attrs_list.append((tid, tx_exon_attr_str))

                all_children = []
                for e in self.db.children(t, order_by="start"):
                    if e.featuretype == "exon":
                        gene_exons.append((e.start, e.end, e.strand))
                        all_children.append(("exon", e.start, e.end))
                    elif e.featuretype == "CDS":
                        all_children.append(("CDS", e.start, e.end))
                        self.other_features[tid].append((e.start, e.end, e.featuretype))
                    elif e.featuretype in self.OTHER_FEATURES:
                        self.other_features[tid].append((e.start, e.end, e.featuretype))

                exons: List[Tuple[int, int]] = []
                n_ch = len(all_children)
                for idx, (ftype, s, e_coord) in enumerate(all_children):
                    if ftype == "exon":
                        exons.append((s, e_coord))
                        continue
                    prev_c = None if idx == 0 else (all_children[idx - 1][1], all_children[idx - 1][2])
                    next_c = None if idx == n_ch - 1 else (all_children[idx + 1][1], all_children[idx + 1][2])
                    if (prev_c and overlaps((s, e_coord), prev_c)) or (next_c and overlaps((s, e_coord), next_c)):
                        continue
                    exons.append((s, e_coord))
                    used_cds = True

                if exons:
                    self.all_isoforms_exons[tid] = exons
                    self.all_isoforms_introns[tid] = junctions_from_blocks(exons)

            if tx_attrs_list and gene_exons:
                for tid, e_attr in tx_attrs_list:
                    for es, ee, estrand in gene_exons:
                        self.feature_attributes[f"{tid}_{es}_{ee}_{estrand}"] = e_attr

        if used_cds:
            logger.warning("Gene %s uses non-overlapping CDS as exons", self.gene_db_list[0].id)

    @classmethod
    def from_iqidx_cluster(
        cls,
        iqidx: MmapAnnotationIndex,
        gene_indices: Sequence[int],
        delta: int = 0,
        prepare_profiles: bool = True,
    ) -> "GeneInfo":
        """Build complete GeneInfo directly from mmap `.iqidx` views without SQLite queries."""
        gi = cls.__new__(cls)
        gi.db = iqidx
        gi.gene_db_list = []
        gi.delta = delta
        gi.reference_region = None
        gi.canonical_sites = {}
        gi.gene_regions = {}
        gi.intron_profiles = FeatureProfiles()
        gi.exon_profiles = FeatureProfiles()
        gi.split_exon_profiles = FeatureProfiles()
        gi.isoform_strands = {}
        gi.gene_strands = {}
        gi.other_features = defaultdict(list)
        gi.sources = {}
        gi.gene_id_map = {}
        gi.feature_attributes = defaultdict(str)
        gi.all_isoforms_introns = {}
        gi.all_isoforms_exons = {}
        gi.exon_property_map = None
        gi.intron_property_map = None

        first_idx = gene_indices[0]
        gi.chr_id = iqidx.chrom_names[int(iqidx.gene_chrom[first_idx])]
        gi.start = min(int(iqidx.gene_start[idx]) for idx in gene_indices)
        gi.end = max(int(iqidx.gene_end[idx]) for idx in gene_indices)
        gi.all_read_region_start = gi.start
        gi.all_read_region_end = gi.end

        g_tx_ptr = iqidx.gene_tx_ptr
        t_ex_ptr = iqidx.tx_exon_ptr
        t_ot_ptr = iqidx.tx_other_ptr
        e_start = iqidx.exon_start
        e_end = iqidx.exon_end
        o_start = iqidx.other_start
        o_end = iqidx.other_end
        o_type = iqidx.other_type

        for g_idx in gene_indices:
            gid = iqidx.get_str(int(iqidx.gene_id_str[g_idx]))
            gi.gene_strands[gid] = _I8_TO_STRAND[int(iqidx.gene_strand[g_idx])]
            gi.sources[gid] = iqidx.get_str(int(iqidx.gene_source_str[g_idx]))
            g_attr = iqidx.get_str(int(iqidx.gene_attr_str[g_idx]))
            if g_attr:
                gi.feature_attributes[gid] = g_attr
            gi.gene_regions[gid] = (int(iqidx.gene_start[g_idx]), int(iqidx.gene_end[g_idx]))

            for t_idx in range(int(g_tx_ptr[g_idx]), int(g_tx_ptr[g_idx + 1])):
                tid = iqidx.get_str(int(iqidx.tx_id_str[t_idx]))
                gi.isoform_strands[tid] = _I8_TO_STRAND[int(iqidx.tx_strand[t_idx])]
                gi.sources[tid] = iqidx.get_str(int(iqidx.tx_source_str[t_idx]))
                gi.gene_id_map[tid] = gid
                t_attr = iqidx.get_str(int(iqidx.tx_attr_str[t_idx]))
                if t_attr:
                    gi.feature_attributes[tid] = t_attr

                ex_0, ex_1 = int(t_ex_ptr[t_idx]), int(t_ex_ptr[t_idx + 1])
                exons = [(int(e_start[i]), int(e_end[i])) for i in range(ex_0, ex_1)]
                if exons:
                    gi.all_isoforms_exons[tid] = exons
                    gi.all_isoforms_introns[tid] = junctions_from_blocks(exons)

                ot_0, ot_1 = int(t_ot_ptr[t_idx]), int(t_ot_ptr[t_idx + 1])
                if ot_1 > ot_0:
                    gi.other_features[tid] = [
                        (int(o_start[i]), int(o_end[i]), OTHER_FEATURE_NAMES[int(o_type[i])])
                        for i in range(ot_0, ot_1)
                    ]

        gi._finalize_profiles_and_properties(prepare_profiles)
        return gi

    def _finalize_profiles_and_properties(self, prepare_profiles: bool) -> None:
        introns = set()
        exons = set()
        for t_introns in self.all_isoforms_introns.values():
            introns.update(t_introns)
        for t_exons in self.all_isoforms_exons.values():
            exons.update(t_exons)
        self.intron_profiles.set_features(sorted(introns))
        self.exon_profiles.set_features(sorted(exons))

        if prepare_profiles:
            self.split_exon_profiles.set_features(self.split_exons(self.exon_profiles.features))
            self.set_junction_profiles(self.all_isoforms_introns, self.all_isoforms_exons)
            self.exon_property_map = self.set_feature_properties(self.all_isoforms_exons, self.exon_profiles)
            self.intron_property_map = self.set_feature_properties(self.all_isoforms_introns, self.intron_profiles)

    @classmethod
    def from_region(cls, chr_id: str, start: int, end: int, delta: int = 0) -> "GeneInfo":
        gi = cls.__new__(cls)
        gi.db = None
        gi.gene_db_list = []
        gi.chr_id = chr_id
        gi.start = start
        gi.end = end
        gi.delta = delta
        gi.all_read_region_start = start
        gi.all_read_region_end = end
        gi.canonical_sites = {}
        gi.gene_regions = {}
        gi.reference_region = None
        gi.intron_profiles = FeatureProfiles()
        gi.exon_profiles = FeatureProfiles()
        gi.split_exon_profiles = FeatureProfiles()
        gi.all_isoforms_introns = {}
        gi.all_isoforms_exons = {}
        gi.isoform_strands = {}
        gi.gene_strands = {}
        gi.other_features = {}
        gi.sources = {}
        gi.gene_id_map = {}
        gi.feature_attributes = {}
        gi.exon_property_map = None
        gi.intron_property_map = None
        return gi

    @classmethod
    def from_model(cls, transcript_model: TranscriptModel, delta: int = 0) -> "GeneInfo":
        return cls.from_models([transcript_model], delta)

    @classmethod
    def from_models(cls, transcript_model_storage: Sequence[TranscriptModel], delta: int = 0) -> "GeneInfo":
        if not transcript_model_storage:
            return cls.from_region("", 0, 0, delta)
        gi = cls.from_region(
            transcript_model_storage[0].chr_id,
            transcript_model_storage[0].get_start(),
            transcript_model_storage[0].get_end(),
            delta,
        )
        introns = set()
        exons = set()
        for tm in transcript_model_storage:
            gi.start = min(gi.start, tm.get_start())
            gi.end = max(gi.end, tm.get_end())
            tid = tm.transcript_id
            gi.all_isoforms_exons[tid] = tm.exon_blocks
            gi.all_isoforms_introns[tid] = junctions_from_blocks(tm.exon_blocks)
            exons.update(gi.all_isoforms_exons[tid])
            introns.update(gi.all_isoforms_introns[tid])
            gi.isoform_strands[tid] = tm.strand
            gi.gene_strands[tm.gene_id] = tm.strand
            gi.sources[tid] = tm.source
            gi.sources[tm.gene_id] = tm.source
            gi.other_features[tid] = tm.other_features
            gi.gene_id_map[tid] = tm.gene_id

        gi.all_read_region_start = gi.start
        gi.all_read_region_end = gi.end
        gi.intron_profiles.set_features(sorted(introns))
        gi.exon_profiles.set_features(sorted(exons))
        gi.split_exon_profiles.set_features(gi.split_exons(gi.exon_profiles.features))
        gi.set_junction_profiles(gi.all_isoforms_introns, gi.all_isoforms_exons)
        return gi

    def empty(self) -> bool:
        return len(self.all_isoforms_exons) == 0

    def get_gene_region(self) -> Tuple[str, int, int]:
        chrom = self.gene_db_list[0].seqid
        return chrom, min(g.start for g in self.gene_db_list), max(g.end for g in self.gene_db_list)

    def get_gene_regions(self) -> Dict[str, Tuple[int, int]]:
        if not self.gene_regions and self.gene_db_list:
            self.gene_regions = {g.id: (g.start, g.end) for g in self.gene_db_list}
        return self.gene_regions

    def transcript_start(self, transcript_id: str) -> int:
        return self.all_isoforms_exons[transcript_id][0][0]

    def transcript_end(self, transcript_id: str) -> int:
        return self.all_isoforms_exons[transcript_id][-1][1]

    def transcript_region(self, transcript_id: str) -> Tuple[int, int]:
        exons = self.all_isoforms_exons[transcript_id]
        return (exons[0][0], exons[-1][1])

    def transcript_exon_count(self, transcript_id: str) -> int:
        return len(self.all_isoforms_exons[transcript_id])

    def total_transcript_length(self, transcript_id: str) -> int:
        return sum(interval_len(e) for e in self.all_isoforms_exons[transcript_id])

    def set_junction_profiles(self, all_isoforms_introns, all_isoforms_exons) -> None:
        eq_exact = partial(equal_ranges, delta=0)
        for t_id in all_isoforms_exons:
            t_region = self.transcript_region(t_id)
            self.intron_profiles.set_profiles(t_id, all_isoforms_introns[t_id], t_region, eq_exact)
            self.exon_profiles.set_profiles(t_id, all_isoforms_exons[t_id], t_region, eq_exact)
            self.split_exon_profiles.set_profiles(t_id, all_isoforms_exons[t_id], t_region, contains)

    @staticmethod
    def split_exons(exons: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
        exon_starts = sorted(x[0] for x in exons)
        exon_ends = sorted(x[1] for x in exons)
        current_state = 0
        starts_pos = 0
        ends_pos = 0
        last_border = -1
        exon_blocks: List[Tuple[int, int]] = []

        while starts_pos < len(exon_starts):
            if exon_starts[starts_pos] <= exon_ends[ends_pos]:
                if starts_pos == 0 or exon_starts[starts_pos] > exon_starts[starts_pos - 1]:
                    cur_border = exon_starts[starts_pos]
                    if last_border != -1 and current_state > 0:
                        exon_blocks.append((last_border, cur_border - 1))
                    last_border = cur_border
                current_state += 1
                starts_pos += 1
            else:
                if ends_pos == 0 or exon_ends[ends_pos] > exon_ends[ends_pos - 1]:
                    cur_border = exon_ends[ends_pos]
                    exon_blocks.append((last_border, cur_border))
                    last_border = cur_border + 1
                current_state -= 1
                ends_pos += 1

        while ends_pos < len(exon_ends):
            if ends_pos == 0 or exon_ends[ends_pos] > exon_ends[ends_pos - 1]:
                cur_border = exon_ends[ends_pos]
                exon_blocks.append((last_border, cur_border))
                last_border = cur_border + 1
            current_state -= 1
            ends_pos += 1

        return exon_blocks

    def set_feature_properties(self, isoforms_to_feature_map, feature_profiles: FeatureProfiles) -> List[FeatureInfo]:
        """Exact O(N + K) active sweep-line replacing O(N^2) nested pair comparison."""
        feats = feature_profiles.features
        similar_features = set()
        contained_features = set()
        active: List[int] = []
        delta = self.delta

        for idx, f1 in enumerate(feats):
            cutoff = f1[0] - delta
            next_active: List[int] = []
            for j in active:
                f2 = feats[j]
                if f2[1] < cutoff:
                    continue
                next_active.append(j)
                if f1 == f2:
                    continue
                if equal_ranges(f1, f2, delta):
                    similar_features.add(f1)
                    similar_features.add(f2)
                if contains(f1, f2):
                    contained_features.add(f2)
                elif contains(f2, f1):
                    contained_features.add(f1)
            next_active.append(idx)
            active = next_active

        feature_to_isoform: Dict[Tuple[int, int], list] = defaultdict(list)
        for t, iso_feats in isoforms_to_feature_map.items():
            if not iso_feats:
                continue
            if len(iso_feats) == 1:
                feature_to_isoform[iso_feats[0]].append((t, "T"))
            else:
                feature_to_isoform[iso_feats[0]].append((t, "T"))
                feature_to_isoform[iso_feats[-1]].append((t, "T"))
                for e in iso_feats[1:-1]:
                    feature_to_isoform[e].append((t, ""))

        feature_properties: List[FeatureInfo] = []
        for feature in feats:
            linked = feature_to_isoform[feature]
            if all(x[1] == "T" for x in linked):
                ftype = "X"
            elif any(x[1] == "T" for x in linked):
                ftype = "T"
            else:
                ftype = "I"

            if feature in similar_features:
                ftype += "S"
            if feature in contained_features:
                ftype += "C"

            strands = sorted({self.isoform_strands[t[0]] for t in linked})
            gene_ids = {self.gene_id_map[t[0]] for t in linked}
            if len(linked) == 1:
                ftype += "U"
            elif len(gene_ids) > 1:
                ftype += "M"

            feature_properties.append(
                FeatureInfo(self.chr_id, feature[0], feature[1], "".join(strands), ftype, list(gene_ids))
            )
        return feature_properties

    @staticmethod
    def _group_overlapping_intervals(intervals: Sequence[Tuple[int, int]]) -> List[list]:
        groups: List[list] = []
        current: Optional[list] = None
        for pos, (s, e) in enumerate(intervals):
            if current is None or s > current[1]:
                current = [s, e, [pos]]
                groups.append(current)
            else:
                current[2].append(pos)
                if e > current[1]:
                    current[1] = e
        return groups

    @property
    def exon_overlap_regions(self) -> List[ExonRegion]:
        regions = getattr(self, "_exon_overlap_regions", None)
        if regions is None:
            exon_property_map = getattr(self, "exon_property_map", None)
            exon_profiles = getattr(self, "exon_profiles", None)
            if not exon_property_map or exon_profiles is None or not exon_profiles.features:
                regions = []
            else:
                regions, _ = self.build_exon_overlap_regions(exon_profiles.features, exon_property_map)
            self._exon_overlap_regions = regions
        return regions

    def build_exon_overlap_regions(self, exon_features, exon_property_map):
        n = len(exon_features)
        region_map = [-1] * n
        regions: List[ExonRegion] = []
        by_strand = defaultdict(list)
        for i in range(n):
            by_strand[exon_property_map[i].strand].append(i)

        for strand, idx_list in by_strand.items():
            intervals = [exon_features[i] for i in idx_list]
            for region_start, region_end, positions in self._group_overlapping_intervals(intervals):
                region = ExonRegion(self.chr_id, region_start, region_end, strand)
                region_index = len(regions)
                regions.append(region)
                for pos in positions:
                    i = idx_list[pos]
                    region.member_exon_indices.append(i)
                    region.gene_ids.update(exon_property_map[i].gene_ids)
                    region_map[i] = region_index
        return regions, region_map

    def build_exon_splice_site_regions(self) -> List[ExonSpliceSiteRegion]:
        regions: List[ExonSpliceSiteRegion] = []
        if not self.all_isoforms_exons:
            return regions

        gene_to_tx = defaultdict(list)
        for t_id, gene_id in self.gene_id_map.items():
            gene_to_tx[gene_id].append(t_id)

        for gene_id, t_ids in gene_to_tx.items():
            strand = self.gene_strands.get(gene_id, ".")
            exon_set = set()
            intron_set = set()
            for t_id in t_ids:
                exon_set.update(self.all_isoforms_exons.get(t_id, ()))
                intron_set.update(self.all_isoforms_introns.get(t_id, ()))
            if not exon_set:
                continue
            introns = sorted(intron_set)
            candidates = [
                e for e in sorted(exon_set)
                if not any(e[0] <= i[0] and e[1] >= i[1] for i in introns)
            ]
            if not candidates:
                continue
            for region_start, region_end, positions in self._group_overlapping_intervals(candidates):
                region = ExonSpliceSiteRegion(self.chr_id, region_start, region_end, strand, gene_id)
                region.candidates = [candidates[p] for p in positions]
                regions.append(region)

        for region in regions:
            left_counts: Dict[int, int] = defaultdict(int)
            right_counts: Dict[int, int] = defaultdict(int)
            for s, e in region.candidates:
                left_counts[s] += 1
                right_counts[e] += 1
            region.left_unique = [left_counts[s] == 1 for s, _ in region.candidates]
            region.right_unique = [right_counts[e] == 1 for _, e in region.candidates]
        return regions

    def set_reference_sequence(self, start: int, end: int, chr_record) -> None:
        self.all_read_region_start = start
        self.all_read_region_end = end
        self.reference_region = str(chr_record[start - 1 : end])
        self.canonical_sites = {}

    def get_ref_seq(self, start: int, end: int) -> str:
        if self.reference_region is None:
            return ""
        return self.reference_region[start - self.all_read_region_start : end - self.all_read_region_start + 1]

    def serialize(self, outfile) -> None:
        write_int(self.delta, outfile)
        gene_ids = [g.id for g in self.gene_db_list] if self.gene_db_list else list(self.gene_strands.keys())
        write_int(len(gene_ids), outfile)
        for gid in gene_ids:
            write_string(gid, outfile)
        write_string(self.chr_id, outfile)
        write_int(self.start, outfile)
        write_int(self.end, outfile)

    @classmethod
    def deserialize(cls, infile, genedb) -> "GeneInfo":
        delta = read_int(infile)
        gene_count = read_int(infile)
        gene_ids = [read_string(infile) for _ in range(gene_count)]
        chr_id = read_string(infile)
        start = read_int(infile)
        end = read_int(infile)
        if not gene_ids or genedb is None:
            gi = cls.from_region(chr_id, start, end, delta)
            gi.db = genedb
            return gi
        if isinstance(genedb, MmapAnnotationIndex):
            idx_map = getattr(genedb, "_gene_id_to_idx_cache", None)
            if idx_map is None:
                idx_map = {genedb.get_str(int(s)): idx for idx, s in enumerate(genedb.gene_id_str)}
                genedb._gene_id_to_idx_cache = idx_map
            gene_indices = [idx_map[gid] for gid in gene_ids if gid in idx_map]
            if not gene_indices:
                gi = cls.from_region(chr_id, start, end, delta)
                gi.db = genedb
                return gi
            return cls.from_iqidx_cluster(genedb, gene_indices, delta=delta, prepare_profiles=True)
        gene_dbs = [genedb[gid] for gid in gene_ids]
        return cls(gene_dbs, genedb, delta=delta, prepare_profiles=True)
