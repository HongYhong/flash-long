"""Region-based exon inclusion/exclusion and splice-site quantification counters."""

from __future__ import annotations

from collections import OrderedDict, defaultdict
from typing import Tuple

from isoquant_next.alignment.read_groups import AbstractReadGrouper
from isoquant_next.assignment.types import MatchEventSubtype
from isoquant_next.core.utils import junctions_from_blocks
from isoquant_next.quantification.counters import (
    AbstractCounter,
    IncrementalDict,
    ProfileFeatureCounter,
)

_TSS_MATCH_LEFT = frozenset(
    {
        MatchEventSubtype.terminal_site_match_left,
        MatchEventSubtype.terminal_site_match_left_precise,
    }
)
_TSS_MATCH_RIGHT = frozenset(
    {
        MatchEventSubtype.terminal_site_match_right,
        MatchEventSubtype.terminal_site_match_right_precise,
    }
)

DEFAULT_EXCLUSION_MARGIN = 50


class JointExonCounter(AbstractCounter):
    """Groups overlapping exons into regions and emits member inclusion or whole-region exclusion."""

    def __init__(self, output_prefix: str, string_pools=None, group_index: int = 0) -> None:
        super().__init__(output_prefix, string_pools is None)
        self.string_pools = string_pools
        self.group_index = group_index
        self.inclusion_counter = defaultdict(lambda: IncrementalDict(int))
        self.exclusion_counter = defaultdict(lambda: IncrementalDict(int))
        self.feature_row_prefix: OrderedDict = OrderedDict()
        self.encountered_group_ids: set[int] = set()

    def add_read_info(self, read_assignment) -> None:
        if not ProfileFeatureCounter.is_valid(read_assignment) or not ProfileFeatureCounter.is_assigned_to_gene(
            read_assignment
        ):
            return
        gene_info = read_assignment.gene_info
        if not getattr(gene_info, "exon_overlap_regions", None):
            return
        read_gene = None
        for m in read_assignment.isoform_matches:
            if m.assigned_gene:
                read_gene = m.assigned_gene
                break
        if read_gene is None:
            return
        group_id = (
            0
            if (self.ignore_read_groups or not read_assignment.read_group_ids)
            else read_assignment.read_group_ids[self.group_index]
        )
        self._classify(read_assignment.exon_gene_profile, read_assignment.strand, gene_info, group_id, read_gene)

    def _classify(self, profile, read_strand: str, gene_info, group_id: int, read_gene: str) -> None:
        self.encountered_group_ids.add(group_id)
        exon_features = gene_info.exon_profiles.features
        property_map = gene_info.exon_property_map
        for region in gene_info.exon_overlap_regions:
            members = region.member_exon_indices
            if not members or read_strand not in property_map[members[0]].strand:
                continue
            states = [int(profile[i]) for i in members]
            if any(s == 0 for s in states):
                continue
            inclusions = [i for i, s in zip(members, states) if s == 1]
            if inclusions:
                for i in inclusions:
                    if read_gene not in property_map[i].gene_ids:
                        continue
                    counter_key = (property_map[i].id, read_gene)
                    self.inclusion_counter[counter_key].inc(group_id)
                    row_key = ("inc", counter_key)
                    if row_key not in self.feature_row_prefix:
                        exon_start, exon_end = exon_features[i]
                        self.feature_row_prefix[row_key] = (
                            f"{region.chr_id}\t{region.start}\t{region.end}\t{region.strand}"
                            f"\t{exon_start}\t{exon_end}\t{read_gene}\tinclusion"
                        )
            else:
                if read_gene not in region.gene_ids:
                    continue
                counter_key = (region.id, read_gene)
                self.exclusion_counter[counter_key].inc(group_id)
                row_key = ("exc", counter_key)
                if row_key not in self.feature_row_prefix:
                    self.feature_row_prefix[row_key] = (
                        f"{region.chr_id}\t{region.start}\t{region.end}\t{region.strand}"
                        f"\t.\t.\t{read_gene}\texclusion"
                    )

    def _get_group_name(self, group_id: int) -> str:
        if self.string_pools is None:
            return AbstractReadGrouper.default_group_id
        return self.string_pools.get_read_group_pool(self.group_index).get_str(group_id)

    def dump(self) -> None:
        with open(self.output_counts_file_name, "w", encoding="utf-8") as f:
            f.write(
                "chr\tregion_start\tregion_end\tstrand\texon_start\texon_end\tgene_id\tfeature_kind\tgroup_id\tcount\n"
            )
            all_group_ids = sorted(self.encountered_group_ids)
            for (kind, fid), prefix in self.feature_row_prefix.items():
                counter = self.inclusion_counter if kind == "inc" else self.exclusion_counter
                for group_id in all_group_ids:
                    count = counter[fid].get(group_id)
                    if count > 0:
                        f.write(f"{prefix}\t{self._get_group_name(group_id)}\t{count}\n")

    def finalize(self, args=None) -> None:
        return

    def add_confirmed_features(self, features) -> None:
        return

    def add_unassigned(self, read_assignment) -> None:
        return

    def add_unaligned(self, n_reads: int = 1) -> None:
        return


class ExonSpliceSiteCounter(AbstractCounter):
    """Per-region full, left-half, right-half exon inclusion and region exclusion/ambiguous counter."""

    def __init__(
        self,
        output_prefix: str,
        string_pools=None,
        group_index: int = 0,
        delta: int = 6,
        exclusion_margin: int = DEFAULT_EXCLUSION_MARGIN,
        emit_read_ids: bool = False,
    ) -> None:
        super().__init__(output_prefix, string_pools is None)
        self.string_pools = string_pools
        self.group_index = group_index
        self.delta = delta
        self.exclusion_margin = exclusion_margin
        self.emit_read_ids = emit_read_ids
        self.regions: dict = {}

    @staticmethod
    def _new_bucket() -> dict:
        return {}

    def _bump(self, bucket: dict, group_id: int, read_id: str) -> None:
        entry = bucket.get(group_id)
        if entry is None:
            entry = [0, []]
            bucket[group_id] = entry
        entry[0] += 1
        if self.emit_read_ids:
            entry[1].append(read_id)

    def _get_region_record(self, region) -> dict:
        key = (region.chr_id, region.start, region.end, region.strand, region.gene_id)
        rec = self.regions.get(key)
        if rec is None:
            rec = {"cands": {}, "excl": self._new_bucket(), "amb": self._new_bucket()}
            self.regions[key] = rec
        return rec

    def _get_cand_bucket(self, rec: dict, cand: Tuple[int, int], side: str) -> dict:
        cand_rec = rec["cands"].get(cand)
        if cand_rec is None:
            cand_rec = {"full": self._new_bucket(), "left": self._new_bucket(), "right": self._new_bucket()}
            rec["cands"][cand] = cand_rec
        return cand_rec[side]

    def _get_anchoring(self, read_assignment, strand: str) -> Tuple[bool, bool]:
        polya_anchored = bool(read_assignment.polyA_found)
        tss_events = _TSS_MATCH_LEFT if strand == "+" else _TSS_MATCH_RIGHT
        tss_anchored = any(
            e.event_type in tss_events
            for m in read_assignment.isoform_matches
            for e in m.match_subclassifications
        )
        return tss_anchored, polya_anchored

    def add_read_info(self, read_assignment) -> None:
        if not ProfileFeatureCounter.is_valid(read_assignment) or not ProfileFeatureCounter.is_assigned_to_gene(
            read_assignment
        ):
            return
        gene_info = read_assignment.gene_info
        read_gene = None
        for m in read_assignment.isoform_matches:
            if m.assigned_gene:
                read_gene = m.assigned_gene
                break
        if read_gene is None:
            return

        regions = getattr(gene_info, "_exon_splice_site_regions", None)
        if regions is None:
            regions = gene_info.build_exon_splice_site_regions()
            gene_info._exon_splice_site_regions = regions
        if not regions:
            return

        blocks = read_assignment.corrected_exons if read_assignment.corrected_exons else read_assignment.exons
        if not blocks:
            return
        introns = junctions_from_blocks(blocks)
        strand = read_assignment.strand
        tss_anchored, polya_anchored = self._get_anchoring(read_assignment, strand)
        read_id = read_assignment.read_id
        group_id = (
            0
            if (self.ignore_read_groups or not read_assignment.read_group_ids)
            else read_assignment.read_group_ids[self.group_index]
        )

        for region in regions:
            if region.gene_id != read_gene or strand != region.strand:
                continue
            self._classify_region(
                region, blocks, introns, strand, tss_anchored, polya_anchored, group_id, read_id
            )

    def _classify_region(
        self, region, blocks, introns, strand: str, tss_anchored: bool, polya_anchored: bool, group_id: int, read_id: str
    ) -> None:
        delta = self.delta
        n = len(blocks)
        cands = region.candidates
        rstart, rend = region.start, region.end

        left_ext_anchored = tss_anchored if strand == "+" else polya_anchored
        right_ext_anchored = polya_anchored if strand == "+" else tss_anchored

        def left_demonstrated(j: int) -> bool:
            return j > 0 or left_ext_anchored

        def right_demonstrated(j: int) -> bool:
            return j < n - 1 or right_ext_anchored

        full_hits = []
        for xs, xe in cands:
            for j, b in enumerate(blocks):
                if abs(b[0] - xs) <= delta and abs(b[1] - xe) <= delta:
                    if ((j > 0) or (j < n - 1)) and left_demonstrated(j) and right_demonstrated(j):
                        full_hits.append((xs, xe))
                        break

        if len(full_hits) == 1:
            rec = self._get_region_record(region)
            self._bump(self._get_cand_bucket(rec, full_hits[0], "full"), group_id, read_id)
            return
        if len(full_hits) > 1:
            self._bump(self._get_region_record(region)["amb"], group_id, read_id)
            return

        margin = self.exclusion_margin
        for intr in introns:
            if intr[0] + margin <= rstart and intr[1] - margin >= rend:
                self._bump(self._get_region_record(region)["excl"], group_id, read_id)
                return

        if not any(b[0] <= rend and b[1] >= rstart for b in blocks):
            return

        half_hits = []
        for k, (xs, xe) in enumerate(cands):
            if region.left_unique[k]:
                for j, b in enumerate(blocks):
                    if abs(b[0] - xs) <= delta and abs(b[1] - xe) > delta and b[1] < xe and left_demonstrated(j):
                        half_hits.append(((xs, xe), "left"))
                        break
            if region.right_unique[k]:
                for j, b in enumerate(blocks):
                    if abs(b[1] - xe) <= delta and abs(b[0] - xs) > delta and b[0] > xs and right_demonstrated(j):
                        half_hits.append(((xs, xe), "right"))
                        break

        rec = self._get_region_record(region)
        if len(half_hits) == 1:
            cand, side = half_hits[0]
            self._bump(self._get_cand_bucket(rec, cand, side), group_id, read_id)
        else:
            self._bump(rec["amb"], group_id, read_id)

    def _get_group_name(self, group_id: int) -> str:
        if self.string_pools is None:
            return AbstractReadGrouper.default_group_id
        return self.string_pools.get_read_group_pool(self.group_index).get_str(group_id)

    def _header(self) -> str:
        cols = ["region_gene_candidate", "n_full", "n_left", "n_right", "group_id"]
        if self.emit_read_ids:
            cols += ["read_ids_full", "read_ids_left", "read_ids_right"]
        return "\t".join(cols) + "\n"

    @staticmethod
    def _fmt_list(items: list) -> str:
        return ";".join(items) if items else "NA"

    def _write_inclusion_rows(self, f, feature_id: str, full: dict, left: dict, right: dict) -> None:
        for gid in sorted(set(full) | set(left) | set(right)):
            fe, le, re_ = full.get(gid), left.get(gid), right.get(gid)
            fields = [
                feature_id,
                str(fe[0] if fe else 0),
                str(le[0] if le else 0),
                str(re_[0] if re_ else 0),
                self._get_group_name(gid),
            ]
            if self.emit_read_ids:
                fields += [
                    self._fmt_list(fe[1] if fe else []),
                    self._fmt_list(le[1] if le else []),
                    self._fmt_list(re_[1] if re_ else []),
                ]
            f.write("\t".join(fields) + "\n")

    def _write_region_rows(self, f, feature_id: str, bucket: dict) -> None:
        for gid in sorted(bucket):
            entry = bucket[gid]
            fields = [feature_id, str(entry[0]), "0", "0", self._get_group_name(gid)]
            if self.emit_read_ids:
                fields += [self._fmt_list(entry[1]), "NA", "NA"]
            f.write("\t".join(fields) + "\n")

    def dump(self) -> None:
        with open(self.output_counts_file_name, "w", encoding="utf-8") as f:
            f.write(self._header())
            for key in sorted(self.regions.keys(), key=lambda k: (k[0], k[1], k[2], k[4])):
                chr_id, rstart, rend, strand, gene_id = key
                rec = self.regions[key]
                region_coord = f"{chr_id}_{rstart}_{rend}_{strand}"
                for cand in sorted(rec["cands"].keys()):
                    buckets = rec["cands"][cand]
                    full, left, right = buckets["full"], buckets["left"], buckets["right"]
                    if not full and not left and not right:
                        continue
                    cand_str = f"{chr_id}_{cand[0]}_{cand[1]}_{strand}"
                    self._write_inclusion_rows(
                        f, f"{region_coord}__{gene_id}__{cand_str}", full, left, right
                    )
                if rec["excl"]:
                    self._write_region_rows(f, f"{region_coord}__{gene_id}__exclusion", rec["excl"])
                if rec["amb"]:
                    self._write_region_rows(f, f"{region_coord}__{gene_id}__ambiguous", rec["amb"])

    def finalize(self, args=None) -> None:
        return

    def add_confirmed_features(self, features) -> None:
        return

    def add_unassigned(self, read_assignment) -> None:
        return

    def add_unaligned(self, n_reads: int = 1) -> None:
        return


__all__ = ["DEFAULT_EXCLUSION_MARGIN", "ExonSpliceSiteCounter", "JointExonCounter"]
