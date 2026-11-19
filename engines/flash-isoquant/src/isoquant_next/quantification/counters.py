"""Core gene, transcript, exon, intron, and intron-retention quantification counters."""

from __future__ import annotations

import logging
from collections import OrderedDict, defaultdict
from typing import Iterable, Optional

import numpy as np

from isoquant_next.alignment.read_groups import AbstractReadGrouper
from isoquant_next.annotation.gene_info import FeatureInfo
from isoquant_next.assignment.types import MatchEventSubtype, ReadAssignmentType
from isoquant_next.config.enums import GroupedOutputFormat, NormalizationMethod
from isoquant_next.core.file_naming import (
    counts_file_name,
    counts_prefix,
    counts_stats_file_name,
    counts_usable_file_name,
    tpm_file_name,
    tpm_prefix,
)
from isoquant_next.quantification.mtx_converter import (
    GROUP_COUNT_CUTOFF,
    convert_profile_to_matrix,
    convert_profile_to_mtx,
    convert_to_matrix,
    convert_to_mtx,
)
from isoquant_next.quantification.weights import (
    GeneAssignmentExtractor,
    ReadWeightCounter,
    TranscriptAssignmentExtractor,
    convert_ungrouped_to_tpm,
)

logger = logging.getLogger("IsoQuant")

INTRON_RETENTION_EVENT_TYPES = frozenset(
    {
        MatchEventSubtype.intron_retention,
        MatchEventSubtype.unspliced_intron_retention,
    }
)


class IncrementalDict:
    """Sparse group counter avoiding eager zero-entry materialization on `.get()`."""

    __slots__ = ("data", "default_type")

    def __init__(self, default_type=float) -> None:
        self.data: dict = {}
        self.default_type = default_type

    def inc(self, key: int, value=1.0) -> None:
        if key not in self.data:
            self.data[key] = self.default_type(value)
        else:
            self.data[key] += self.default_type(value)

    def get(self, key: int):
        if key not in self.data:
            return self.default_type()
        return self.data[key]


class AbstractCounter:
    """Base interface for IsoQuant per-chromosome and sample-merged feature counters."""

    def __init__(self, output_prefix: str, ignore_read_groups: bool = False) -> None:
        self.ignore_read_groups = ignore_read_groups
        self.output_prefix = output_prefix
        self.output_counts_prefix = counts_prefix(output_prefix)
        self.output_tpm_prefix = tpm_prefix(output_prefix)
        self.output_counts_file_name = counts_file_name(
            self.output_counts_prefix, linear=not ignore_read_groups
        )
        self.output_tpm_file_name = tpm_file_name(self.output_tpm_prefix)
        self.output_file = self.output_counts_file_name
        with open(self.output_file, "w", encoding="utf-8"):
            pass
        self.output_stats_file_name: Optional[str] = None
        self.usable_file_name: Optional[str] = None

    def get_output_file_handler(self):
        return open(self.output_file, "a", encoding="utf-8")

    def add_read_info(self, read_assignment) -> None:
        raise NotImplementedError()

    def add_confirmed_features(self, features) -> None:
        raise NotImplementedError()

    def dump(self) -> None:
        raise NotImplementedError()

    def flush(self) -> None:
        pass

    def finalize(self, args=None) -> None:
        raise NotImplementedError()

    def add_unassigned(self, read_assignment) -> None:
        raise NotImplementedError()

    def add_unaligned(self, n_reads: int = 1) -> None:
        raise NotImplementedError()


class CompositeCounter:
    """Fan-out delegator coordinating multiple feature counters in one pass."""

    def __init__(self, counters: Optional[list] = None) -> None:
        self.counters: list = list(counters) if counters is not None else []

    def add_counters(self, counters: Iterable) -> None:
        self.counters.extend(counters)

    def add_counter(self, counter) -> None:
        self.counters.append(counter)

    def add_read_info(self, read_assignment) -> None:
        for p in self.counters:
            p.add_read_info(read_assignment)

    def add_confirmed_features(self, features) -> None:
        for p in self.counters:
            p.add_confirmed_features(features)

    def dump(self) -> None:
        for p in self.counters:
            p.dump()

    def flush(self) -> None:
        for p in self.counters:
            p.flush()

    def finalize(self, args=None) -> None:
        for p in self.counters:
            p.finalize(args)

    def add_unassigned(self, read_assignment) -> None:
        for p in self.counters:
            p.add_unassigned(read_assignment)

    def add_unaligned(self, n_reads: int = 1) -> None:
        for p in self.counters:
            p.add_unaligned(n_reads)


class AssignedFeatureCounter(AbstractCounter):
    """Gene or transcript abundance counter using strategy weights and StringPool group IDs."""

    def __init__(
        self,
        output_prefix: str,
        assignment_extractor,
        string_pools,
        read_counter: ReadWeightCounter,
        all_features=None,
        group_index: int = 0,
    ) -> None:
        super().__init__(output_prefix, string_pools is None)
        self.assignment_extractor = assignment_extractor
        self.string_pools = string_pools
        self.all_features: set[str] = set(all_features) if all_features is not None else set()
        self.group_index = group_index
        self.read_counter = read_counter

        self.ambiguous_reads = 0
        self.reads_for_tpm: dict[int, int] = defaultdict(int)
        self.not_assigned_reads = 0
        self.not_aligned_reads = 0
        self.feature_counter: dict[str, IncrementalDict] = defaultdict(IncrementalDict)
        self.confirmed_features: set[str] = set()
        self.output_stats_file_name = counts_stats_file_name(self.output_counts_file_name)
        self.usable_file_name = counts_usable_file_name(self.output_counts_file_name)

    def _resolve_group_id(self, read_assignment) -> int:
        if self.ignore_read_groups or not read_assignment.read_group_ids:
            return 0
        return read_assignment.read_group_ids[self.group_index]

    def add_read_info(self, read_assignment=None) -> None:
        if not read_assignment:
            self.not_aligned_reads += 1
            return

        group_id = self._resolve_group_id(read_assignment)
        assignment_type = self.assignment_extractor.get_assignment_type(read_assignment)
        if assignment_type.is_unassigned() or not read_assignment.isoform_matches:
            self.not_assigned_reads += 1
            self.reads_for_tpm[group_id] += 1
            return

        feature_ids = self.assignment_extractor.get_features(read_assignment)
        if not feature_ids:
            self.not_assigned_reads += 1
            self.reads_for_tpm[group_id] += 1
            return

        self.reads_for_tpm[group_id] += 1
        if assignment_type == ReadAssignmentType.ambiguous:
            count_value = self.read_counter.process_ambiguous(len(feature_ids))
            for feature_id in feature_ids:
                self.feature_counter[feature_id].inc(group_id, count_value)
                if count_value > 0.0:
                    self.all_features.add(feature_id)
            self.ambiguous_reads += 1
        elif assignment_type.is_inconsistent():
            count_value = self.read_counter.process_inconsistent(assignment_type, len(feature_ids))
            if count_value > 0.0:
                for feature_id in feature_ids:
                    self.feature_counter[feature_id].inc(group_id, count_value)
                    self.all_features.add(feature_id)
        elif assignment_type.is_unique():
            feature_id = next(iter(feature_ids))
            count_value = self.read_counter.process_unique(read_assignment)
            if count_value > 0.0:
                self.feature_counter[feature_id].inc(group_id, count_value)
                self.all_features.add(feature_id)
                if self.assignment_extractor.confirms_feature(read_assignment):
                    self.confirmed_features.add(feature_id)

    def add_unassigned(self, read_assignment) -> None:
        group_id = self._resolve_group_id(read_assignment)
        self.not_assigned_reads += 1
        self.reads_for_tpm[group_id] += 1

    def add_unaligned(self, n_reads: int = 1) -> None:
        self.not_aligned_reads += n_reads

    def format_header(self, all_groups: list[str], value_name: str = "count") -> str:
        if self.ignore_read_groups:
            return f"feature_id\t{value_name}\n"
        return "feature_id\t" + "\t".join(all_groups) + "\n"

    def add_confirmed_features(self, features: Iterable[str]) -> None:
        self.confirmed_features.update(features)

    def dump(self) -> None:
        all_features = sorted(x for x in self.all_features if x is not None)
        for feature_id in all_features:
            if feature_id in self.confirmed_features:
                continue
            f_data = self.feature_counter[feature_id].data
            for g in f_data:
                f_data[g] = 0.0

        if self.ignore_read_groups:
            self.dump_ungrouped(all_features)
        else:
            self.dump_grouped(all_features)
        self.dump_usable()

    def _get_ordered_groups(self) -> list[str]:
        if self.string_pools is None:
            return [AbstractReadGrouper.default_group_id]
        pool = self.string_pools.get_read_group_pool(self.group_index)
        return [pool.get_str(i) for i in range(len(pool))]

    def _get_group_name(self, group_id: int) -> str:
        if self.string_pools is None:
            return AbstractReadGrouper.default_group_id
        return self.string_pools.get_read_group_pool(self.group_index).get_str(group_id)

    def _get_num_groups(self) -> int:
        if self.string_pools is None:
            return 1
        return len(self.string_pools.get_read_group_pool(self.group_index))

    def dump_ungrouped(self, all_features: list[str]) -> None:
        with self.get_output_file_handler() as output_file:
            output_file.write(self.format_header(self._get_ordered_groups()))
            lines = [f"{fid}\t{self.feature_counter[fid].get(0):.2f}\n" for fid in all_features]
            output_file.writelines(lines)

        with open(self.output_stats_file_name, "w", encoding="utf-8") as f:
            f.write(f"__ambiguous\t{self.ambiguous_reads}\n")
            f.write(f"__no_feature\t{self.not_assigned_reads}\n")
            f.write(f"__not_aligned\t{self.not_aligned_reads}\n")

    def dump_grouped(self, all_features: list[str]) -> None:
        with self.get_output_file_handler() as output_file:
            output_file.write("feature_id\tgroup_id\tcount\n")
            for feature_id in all_features:
                for group_id, count in self.feature_counter[feature_id].data.items():
                    if count != 0.0:
                        group_name = self._get_group_name(group_id)
                        output_file.write(f"{feature_id}\t{group_name}\t{count:.2f}\n")

    def dump_usable(self) -> None:
        with open(self.usable_file_name, "w", encoding="utf-8") as f:
            for group_id in range(self._get_num_groups()):
                f.write(f"{self._get_group_name(group_id)}\t{self.reads_for_tpm[group_id]}\n")

    def load_usable(self, usable_file_name: str) -> None:
        with open(usable_file_name, "r", encoding="utf-8") as fh:
            if self.string_pools is None:
                for line in fh:
                    v = line.strip().split("\t")
                    if len(v) >= 2:
                        self.reads_for_tpm[0] += int(v[1])
            else:
                pool = self.string_pools.get_read_group_pool(self.group_index)
                for line in fh:
                    v = line.strip().split("\t")
                    if len(v) >= 2 and v[0] in pool:
                        self.reads_for_tpm[pool.get_int(v[0])] += int(v[1])

    def finalize(self, args=None) -> None:
        if not args or not getattr(args, "normalization_method", None):
            normalization_method = NormalizationMethod.none
        elif args.normalization_method not in NormalizationMethod.__members__:
            logger.warning(
                "%s in not a valid normalization method, will use simple normalization",
                args.normalization_method,
            )
            normalization_method = NormalizationMethod.simple
        else:
            normalization_method = NormalizationMethod[args.normalization_method]

        if self.ignore_read_groups:
            if normalization_method == NormalizationMethod.none:
                return
            convert_ungrouped_to_tpm(
                self.output_file,
                self.output_tpm_file_name,
                normalization_method,
                self.reads_for_tpm[0],
            )
            return

        num_groups = self._get_num_groups()
        counts_format = {GroupedOutputFormat[f] for f in args.counts_format}
        explicit_matrix = GroupedOutputFormat.matrix in counts_format
        default_format = GroupedOutputFormat.default in counts_format
        actual_num_groups = num_groups

        if explicit_matrix or default_format:
            max_groups = GROUP_COUNT_CUTOFF if (default_format and not explicit_matrix) else 0
            if self.string_pools is not None and max_groups > 0 and num_groups > max_groups:
                logger.info(
                    "Skipping full matrix conversion: %d groups exceeds limit of %d. "
                    "Use '--counts_format matrix' to force conversion.",
                    num_groups,
                    max_groups,
                )
            else:
                actual_num_groups = (
                    convert_to_matrix(self.output_file, self.output_counts_prefix, max_groups=max_groups)
                    or num_groups
                )
                if normalization_method != NormalizationMethod.none:
                    reads_for_tpm = None
                    if normalization_method == NormalizationMethod.usable_reads:
                        reads_for_tpm = {
                            self._get_group_name(gid): self.reads_for_tpm[gid] for gid in range(num_groups)
                        }
                    convert_to_matrix(
                        self.output_file,
                        self.output_tpm_prefix,
                        convert_to_tpm=True,
                        usable_reads_per_group=reads_for_tpm,
                        max_groups=max_groups,
                    )

        if GroupedOutputFormat.mtx in counts_format or (
            default_format and actual_num_groups > GROUP_COUNT_CUTOFF
        ):
            convert_to_mtx(self.output_file, self.output_counts_prefix)
            if normalization_method != NormalizationMethod.none:
                reads_for_tpm = None
                if normalization_method == NormalizationMethod.usable_reads:
                    reads_for_tpm = {
                        self._get_group_name(gid): self.reads_for_tpm[gid] for gid in range(num_groups)
                    }
                convert_to_mtx(
                    self.output_file,
                    self.output_tpm_prefix,
                    convert_to_tpm=True,
                    usable_reads_per_group=reads_for_tpm,
                )


def create_gene_counter(
    output_file_name: str,
    strategy: str,
    complete_feature_list=None,
    string_pools=None,
    group_index: int = 0,
) -> AssignedFeatureCounter:
    return AssignedFeatureCounter(
        output_file_name,
        GeneAssignmentExtractor,
        string_pools,
        ReadWeightCounter(strategy),
        complete_feature_list,
        group_index,
    )


def create_transcript_counter(
    output_file_name: str,
    strategy: str,
    complete_feature_list=None,
    string_pools=None,
    group_index: int = 0,
) -> AssignedFeatureCounter:
    return AssignedFeatureCounter(
        output_file_name,
        TranscriptAssignmentExtractor,
        string_pools,
        ReadWeightCounter(strategy),
        complete_feature_list,
        group_index,
    )


class ProfileFeatureCounter(AbstractCounter):
    """Include/exclude feature counter with sparse nonzero profile traversal."""

    def __init__(self, output_prefix: str, string_pools=None, group_index: int = 0) -> None:
        super().__init__(output_prefix, string_pools is None)
        self.string_pools = string_pools
        self.group_index = group_index
        self.inclusion_feature_counter: dict = defaultdict(lambda: IncrementalDict(int))
        self.exclusion_feature_counter: dict = defaultdict(lambda: IncrementalDict(int))
        self.feature_name_dict: OrderedDict = OrderedDict()
        self.encountered_group_ids: set[int] = set()

    def add_read_info_from_profile(
        self, gene_feature_profile, read_assigned_strand: str, feature_property_map, group_id: int = 0
    ) -> None:
        self.encountered_group_ids.add(group_id)
        # Skip leading/trailing zero spans quickly when Int8ProfileView / NumPy array is supplied
        raw_arr = getattr(gene_feature_profile, "array", None)
        arr = raw_arr if raw_arr is not None else np.asarray(gene_feature_profile, dtype=np.int8)
        nz_indices = np.flatnonzero(arr)
        for idx in nz_indices:
            i = int(idx)
            prop = feature_property_map[i]
            if read_assigned_strand not in prop.strand:
                continue
            val = int(arr[i])
            feature_id = prop.id
            if val == 1:
                self.inclusion_feature_counter[feature_id].inc(group_id)
            elif val == -1:
                self.exclusion_feature_counter[feature_id].inc(group_id)
            else:
                continue
            if feature_id not in self.feature_name_dict:
                self.feature_name_dict[feature_id] = prop.to_str()

    def _get_group_name(self, group_id: int) -> str:
        if self.string_pools is None:
            return AbstractReadGrouper.default_group_id
        return self.string_pools.get_read_group_pool(self.group_index).get_str(group_id)

    def dump(self) -> None:
        with open(self.output_counts_file_name, "w", encoding="utf-8") as f:
            f.write(FeatureInfo.header() + "\tgroup_id\tinclude_counts\texclude_counts\n")
            all_group_ids = sorted(self.encountered_group_ids)
            for feature_id, feature_name in self.feature_name_dict.items():
                inc_dict = self.inclusion_feature_counter[feature_id]
                exc_dict = self.exclusion_feature_counter[feature_id]
                for group_id in all_group_ids:
                    incl_count = inc_dict.get(group_id)
                    excl_count = exc_dict.get(group_id)
                    if incl_count > 0 or excl_count > 0:
                        f.write(
                            f"{feature_name}\t{self._get_group_name(group_id)}\t{int(incl_count)}\t{int(excl_count)}\n"
                        )

    def _feature_type(self) -> str:
        return "feature"

    def finalize(self, args=None) -> None:
        if not args or self.ignore_read_groups:
            return
        counts_format = {GroupedOutputFormat[f] for f in args.counts_format}
        explicit_matrix = GroupedOutputFormat.matrix in counts_format
        default_format = GroupedOutputFormat.default in counts_format
        num_groups = (
            len(self.string_pools.get_read_group_pool(self.group_index))
            if self.string_pools is not None
            else 1
        )
        actual_num_groups = num_groups
        if explicit_matrix or default_format:
            max_groups = GROUP_COUNT_CUTOFF if (default_format and not explicit_matrix) else 0
            if self.string_pools is not None and max_groups > 0 and num_groups > max_groups:
                logger.info(
                    "Skipping full matrix conversion: %d groups exceeds limit of %d. "
                    "Use '--counts_format matrix' to force conversion.",
                    num_groups,
                    max_groups,
                )
            else:
                actual_num_groups = (
                    convert_profile_to_matrix(
                        self.output_file,
                        self.output_counts_prefix,
                        feature_type=self._feature_type(),
                        max_groups=max_groups,
                    )
                    or num_groups
                )
        if GroupedOutputFormat.mtx in counts_format or (
            default_format and actual_num_groups > GROUP_COUNT_CUTOFF
        ):
            convert_profile_to_mtx(self.output_file, self.output_counts_prefix)

    @staticmethod
    def is_valid(assignment) -> bool:
        return (
            assignment is not None
            and getattr(assignment, "exon_gene_profile", None) is not None
            and getattr(assignment, "intron_gene_profile", None) is not None
            and getattr(assignment, "gene_info", None) is not None
        )

    @staticmethod
    def is_assigned_to_gene(assignment) -> bool:
        return (
            not assignment.gene_assignment_type.is_unassigned()
            and not assignment.gene_assignment_type.is_ambiguous()
        )


class ExonCounter(ProfileFeatureCounter):
    def _feature_type(self) -> str:
        return "exon"

    def add_read_info(self, read_assignment) -> None:
        if not ProfileFeatureCounter.is_valid(read_assignment) or not ProfileFeatureCounter.is_assigned_to_gene(
            read_assignment
        ):
            return
        group_id = (
            0
            if (self.ignore_read_groups or not read_assignment.read_group_ids)
            else read_assignment.read_group_ids[self.group_index]
        )
        self.add_read_info_from_profile(
            read_assignment.exon_gene_profile,
            read_assignment.strand,
            read_assignment.gene_info.exon_property_map,
            group_id,
        )


class IntronCounter(ProfileFeatureCounter):
    def _feature_type(self) -> str:
        return "intron"

    def add_read_info(self, read_assignment) -> None:
        if not ProfileFeatureCounter.is_valid(read_assignment) or not ProfileFeatureCounter.is_assigned_to_gene(
            read_assignment
        ):
            return
        group_id = (
            0
            if (self.ignore_read_groups or not read_assignment.read_group_ids)
            else read_assignment.read_group_ids[self.group_index]
        )
        self.add_read_info_from_profile(
            read_assignment.intron_gene_profile,
            read_assignment.strand,
            read_assignment.gene_info.intron_property_map,
            group_id,
        )


class IntronRetentionCounter(ProfileFeatureCounter):
    def add_read_info(self, read_assignment) -> None:
        if not ProfileFeatureCounter.is_valid(read_assignment) or not ProfileFeatureCounter.is_assigned_to_gene(
            read_assignment
        ):
            return
        if not read_assignment.isoform_matches:
            return

        group_id = (
            0
            if (self.ignore_read_groups or not read_assignment.read_group_ids)
            else read_assignment.read_group_ids[self.group_index]
        )
        gene_info = read_assignment.gene_info
        intron_features = gene_info.intron_profiles.features
        property_map = gene_info.intron_property_map
        if not property_map:
            return

        feature_index = getattr(gene_info, "_intron_feature_index", None)
        if feature_index is None:
            feature_index = {coords: i for i, coords in enumerate(intron_features)}
            gene_info._intron_feature_index = feature_index

        for match in read_assignment.isoform_matches:
            if match.assigned_transcript is None:
                continue
            isoform_introns = gene_info.all_isoforms_introns.get(match.assigned_transcript, [])
            for event in match.match_subclassifications:
                if event.event_type not in INTRON_RETENTION_EVENT_TYPES:
                    continue
                for idx in range(event.isoform_region[0], event.isoform_region[1] + 1):
                    if idx >= len(isoform_introns):
                        continue
                    feature_idx = feature_index.get(isoform_introns[idx])
                    if feature_idx is None:
                        continue
                    prop = property_map[feature_idx]
                    feature_id = prop.id
                    self.inclusion_feature_counter[feature_id].inc(group_id)
                    self.encountered_group_ids.add(group_id)
                    if feature_id not in self.feature_name_dict:
                        self.feature_name_dict[feature_id] = prop.to_str()


__all__ = [
    "AbstractCounter",
    "AssignedFeatureCounter",
    "CompositeCounter",
    "ExonCounter",
    "INTRON_RETENTION_EVENT_TYPES",
    "IncrementalDict",
    "IntronCounter",
    "IntronRetentionCounter",
    "ProfileFeatureCounter",
    "create_gene_counter",
    "create_transcript_counter",
]
