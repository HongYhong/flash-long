"""Explicit per-chromosome discovery state replacing ModelStore class-level mutable globals."""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Optional

from isoquant_next.annotation.gene_info import StrandDetector, TranscriptModel
from isoquant_next.assignment.assigner import LongReadAssigner
from isoquant_next.assignment.profiles import CombinedProfileConstructor
from isoquant_next.assignment.types import ReadAssignment
from isoquant_next.config.enums import StrandnessReportingLevel
from isoquant_next.discovery.intron_graph import (
    IntronGraph,
    IntronPathProcessor,
    IntronPathStorage,
)


def reference_unique_known_isoform(assignment: ReadAssignment) -> Optional[str]:
    if assignment.assignment_type.is_unique() and assignment.isoform_matches:
        return assignment.isoform_matches[0].assigned_transcript
    return None


class ChromosomeDiscoveryState:
    """Explicit chromosome/worker-scoped dedup sets (replaces legacy ModelStore class globals)."""

    __slots__ = ("detected_known_isoforms", "extended_transcript_ids")

    def __init__(self) -> None:
        self.detected_known_isoforms: set[str] = set()
        self.extended_transcript_ids: set[str] = set()


class ModelStore:
    """Single source of truth for constructed models and read bindings within a locus."""

    def __init__(self, gene_info, chr_state: Optional[ChromosomeDiscoveryState] = None) -> None:
        self.gene_info = gene_info
        self.chr_state = chr_state if chr_state is not None else ChromosomeDiscoveryState()
        self.transcript_model_storage: List[TranscriptModel] = []
        self.transcript_read_ids: Dict[str, List[ReadAssignment]] = defaultdict(list)
        self.internal_counter: Dict[str, int] = defaultdict(int)
        self.read_assignment_counts: Dict[str, int] = defaultdict(int)
        self.construction_assignment: Dict[str, str] = {}
        self.nonfl_known_reads: Dict[str, List[ReadAssignment]] = {}

    @property
    def detected_known_isoforms(self) -> set[str]:
        return self.chr_state.detected_known_isoforms

    @property
    def extended_transcript_ids(self) -> set[str]:
        return self.chr_state.extended_transcript_ids

    def save_assigned_read(self, read_assignment: ReadAssignment, transcript_id: str) -> None:
        read_id = read_assignment.read_id
        self.transcript_read_ids[transcript_id].append(read_assignment)
        self.internal_counter[transcript_id] += 1
        self.read_assignment_counts[read_id] += 1
        self.construction_assignment[read_id] = transcript_id

    def move_read_to_model(self, read_assignment: ReadAssignment, source_id: str, target_id: str) -> None:
        self.transcript_read_ids[source_id].remove(read_assignment)
        self.internal_counter[source_id] -= 1
        self.read_assignment_counts[read_assignment.read_id] -= 1
        self.save_assigned_read(read_assignment, target_id)

    def delete_from_storage(self, transcript_id: str) -> None:
        for a in self.transcript_read_ids[transcript_id]:
            self.read_assignment_counts[a.read_id] -= 1
        del self.transcript_read_ids[transcript_id]
        del self.internal_counter[transcript_id]

    def transcript_from_reference(self, isoform_id: str) -> TranscriptModel:
        return TranscriptModel.from_reference_transcript(self.gene_info, isoform_id)

    def mapping_quality(self, transcript_id: str) -> float:
        reads = self.transcript_read_ids[transcript_id]
        return sum(a.mapping_quality for a in reads) / len(reads)


class ModelConstructionContext:
    """Read-only inputs + IntronGraph + explicit ModelStore shared across all 5 stages."""

    def __init__(
        self,
        gene_info,
        chr_record,
        args,
        id_distributor,
        string_pools,
        polya_predictions,
        tss_predictions,
        store: ModelStore,
        grouping_strategy_names=None,
        use_technical_replicas: bool = False,
    ) -> None:
        self.gene_info = gene_info
        self.chr_record = chr_record
        self.args = args
        self.id_distributor = id_distributor
        self.string_pools = string_pools
        self.polya_predictions = polya_predictions
        self.tss_predictions = tss_predictions
        self.store = store

        self.grouping_strategy_names = list(grouping_strategy_names) if grouping_strategy_names else []
        self.use_technical_replicas = use_technical_replicas
        self.file_name_group_idx = (
            self.grouping_strategy_names.index("file_name") if "file_name" in self.grouping_strategy_names else -1
        )

        self.create_alt_ts_nics = bool(getattr(args, "genedb", None))
        self.use_tss_model = bool(getattr(args, "fl_data", False))
        self.novel_apa = bool(getattr(args, "novel_apa", False))

        self.strand_detector = StrandDetector(self.chr_record)
        self.intron_genes = defaultdict(set)
        self.set_gene_properties()

        self.assigner = LongReadAssigner(self.gene_info, self.args, self.string_pools)
        self.profile_constructor = CombinedProfileConstructor(self.gene_info, self.args)

        self.intron_graph: Optional[IntronGraph] = None
        self.path_processor: Optional[IntronPathProcessor] = None
        self.path_storage: Optional[IntronPathStorage] = None
        self.known_isoforms_in_graph: Dict[tuple, str] = {}
        self.known_introns: set = set()
        self.known_isoforms_in_graph_ids: Dict[str, tuple] = {}

    def get_transcript_id(self):
        return self.id_distributor.increment()

    def set_gene_properties(self) -> None:
        intron_strands_dicts = defaultdict(lambda: defaultdict(int))
        self.intron_genes = defaultdict(set)
        for t_id, introns in self.gene_info.all_isoforms_introns.items():
            strand = self.gene_info.isoform_strands[t_id]
            gene_id = self.gene_info.gene_id_map[t_id]
            for intron in introns:
                intron_strands_dicts[intron][strand] += 1
                self.intron_genes[intron].add(gene_id)

        for intron, s_counts in intron_strands_dicts.items():
            if len(s_counts) == 1:
                self.strand_detector.set_strand(intron, next(iter(s_counts.keys())))
            else:
                self.strand_detector.set_strand(intron)

    def select_reference_gene(self, transcript_introns, transcript_range, transcript_strand: str) -> Optional[str]:
        if self.gene_info.empty():
            return None
        gene_counts = defaultdict(int)
        for intron in transcript_introns:
            if intron in self.intron_genes:
                for g_id in self.intron_genes[intron]:
                    gene_counts[g_id] += 1
        ordered = sorted(gene_counts.items(), key=lambda x: (x[1], x[0]), reverse=True)
        for gene_id, _ in ordered:
            if transcript_strand == "." or self.gene_info.gene_strands[gene_id] == transcript_strand:
                return gene_id
        return None

    def get_known_spliced_isoforms(self, gene_info) -> Dict[tuple, str]:
        known = {}
        for isoform_id, isoform_introns in gene_info.all_isoforms_introns.items():
            intron_path = self.path_processor.thread_introns(isoform_introns)
            if intron_path:
                known[tuple(intron_path)] = isoform_id
        return known

    def build_graph(self, read_assignment_storage) -> None:
        self.intron_graph = IntronGraph(
            self.args,
            self.gene_info,
            read_assignment_storage,
            polya_predictions=self.polya_predictions,
            tss_predictions=self.tss_predictions,
        )
        self.path_processor = IntronPathProcessor(self.args, self.intron_graph)
        self.path_storage = IntronPathStorage(self.args, self.path_processor)
        self.path_storage.fill(read_assignment_storage)
        self.known_isoforms_in_graph = self.get_known_spliced_isoforms(self.gene_info)
        self.known_introns = set(self.gene_info.intron_profiles.features)
        for intron_path, isoform_id in self.known_isoforms_in_graph.items():
            self.known_isoforms_in_graph_ids[isoform_id] = intron_path
