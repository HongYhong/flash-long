"""Enums and strategy alias tables matching IsoQuant 4.0.0 semantics."""

from __future__ import annotations

from enum import Enum

PACBIO_CCS_DATA = "pacbio_ccs"
NANOPORE_DATA = "nanopore"
ASSEMBLY = "assembly"

DATA_TYPE_ALIASES: dict[str, str] = {
    PACBIO_CCS_DATA: PACBIO_CCS_DATA,
    "pacbio": PACBIO_CCS_DATA,
    NANOPORE_DATA: NANOPORE_DATA,
    "ont": NANOPORE_DATA,
    ASSEMBLY: ASSEMBLY,
    "transcripts": ASSEMBLY,
}

DATATYPE_TO_ALIGNER: dict[str, str] = {
    ASSEMBLY: "minimap2",
    PACBIO_CCS_DATA: "minimap2",
    NANOPORE_DATA: "minimap2",
}


class IsoQuantMode(Enum):
    """Run modes for bulk, single-cell, and spatial protocols."""

    bulk = 1
    tenX_v3 = 2
    tenX_v2 = 22
    tenX_v3_split = 23
    tenX_v2_split = 24
    curio = 3
    stereoseq_nosplit = 4
    stereoseq = 5
    visium_hd = 61
    visium_5prime = 7
    custom_sc = 10

    def needs_barcode_calling(self) -> bool:
        return self != IsoQuantMode.bulk

    def needs_pcr_deduplication(self) -> bool:
        return self != IsoQuantMode.bulk

    def produces_new_fasta(self) -> bool:
        return self in (
            IsoQuantMode.stereoseq,
            IsoQuantMode.tenX_v3_split,
            IsoQuantMode.tenX_v2_split,
        )

    def needs_barcode_iterator(self) -> bool:
        return self in (IsoQuantMode.stereoseq, IsoQuantMode.stereoseq_nosplit)

    def enforces_single_thread(self) -> bool:
        return False


ISOQUANT_MODES: list[str] = [x.name for x in IsoQuantMode]


class AnalysisType(Enum):
    quantification = 1
    exon_quantification = 2
    fusion = 3
    transcript_discovery = 4


ANALYSIS_ALIASES: dict[str, AnalysisType] = {
    "quantification": AnalysisType.quantification,
    "quant": AnalysisType.quantification,
    "exon_quantification": AnalysisType.exon_quantification,
    "ex_quant": AnalysisType.exon_quantification,
    "fusion": AnalysisType.fusion,
    "transcript_discovery": AnalysisType.transcript_discovery,
    "td": AnalysisType.transcript_discovery,
}

ANALYSIS_CHOICES: list[str] = list(ANALYSIS_ALIASES.keys())


class PolyAUsageStrategies(Enum):
    auto = 1
    never = 2
    always = 3


def set_polya_requirement_strategy(
    flag: bool, polya_requirement_strategy: PolyAUsageStrategies
) -> bool:
    if polya_requirement_strategy == PolyAUsageStrategies.auto:
        return flag
    if polya_requirement_strategy == PolyAUsageStrategies.never:
        return False
    return True


class PolyATrimmed(Enum):
    none = 1
    all = 2
    stranded = 3


class StrandnessReportingLevel(Enum):
    only_canonical = 1
    only_stranded = 2
    all = 3
    auto = 10


class AmbiguityResolvingMethod(Enum):
    none = -10
    monoexon_only = 10
    monoexon_and_fsm = 20
    all = 30
    minimal_score = -0.5
    top_scored_factor = 1.5


class MultimapResolvingStrategy(Enum):
    ignore_multimapper = 0
    merge = 1
    take_best = 2


class CountingStrategy(Enum):
    unique_only = 1
    with_ambiguous = 2
    fsm_only = 3
    unique_splicing_consistent = 10
    unique_inconsistent = 11
    all = 20

    def no_inconsistent(self) -> bool:
        return self in (CountingStrategy.unique_only, CountingStrategy.with_ambiguous)

    def ambiguous(self) -> bool:
        return self in (CountingStrategy.all, CountingStrategy.with_ambiguous)

    def inconsistent_minor(self) -> bool:
        return self in (
            CountingStrategy.unique_splicing_consistent,
            CountingStrategy.unique_inconsistent,
            CountingStrategy.all,
        )

    def inconsistent(self) -> bool:
        return self in (CountingStrategy.unique_inconsistent, CountingStrategy.all)


COUNTING_STRATEGIES: list[str] = [
    CountingStrategy.unique_only.name,
    CountingStrategy.with_ambiguous.name,
    CountingStrategy.fsm_only.name,
    CountingStrategy.unique_splicing_consistent.name,
    CountingStrategy.unique_inconsistent.name,
    CountingStrategy.all.name,
]


class GroupedOutputFormat(Enum):
    none = 0
    matrix = 1
    mtx = 2
    default = 3


class NormalizationMethod(Enum):
    none = 0
    simple = 1
    usable_reads = 2
