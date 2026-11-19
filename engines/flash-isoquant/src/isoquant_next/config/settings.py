"""Immutable strongly-typed configuration contracts and JSON manifest persistence."""

from __future__ import annotations

import hashlib
import json
import os
from dataclasses import asdict, dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Mapping, Optional, Tuple

from isoquant_next import __version__
from isoquant_next.config.enums import (
    ASSEMBLY,
    AmbiguityResolvingMethod,
    AnalysisType,
    CountingStrategy,
    DATA_TYPE_ALIASES,
    IsoQuantMode,
    MultimapResolvingStrategy,
    NANOPORE_DATA,
    PACBIO_CCS_DATA,
    PolyATrimmed,
    PolyAUsageStrategies,
    StrandnessReportingLevel,
)
from isoquant_next.core.error_codes import IsoQuantError, IsoQuantExitCode

MANIFEST_SCHEMA_VERSION = "1.0"


@dataclass(frozen=True, slots=True)
class MatchingConfig:
    strategy_name: str = "precise"
    delta: int = 4
    minor_exon_extension: int = 50
    major_exon_extension: int = 300
    max_intron_shift: int = 30
    max_missed_exon_len: int = 50
    max_fake_terminal_exon_len: int = 20
    max_suspicious_intron_abs_len: int = 0
    max_suspicious_intron_rel_len: float = 0.0
    min_abs_exon_overlap: int = 10
    min_rel_exon_overlap: float = 0.2
    micro_intron_length: int = 50
    max_intron_abs_diff: int = 30
    max_intron_rel_diff: float = 0.2
    apa_delta: int = 50
    minimal_exon_overlap: int = 5
    minimal_intron_absence_overlap: int = 20
    polya_window: int = 16
    polya_fraction: float = 0.75
    resolve_ambiguous: AmbiguityResolvingMethod = AmbiguityResolvingMethod.monoexon_and_fsm
    correct_minor_errors: bool = True
    multimap_strategy: MultimapResolvingStrategy = MultimapResolvingStrategy.take_best
    simple_models_mapq_cutoff: int = 30
    polya_percentage_threshold: float = 0.7
    low_polya_percentage_threshold: float = 0.1
    indel_near_splice_site_dist: int = 10
    upstream_region_len: int = 20


MATCHING_STRATEGY_PRESETS: dict[str, MatchingConfig] = {
    "exact": MatchingConfig(
        strategy_name="exact",
        delta=0,
        max_intron_shift=0,
        max_intron_abs_diff=0,
        max_missed_exon_len=0,
        max_fake_terminal_exon_len=0,
        max_suspicious_intron_abs_len=0,
        max_suspicious_intron_rel_len=0.0,
        resolve_ambiguous=AmbiguityResolvingMethod.monoexon_only,
        correct_minor_errors=False,
    ),
    "precise": MatchingConfig(
        strategy_name="precise",
        delta=4,
        max_intron_shift=30,
        max_intron_abs_diff=30,
        max_missed_exon_len=50,
        max_fake_terminal_exon_len=20,
        max_suspicious_intron_abs_len=0,
        max_suspicious_intron_rel_len=0.0,
        resolve_ambiguous=AmbiguityResolvingMethod.monoexon_and_fsm,
        correct_minor_errors=True,
    ),
    "default": MatchingConfig(
        strategy_name="default",
        delta=6,
        max_intron_shift=60,
        max_intron_abs_diff=30,
        max_missed_exon_len=100,
        max_fake_terminal_exon_len=40,
        max_suspicious_intron_abs_len=60,
        max_suspicious_intron_rel_len=1.0,
        resolve_ambiguous=AmbiguityResolvingMethod.monoexon_and_fsm,
        correct_minor_errors=True,
    ),
    "loose": MatchingConfig(
        strategy_name="loose",
        delta=12,
        max_intron_shift=60,
        max_intron_abs_diff=30,
        max_missed_exon_len=100,
        max_fake_terminal_exon_len=40,
        max_suspicious_intron_abs_len=60,
        max_suspicious_intron_rel_len=1.0,
        resolve_ambiguous=AmbiguityResolvingMethod.all,
        correct_minor_errors=True,
    ),
}


@dataclass(frozen=True, slots=True)
class SpliceCorrectionConfig:
    strategy_name: str = "default_pacbio"
    correct_fuzzy_junctions: bool = True
    correct_intron_shifts: bool = False
    correct_skipped_exons: bool = True
    correct_terminal_exons: bool = False
    correct_fake_terminal_exons: bool = False
    correct_microintron_retention: bool = True


SPLICE_CORRECTION_PRESETS: dict[str, SpliceCorrectionConfig] = {
    "none": SpliceCorrectionConfig("none", False, False, False, False, False, False),
    "default_pacbio": SpliceCorrectionConfig("default_pacbio", True, False, True, False, False, True),
    "conservative_ont": SpliceCorrectionConfig("conservative_ont", True, False, True, False, False, False),
    "default_ont": SpliceCorrectionConfig("default_ont", True, False, True, False, True, True),
    "all": SpliceCorrectionConfig("all", True, True, True, True, True, True),
    "assembly": SpliceCorrectionConfig("assembly", False, False, True, False, False, False),
}


@dataclass(frozen=True, slots=True)
class ModelConstructionConfig:
    strategy_name: str = "default_pacbio"
    min_novel_intron_count: int = 1
    graph_clustering_ratio: float = 0.5
    graph_clustering_distance: int = 10
    min_novel_isolated_intron_abs: int = 2
    min_novel_isolated_intron_rel: float = 0.02
    terminal_position_abs: int = 1
    terminal_position_rel: float = 0.05
    terminal_internal_position_rel: float = 0.05
    min_known_count: int = 1
    min_nonfl_count: int = 2
    min_novel_count: int = 2
    min_novel_count_rel: float = 0.02
    min_mono_count_rel: float = 0.005
    singleton_adjacent_cov: int = 100
    fl_only: bool = False
    report_novel_unspliced: bool = True
    min_mono_exon_coverage: float = 0.75
    require_monointronic_polya: bool = False
    require_monoexonic_polya: bool = True
    polya_requirement_strategy: PolyAUsageStrategies = PolyAUsageStrategies.auto
    polya_trimmed: PolyATrimmed = PolyATrimmed.none
    report_canonical: StrandnessReportingLevel = StrandnessReportingLevel.only_canonical


MODEL_CONSTRUCTION_PRESETS: dict[str, ModelConstructionConfig] = {
    "reliable": ModelConstructionConfig(
        "reliable", 2, 0.5, 20, 5, 0.05, 1, 0.1, 0.1, 2, 4, 8, 0.05, 0.05, 50,
        True, False, 0.75, True, True, PolyAUsageStrategies.auto, PolyATrimmed.none,
        StrandnessReportingLevel.only_canonical,
    ),
    "default_pacbio": ModelConstructionConfig(
        "default_pacbio", 1, 0.5, 10, 2, 0.02, 1, 0.05, 0.05, 1, 2, 2, 0.02, 0.005, 100,
        False, True, 0.75, False, True, PolyAUsageStrategies.auto, PolyATrimmed.none,
        StrandnessReportingLevel.only_canonical,
    ),
    "sensitive_pacbio": ModelConstructionConfig(
        "sensitive_pacbio", 1, 0.5, 5, 2, 0.005, 1, 0.01, 0.02, 1, 2, 2, 0.005, 0.001, 100,
        False, True, 0.75, False, False, PolyAUsageStrategies.auto, PolyATrimmed.none,
        StrandnessReportingLevel.only_stranded,
    ),
    "pacbio_all": ModelConstructionConfig(
        "pacbio_all", 1, 0.5, 5, 2, 0.005, 1, 0.01, 0.02, 1, 2, 1, 0.002, 0.001, 100,
        False, True, 0.75, False, False, PolyAUsageStrategies.auto, PolyATrimmed.none,
        StrandnessReportingLevel.only_stranded,
    ),
    "default_ont": ModelConstructionConfig(
        "default_ont", 1, 0.5, 20, 3, 0.02, 1, 0.05, 0.05, 1, 3, 3, 0.02, 0.02, 10,
        False, False, 0.75, True, True, PolyAUsageStrategies.auto, PolyATrimmed.none,
        StrandnessReportingLevel.only_canonical,
    ),
    "sensitive_ont": ModelConstructionConfig(
        "sensitive_ont", 1, 0.5, 20, 3, 0.005, 1, 0.01, 0.02, 1, 2, 3, 0.005, 0.005, 10,
        False, True, 0.75, False, False, PolyAUsageStrategies.auto, PolyATrimmed.none,
        StrandnessReportingLevel.only_stranded,
    ),
    "fl_pacbio": ModelConstructionConfig(
        "fl_pacbio", 1, 0.5, 10, 2, 0.02, 1, 0.05, 0.01, 1, 2, 3, 0.02, 0.005, 100,
        True, True, 0.75, False, False, PolyAUsageStrategies.auto, PolyATrimmed.none,
        StrandnessReportingLevel.only_canonical,
    ),
    "all": ModelConstructionConfig(
        "all", 0, 0.3, 5, 1, 0.002, 1, 0.01, 0.01, 1, 1, 1, 0.002, 0.001, 500,
        False, True, 0.75, False, False, PolyAUsageStrategies.auto, PolyATrimmed.none,
        StrandnessReportingLevel.all,
    ),
    "assembly": ModelConstructionConfig(
        "assembly", 0, 0.3, 5, 1, 0.05, 1, 0.01, 0.02, 1, 1, 1, 0.05, 0.01, 50,
        False, True, 0.75, False, False, PolyAUsageStrategies.auto, PolyATrimmed.none,
        StrandnessReportingLevel.only_stranded,
    ),
}


@dataclass(frozen=True, slots=True)
class QuantificationConfig:
    transcript_strategy: CountingStrategy = CountingStrategy.unique_only
    gene_strategy: CountingStrategy = CountingStrategy.unique_only
    count_exons: bool = False
    count_intron_retentions: bool = False
    counts_format: str = "default"
    normalization: str = "none"


@dataclass(frozen=True, slots=True)
class BarcodeConfig:
    mode: IsoQuantMode = IsoQuantMode.bulk
    barcoded_reads: bool = False
    barcode2spot: Optional[str] = None
    barcode2barcode: Optional[str] = None
    umi_length: int = 0
    read_group: Tuple[str, ...] = ()
    large_output: Tuple[str, ...] = ()
    bam_tags: Tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class RuntimeConfig:
    threads: int = 16
    output_dir: str = "./isoquant_output"
    prefix: str = "OUT"
    clean_start: bool = False
    resume: bool = False
    low_memory: bool = True
    mp_context: str = "forkserver"
    pin_math_threads: bool = True
    keep_tmp: bool = False
    use_mmap_index: bool = True


@dataclass(frozen=True, slots=True)
class PipelineConfig:
    """Top-level immutable pipeline configuration replacing dynamic Namespace mutations."""

    data_type: str
    fl_data: bool = False
    reference: Optional[str] = None
    genedb: Optional[str] = None
    complete_genedb: bool = False
    analyses: Tuple[AnalysisType, ...] = (AnalysisType.quantification, AnalysisType.transcript_discovery)
    matching: MatchingConfig = field(default_factory=MatchingConfig)
    splice_correction: SpliceCorrectionConfig = field(default_factory=SpliceCorrectionConfig)
    model_construction: ModelConstructionConfig = field(default_factory=ModelConstructionConfig)
    quantification: QuantificationConfig = field(default_factory=QuantificationConfig)
    barcodes: BarcodeConfig = field(default_factory=BarcodeConfig)
    runtime: RuntimeConfig = field(default_factory=RuntimeConfig)

    @classmethod
    def resolve(
        cls,
        *,
        data_type: str,
        mode: IsoQuantMode = IsoQuantMode.bulk,
        fl_data: bool = False,
        reference: Optional[str] = None,
        genedb: Optional[str] = None,
        complete_genedb: bool = False,
        matching_strategy: Optional[str] = None,
        delta: Optional[int] = None,
        splice_correction_strategy: Optional[str] = None,
        model_construction_strategy: Optional[str] = None,
        graph_clustering_distance: Optional[int] = None,
        report_novel_unspliced: Optional[bool] = None,
        polya_requirement: str = "auto",
        polya_trimmed: str = "none",
        report_canonical: str = "auto",
        read_group: Optional[list[str]] = None,
        has_replicas: bool = False,
        barcode2spot: Optional[str] = None,
        barcode2barcode: Optional[str] = None,
        bam_tags: Optional[str] = None,
        large_output: Optional[list[str]] = None,
        analyses: Optional[Tuple[AnalysisType, ...]] = None,
        threads: int = 16,
        output_dir: str = "./isoquant_output",
        prefix: str = "OUT",
        resume: bool = False,
        clean_start: bool = False,
    ) -> "PipelineConfig":
        canonical_dt = DATA_TYPE_ALIASES.get(data_type)
        if canonical_dt is None:
            raise IsoQuantError(
                IsoQuantExitCode.INVALID_PARAMETER,
                f"Unsupported data_type {data_type!r}; choose from {sorted(DATA_TYPE_ALIASES)}",
            )

        if analyses is None:
            if mode.needs_barcode_calling():
                analyses = (AnalysisType.quantification,)
            elif genedb:
                analyses = (AnalysisType.quantification, AnalysisType.transcript_discovery)
            else:
                analyses = (AnalysisType.transcript_discovery,)

        default_matching = {
            ASSEMBLY: "precise",
            PACBIO_CCS_DATA: "precise",
            NANOPORE_DATA: "default",
        }[canonical_dt]
        ms_name = matching_strategy or default_matching
        if ms_name not in MATCHING_STRATEGY_PRESETS:
            raise IsoQuantError(IsoQuantExitCode.INVALID_PARAMETER, f"Unknown matching_strategy: {ms_name}")
        m_cfg = MATCHING_STRATEGY_PRESETS[ms_name]
        if delta is not None:
            if delta < 0:
                raise IsoQuantError(IsoQuantExitCode.INVALID_PARAMETER, "--delta can not be negative")
            m_cfg = MatchingConfig(**{**asdict(m_cfg), "delta": delta})

        default_sc = {
            ASSEMBLY: "assembly",
            PACBIO_CCS_DATA: "default_pacbio",
            NANOPORE_DATA: "default_ont",
        }[canonical_dt]
        sc_name = splice_correction_strategy or default_sc
        if sc_name not in SPLICE_CORRECTION_PRESETS:
            raise IsoQuantError(IsoQuantExitCode.INVALID_PARAMETER, f"Unknown splice_correction_strategy: {sc_name}")
        sc_cfg = SPLICE_CORRECTION_PRESETS[sc_name]

        default_mc = {
            ASSEMBLY: "assembly",
            PACBIO_CCS_DATA: "fl_pacbio" if fl_data else "default_pacbio",
            NANOPORE_DATA: "default_ont",
        }[canonical_dt]
        mc_name = model_construction_strategy or default_mc
        if mc_name not in MODEL_CONSTRUCTION_PRESETS:
            raise IsoQuantError(IsoQuantExitCode.INVALID_PARAMETER, f"Unknown model_construction_strategy: {mc_name}")
        mc_preset = MODEL_CONSTRUCTION_PRESETS[mc_name]
        if graph_clustering_distance is not None and graph_clustering_distance < 0:
            raise IsoQuantError(
                IsoQuantExitCode.INVALID_PARAMETER,
                "--graph_clustering_distance can not be negative",
            )
        canon_enum = StrandnessReportingLevel[report_canonical]
        if canon_enum == StrandnessReportingLevel.auto:
            canon_enum = mc_preset.report_canonical
        mc_cfg = ModelConstructionConfig(
            **{
                **asdict(mc_preset),
                "graph_clustering_distance": (
                    mc_preset.graph_clustering_distance
                    if graph_clustering_distance is None
                    else graph_clustering_distance
                ),
                "report_novel_unspliced": (
                    mc_preset.report_novel_unspliced
                    if report_novel_unspliced is None
                    else report_novel_unspliced
                ),
                "polya_requirement_strategy": PolyAUsageStrategies[polya_requirement],
                "polya_trimmed": PolyATrimmed[polya_trimmed],
                "report_canonical": canon_enum,
            }
        )

        # Resolve read groups & large output matching isoquant.py set_data_dependent_options
        rg = list(read_group) if read_group is not None else None
        lo = list(large_output) if large_output is not None else None
        if rg is not None and "none" in rg:
            rg_final: Tuple[str, ...] = ()
        else:
            rg_work = [] if rg is None else rg
            if has_replicas and "file_name" not in rg_work and "no_auto" not in rg_work:
                rg_work.append("file_name")
            if barcode2spot and "barcode_spot" not in rg_work and "no_auto" not in rg_work:
                rg_work.append("barcode_spot")
            if barcode2barcode:
                if "barcode_barcode" not in rg_work and "no_auto" not in rg_work:
                    rg_work.append("barcode_barcode")
                if lo is None:
                    lo = ["allinfo"]
                elif "none" not in lo and "allinfo" not in lo:
                    lo.append("allinfo")
            if mode.needs_barcode_calling():
                bc_keys = {"barcode", "barcode_spot", "barcode_barcode", "no_auto"}
                if not any(item in bc_keys for item in rg_work):
                    rg_work.append("barcode")
            rg_final = tuple(rg_work)
        lo_final = () if lo is None or "none" in lo else tuple(lo)

        tags_tuple = tuple(t for t in bam_tags.split(",") if t) if bam_tags else ()
        bc_cfg = BarcodeConfig(
            mode=mode,
            barcoded_reads=mode.needs_barcode_calling(),
            barcode2spot=barcode2spot,
            barcode2barcode=barcode2barcode,
            read_group=rg_final,
            large_output=lo_final,
            bam_tags=tags_tuple,
        )
        rt_cfg = RuntimeConfig(
            threads=max(1, int(threads)),
            output_dir=output_dir,
            prefix=prefix,
            clean_start=clean_start,
            resume=resume,
        )
        return cls(
            data_type=canonical_dt,
            fl_data=fl_data,
            reference=reference,
            genedb=genedb,
            complete_genedb=complete_genedb,
            analyses=tuple(analyses),
            matching=m_cfg,
            splice_correction=sc_cfg,
            model_construction=mc_cfg,
            barcodes=bc_cfg,
            runtime=rt_cfg,
        )

    def command_fingerprint(self) -> str:
        cfg_dict = dict(self.to_manifest_dict()["config"])
        rt = dict(cfg_dict.get("runtime", {}))
        for volatile_key in ("threads", "clean_start", "resume", "keep_tmp"):
            rt.pop(volatile_key, None)
        cfg_dict["runtime"] = rt
        payload = json.dumps(cfg_dict, sort_keys=True, separators=(",", ":"))
        return hashlib.blake2b(payload.encode("utf-8"), digest_size=16).hexdigest()

    def to_manifest_dict(self) -> dict[str, Any]:
        def _encode(val: Any) -> Any:
            if isinstance(val, Enum):
                return val.name
            if isinstance(val, tuple):
                return [_encode(v) for v in val]
            if isinstance(val, dict):
                return {k: _encode(v) for k, v in val.items()}
            if isinstance(val, list):
                return [_encode(v) for v in val]
            return val

        raw = _encode(asdict(self))
        return {
            "schema_version": MANIFEST_SCHEMA_VERSION,
            "isoquant_next_version": __version__,
            "fingerprint": hashlib.blake2b(
                json.dumps(raw, sort_keys=True, separators=(",", ":")).encode("utf-8"),
                digest_size=16,
            ).hexdigest(),
            "config": raw,
        }

    def save_manifest(self, path: str | Path) -> Path:
        out = Path(path)
        out.parent.mkdir(parents=True, exist_ok=True)
        tmp = out.with_suffix(out.suffix + f".tmp.{os.getpid()}")
        tmp.write_text(json.dumps(self.to_manifest_dict(), indent=2) + "\n", encoding="utf-8")
        os.replace(tmp, out)
        return out

    @classmethod
    def load_manifest(cls, path: str | Path) -> "PipelineConfig":
        in_path = Path(path)
        if not in_path.exists():
            raise IsoQuantError(
                IsoQuantExitCode.RESUME_CONFIG_NOT_FOUND,
                f"Manifest not found: {in_path}",
            )
        doc = json.loads(in_path.read_text(encoding="utf-8"))
        if str(doc.get("schema_version")) != MANIFEST_SCHEMA_VERSION:
            raise IsoQuantError(
                IsoQuantExitCode.FILE_FORMAT_ERROR,
                f"Unsupported manifest schema_version: {doc.get('schema_version')}",
            )
        c: Mapping[str, Any] = doc["config"]

        m = dict(c["matching"])
        m["resolve_ambiguous"] = AmbiguityResolvingMethod[m["resolve_ambiguous"]]
        m["multimap_strategy"] = MultimapResolvingStrategy[m["multimap_strategy"]]

        mc = dict(c["model_construction"])
        mc["polya_requirement_strategy"] = PolyAUsageStrategies[mc["polya_requirement_strategy"]]
        mc["polya_trimmed"] = PolyATrimmed[mc["polya_trimmed"]]
        mc["report_canonical"] = StrandnessReportingLevel[mc["report_canonical"]]

        q = dict(c["quantification"])
        q["transcript_strategy"] = CountingStrategy[q["transcript_strategy"]]
        q["gene_strategy"] = CountingStrategy[q["gene_strategy"]]

        b = dict(c["barcodes"])
        b["mode"] = IsoQuantMode[b["mode"]]
        b["read_group"] = tuple(b["read_group"])
        b["large_output"] = tuple(b["large_output"])
        b["bam_tags"] = tuple(b["bam_tags"])

        return cls(
            data_type=c["data_type"],
            fl_data=bool(c["fl_data"]),
            reference=c["reference"],
            genedb=c["genedb"],
            complete_genedb=bool(c["complete_genedb"]),
            analyses=tuple(AnalysisType[a] for a in c["analyses"]),
            matching=MatchingConfig(**m),
            splice_correction=SpliceCorrectionConfig(**c["splice_correction"]),
            model_construction=ModelConstructionConfig(**mc),
            quantification=QuantificationConfig(**q),
            barcodes=BarcodeConfig(**b),
            runtime=RuntimeConfig(**c["runtime"]),
        )
