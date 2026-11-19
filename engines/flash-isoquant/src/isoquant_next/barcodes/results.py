"""Typed detection result containers and statistics aggregators implementing BarcodeResult."""

from __future__ import annotations

from collections import defaultdict
from typing import ClassVar, Dict, Iterable, List, Optional, Protocol, Union, runtime_checkable

from isoquant_next.barcodes.molecule_structure import ElementType, MoleculeStructure


def increase_if_valid(val: Optional[int], delta: int) -> Optional[int]:
    if val and val != -1:
        return val + delta
    return val


@runtime_checkable
class BarcodeResult(Protocol):
    NOSEQ: ClassVar[str]
    read_id: str
    strand: str

    def get_barcode(self) -> str: ...
    def get_umi(self) -> str: ...
    def is_valid(self) -> bool: ...
    def has_barcode(self) -> bool: ...
    def has_umi(self) -> bool: ...


class BarcodeDetectionResult:
    NOSEQ = "*"

    def __init__(
        self,
        read_id: str,
        barcode: str = NOSEQ,
        UMI: str = NOSEQ,
        BC_score: int = -1,
        UMI_good: bool = False,
        strand: str = ".",
    ) -> None:
        self.read_id: str = read_id
        self._barcode: str = barcode if barcode else self.NOSEQ
        self._umi: str = UMI if UMI else self.NOSEQ
        self.BC_score: int = BC_score
        self.UMI_good: bool = UMI_good
        self.strand: str = strand

    def get_barcode(self) -> str:
        return self._barcode

    def get_umi(self) -> str:
        return self._umi

    @property
    def barcode(self) -> str:
        return self._barcode

    @barcode.setter
    def barcode(self, value: str) -> None:
        self._barcode = value if value else self.NOSEQ

    @property
    def UMI(self) -> str:
        return self._umi

    @UMI.setter
    def UMI(self, value: str) -> None:
        self._umi = value if value else self.NOSEQ

    def is_valid(self) -> bool:
        return self._barcode != self.NOSEQ

    def has_barcode(self) -> bool:
        return self._barcode != self.NOSEQ

    def has_umi(self) -> bool:
        return self._umi != self.NOSEQ

    def update_coordinates(self, delta: int) -> None:
        pass

    def more_informative_than(self, that: "BarcodeDetectionResult") -> bool:
        raise NotImplementedError()

    def get_additional_attributes(self) -> List[str]:
        return []

    def set_strand(self, strand: str) -> None:
        self.strand = strand

    def __str__(self) -> str:
        return f"{self.read_id}\t{self._barcode}\t{self._umi}\t{self.BC_score}\t{self.UMI_good}\t{self.strand}"

    @staticmethod
    def header() -> str:
        return "read_id\tbarcode\tUMI\tBC_score\tvalid_UMI\tstrand"


class LinkerBarcodeDetectionResult(BarcodeDetectionResult):
    def __init__(
        self,
        read_id: str,
        barcode: str = BarcodeDetectionResult.NOSEQ,
        UMI: str = BarcodeDetectionResult.NOSEQ,
        BC_score: int = -1,
        UMI_good: bool = False,
        strand: str = ".",
        polyT: int = -1,
        primer: int = -1,
        linker_start: int = -1,
        linker_end: int = -1,
    ) -> None:
        super().__init__(read_id, barcode, UMI, BC_score, UMI_good, strand)
        self.primer: int = primer
        self.linker_start: int = linker_start
        self.linker_end: int = linker_end
        self.polyT: int = polyT

    def update_coordinates(self, delta: int) -> None:
        self.primer = increase_if_valid(self.primer, delta)
        self.linker_start = increase_if_valid(self.linker_start, delta)
        self.linker_end = increase_if_valid(self.linker_end, delta)
        self.polyT = increase_if_valid(self.polyT, delta)

    def more_informative_than(self, that: "LinkerBarcodeDetectionResult") -> bool:
        if self.BC_score != that.BC_score:
            return self.BC_score > that.BC_score
        if self.linker_start != that.linker_start:
            return self.linker_start > that.linker_start
        if self.primer != that.primer:
            return self.primer > that.primer
        return self.polyT > that.polyT

    def get_additional_attributes(self) -> List[str]:
        attr = []
        if self.polyT != -1:
            attr.append("PolyT detected")
        if self.primer != -1:
            attr.append("Primer detected")
        if self.linker_start != -1:
            attr.append("Linker detected")
        return attr

    def __str__(self) -> str:
        return f"{super().__str__()}\t{self.polyT}\t{self.primer}\t{self.linker_start}\t{self.linker_end}"

    @staticmethod
    def header() -> str:
        return BarcodeDetectionResult.header() + "\tpolyT_start\tprimer_end\tlinker_start\tlinker_end"


class TSOBarcodeDetectionResult(LinkerBarcodeDetectionResult):
    def __init__(
        self,
        read_id: str,
        barcode: str = BarcodeDetectionResult.NOSEQ,
        UMI: str = BarcodeDetectionResult.NOSEQ,
        BC_score: int = -1,
        UMI_good: bool = False,
        strand: str = ".",
        polyT: int = -1,
        primer: int = -1,
        linker_start: int = -1,
        linker_end: int = -1,
        tso: int = -1,
    ) -> None:
        super().__init__(read_id, barcode, UMI, BC_score, UMI_good, strand, polyT, primer, linker_start, linker_end)
        self.tso5: int = tso

    def update_coordinates(self, delta: int) -> None:
        self.tso5 = increase_if_valid(self.tso5, delta)
        super().update_coordinates(delta)

    def __str__(self) -> str:
        return f"{super().__str__()}\t{self.tso5}"

    def get_additional_attributes(self) -> List[str]:
        attr = super().get_additional_attributes()
        if self.tso5 != -1:
            attr.append("TSO detected")
        return attr

    def get_fasta_segment_start(self) -> int:
        return max(0, self.primer - 25, self.polyT - 75)

    def get_fasta_segment_end(self, seq_len: int) -> int:
        if self.tso5 == -1:
            return seq_len
        return min(seq_len, self.tso5 + 25)

    def get_tso_position(self) -> int:
        return self.tso5

    @staticmethod
    def header() -> str:
        return LinkerBarcodeDetectionResult.header() + "\tTSO5"


class TenXBarcodeDetectionResult(BarcodeDetectionResult):
    def __init__(
        self,
        read_id: str,
        barcode: str = BarcodeDetectionResult.NOSEQ,
        UMI: str = BarcodeDetectionResult.NOSEQ,
        BC_score: int = -1,
        UMI_good: bool = False,
        strand: str = ".",
        polyT: int = -1,
        r1: int = -1,
    ) -> None:
        super().__init__(read_id, barcode, UMI, BC_score, UMI_good, strand)
        self.r1: int = r1
        self.polyT: int = polyT

    def update_coordinates(self, delta: int) -> None:
        self.r1 = increase_if_valid(self.r1, delta)
        self.polyT = increase_if_valid(self.polyT, delta)

    def more_informative_than(self, that: "TenXBarcodeDetectionResult") -> bool:
        if self.polyT != that.polyT:
            return self.polyT > that.polyT
        if self.r1 != that.r1:
            return self.r1 > that.r1
        return self.BC_score > that.BC_score

    def get_additional_attributes(self) -> List[str]:
        attr = []
        if self.polyT != -1:
            attr.append("PolyT detected")
        if self.r1 != -1:
            attr.append("R1 detected")
        return attr

    def __str__(self) -> str:
        return f"{super().__str__()}\t{self.polyT}\t{self.r1}"

    @staticmethod
    def header() -> str:
        return BarcodeDetectionResult.header() + "\tpolyT_start\tR1_end"


class TenXSplitBarcodeDetectionResult(TenXBarcodeDetectionResult):
    def __init__(
        self,
        read_id: str,
        barcode: str = BarcodeDetectionResult.NOSEQ,
        UMI: str = BarcodeDetectionResult.NOSEQ,
        BC_score: int = -1,
        UMI_good: bool = False,
        strand: str = ".",
        polyT: int = -1,
        r1: int = -1,
        tso: int = -1,
    ) -> None:
        super().__init__(read_id, barcode, UMI, BC_score, UMI_good, strand, polyT, r1)
        self.tso: int = tso

    def update_coordinates(self, delta: int) -> None:
        super().update_coordinates(delta)
        if self.tso != -1:
            self.tso += delta

    def get_fasta_segment_start(self) -> int:
        return max(0, self.r1 - 10) if self.r1 != -1 else 0

    def get_fasta_segment_end(self, seq_len: int) -> int:
        if self.tso == -1:
            return seq_len
        return min(seq_len, self.tso + 35)

    def get_tso_position(self) -> int:
        return self.tso

    def __str__(self) -> str:
        return f"{super().__str__()}\t{self.tso}"

    @staticmethod
    def header() -> str:
        return TenXBarcodeDetectionResult.header() + "\ttso_start"


class SplittingBarcodeDetectionResult:
    NOSEQ = BarcodeDetectionResult.NOSEQ

    def __init__(self, read_id: str) -> None:
        self.read_id: str = read_id
        self.strand: str = "."
        self.detected_patterns: List[Union[TSOBarcodeDetectionResult, TenXSplitBarcodeDetectionResult]] = []

    def append(self, barcode_detection_result: Union[TSOBarcodeDetectionResult, TenXSplitBarcodeDetectionResult]) -> None:
        self.detected_patterns.append(barcode_detection_result)

    def empty(self) -> bool:
        return not self.detected_patterns

    def filter(self) -> None:
        if not self.detected_patterns:
            return
        barcoded_results = [r for r in self.detected_patterns if r.is_valid()]
        self.detected_patterns = barcoded_results if barcoded_results else [self.detected_patterns[0]]

    def get_barcode(self) -> str:
        for r in self.detected_patterns:
            if r.is_valid():
                return r.get_barcode()
        return self.NOSEQ

    def get_umi(self) -> str:
        for r in self.detected_patterns:
            if r.is_valid():
                return r.get_umi()
        return self.NOSEQ

    def is_valid(self) -> bool:
        return any(r.is_valid() for r in self.detected_patterns)

    def has_barcode(self) -> bool:
        return any(r.has_barcode() for r in self.detected_patterns)

    def has_umi(self) -> bool:
        return any(r.has_umi() for r in self.detected_patterns)

    def set_strand(self, strand: str) -> None:
        self.strand = strand
        for r in self.detected_patterns:
            r.set_strand(strand)

    def update_coordinates(self, delta: int) -> None:
        for r in self.detected_patterns:
            r.update_coordinates(delta)

    def more_informative_than(self, other: "SplittingBarcodeDetectionResult") -> bool:
        self_valid = sum(1 for r in self.detected_patterns if r.is_valid())
        other_valid = sum(1 for r in other.detected_patterns if r.is_valid())
        return self_valid > other_valid

    def get_additional_attributes(self) -> List[str]:
        attrs = set()
        for r in self.detected_patterns:
            attrs.update(r.get_additional_attributes())
        return list(attrs)

    def __str__(self) -> str:
        return "\n".join(str(r) for r in self.detected_patterns)

    @staticmethod
    def header() -> str:
        return TSOBarcodeDetectionResult.header()


class DetectedElement:
    def __init__(self, start: int = -1, end: int = -1, score: int = -1, seq: Optional[str] = None) -> None:
        self.start = start
        self.end = end
        self.score = score
        self.seq = seq


class ExtractionResult:
    NOSEQ = "*"

    def __init__(
        self,
        molecule_structure: MoleculeStructure,
        read_id: str,
        strand: str,
        detected_results: Dict[str, DetectedElement],
    ) -> None:
        self.molecule_structure = molecule_structure
        self.read_id = read_id
        self.strand = strand
        self.detected_results = detected_results

    @property
    def _barcode_elements(self) -> List[str]:
        return self.molecule_structure.barcode_elements

    @property
    def _umi_elements(self) -> List[str]:
        return self.molecule_structure.umi_elements

    def get_barcode(self) -> str:
        parts = []
        for el_name in self._barcode_elements:
            det = self.detected_results.get(el_name)
            if det and det.seq and det.seq != self.NOSEQ:
                parts.append(det.seq)
            else:
                return self.NOSEQ
        return "|".join(parts) if parts else self.NOSEQ

    def get_umi(self) -> str:
        parts = []
        for el_name in self._umi_elements:
            det = self.detected_results.get(el_name)
            if det and det.seq and det.seq != self.NOSEQ:
                parts.append(det.seq)
            else:
                return self.NOSEQ
        return "|".join(parts) if parts else self.NOSEQ

    def is_valid(self) -> bool:
        return self.get_barcode() != self.NOSEQ

    def has_barcode(self) -> bool:
        return self.is_valid()

    def has_umi(self) -> bool:
        return self.get_umi() != self.NOSEQ

    def get_barcode_score(self) -> int:
        max_score = -1
        for el_name in self._barcode_elements:
            det = self.detected_results.get(el_name)
            if det is not None and det.score > max_score:
                max_score = det.score
        return max_score

    def set_strand(self, strand: str) -> None:
        self.strand = strand

    def update_coordinates(self, delta: int) -> None:
        for det in self.detected_results.values():
            if det.start != -1:
                det.start += delta
            if det.end != -1:
                det.end += delta

    def more_informative_than(self, other: "ExtractionResult") -> bool:
        self_count = sum(1 for d in self.detected_results.values() if d.start != -1)
        other_count = sum(1 for d in other.detected_results.values() if d.start != -1)
        if self_count != other_count:
            return self_count > other_count
        self_score = sum(d.score for d in self.detected_results.values() if d.score > 0)
        other_score = sum(d.score for d in other.detected_results.values() if d.score > 0)
        return self_score > other_score

    def get_additional_attributes(self) -> List[str]:
        barcode_umi_set = set(self._barcode_elements + self._umi_elements)
        return [
            f"{el_name} detected"
            for el_name, d in self.detected_results.items()
            if d.start != -1 and el_name not in barcode_umi_set
        ]

    def _get_supplementary_elements(self) -> List[str]:
        excluded = set(self._barcode_elements) | set(self._umi_elements)
        res = []
        for el in self.molecule_structure.ordered_elements:
            if el.element_type == ElementType.cDNA:
                continue
            if el.element_name in self.molecule_structure.elements_to_concatenate:
                continue
            if el.element_name in self.molecule_structure.duplicated_elements:
                continue
            if el.element_name not in excluded and el.element_name not in res:
                res.append(el.element_name)
        return res

    def header(self) -> str:
        header = "read_id\tbarcode\tUMI\tBC_score\tvalid_UMI\tstrand"
        barcode_umi_set = set(self._barcode_elements + self._umi_elements)
        for el in self.molecule_structure:
            if el.element_type == ElementType.cDNA:
                continue
            if el.element_name in self.molecule_structure.elements_to_concatenate:
                continue
            if el.element_name in self.molecule_structure.duplicated_elements:
                continue
            if el.element_name in barcode_umi_set:
                continue
            if el.element_type == ElementType.PolyT:
                header += "\tpolyT_start\tpolyT_end"
            elif el.element_type.needs_only_coordinates():
                header += f"\t{el.element_name}_start\t{el.element_name}_end"
            elif el.element_type.needs_sequence_extraction() or el.element_type.needs_correction():
                header += f"\t{el.element_name}_start\t{el.element_name}_end\t{el.element_name}_seq"
        return header

    def __str__(self) -> str:
        parts = [
            f"{self.read_id}\t{self.get_barcode()}\t{self.get_umi()}\t{self.get_barcode_score()}\t{self.has_umi()}\t{self.strand}"
        ]
        barcode_umi_set = set(self._barcode_elements + self._umi_elements)
        for el in self.molecule_structure:
            if el.element_type == ElementType.cDNA:
                continue
            if el.element_name in self.molecule_structure.elements_to_concatenate:
                continue
            if el.element_name in self.molecule_structure.duplicated_elements:
                continue
            if el.element_name in barcode_umi_set:
                continue
            if el.element_type == ElementType.PolyT or el.element_type.is_constant():
                det = self.detected_results.get(
                    ElementType.PolyT.name if el.element_type == ElementType.PolyT else el.element_name
                )
                if det is None:
                    parts.append("\t-1\t-1")
                else:
                    parts.append(f"\t{det.start}\t{det.end}")
            elif el.element_type.needs_sequence_extraction() or el.element_type.needs_correction():
                det = self.detected_results.get(el.element_name)
                if det is None:
                    parts.append(f"\t-1\t-1\t{self.NOSEQ}")
                else:
                    seq = det.seq if det.seq else self.NOSEQ
                    parts.append(f"\t{det.start}\t{det.end}\t{seq}")
        return "".join(parts)


class ExtractionReadStats:
    """Statistics tracker matching upstream ExtractionResult format."""

    def __init__(self) -> None:
        self.read_count: int = 0
        self.bc_count: int = 0
        self.umi_count: int = 0
        self.pattern_counts: Dict[str, int] = defaultdict(int)

    def add_read(self, barcode_detection_result: ExtractionResult) -> None:
        self.read_count += 1
        if barcode_detection_result.has_barcode():
            self.bc_count += 1
        if barcode_detection_result.has_umi():
            self.umi_count += 1
        for el, det in barcode_detection_result.detected_results.items():
            if det.start != -1:
                self.pattern_counts[el] += 1

    def __str__(self) -> str:
        lines = [
            f"Total reads:\t{self.read_count}",
            f"Barcode detected:\t{self.bc_count}",
            f"UMI detected:\t{self.umi_count}",
        ]
        for a, cnt in self.pattern_counts.items():
            lines.append(f"{a}:\t{cnt}")
        return "\n".join(lines) + "\n"


class ReadStats:
    def __init__(self) -> None:
        self.read_count: int = 0
        self.bc_count: int = 0
        self.umi_count: int = 0
        self.additional_attributes_counts: Dict[str, int] = defaultdict(int)

    def add_read(self, barcode_detection_result: BarcodeResult) -> None:
        self.read_count += 1
        if hasattr(barcode_detection_result, "get_additional_attributes"):
            for a in barcode_detection_result.get_additional_attributes():
                self.additional_attributes_counts[a] += 1
        if barcode_detection_result.has_barcode():
            self.bc_count += 1
        if barcode_detection_result.has_umi():
            self.umi_count += 1

    def add_custom_stats(self, stat_name: str, val: int) -> None:
        self.additional_attributes_counts[stat_name] += val

    def __str__(self) -> str:
        out = [
            f"Total reads\t{self.read_count}",
            f"Barcode detected\t{self.bc_count}",
            f"Reliable UMI\t{self.umi_count}",
        ]
        for a, cnt in self.additional_attributes_counts.items():
            out.append(f"{a}\t{cnt}")
        return "\n".join(out) + "\n"

    def __iter__(self) -> Iterable[str]:
        yield f"Total reads: {self.read_count}"
        yield f"Barcode detected: {self.bc_count}"
        yield f"Reliable UMI: {self.umi_count}"
        for a, cnt in self.additional_attributes_counts.items():
            yield f"{a}: {cnt}"
