"""Universal single-cell molecule structure extractor using Cython SSW & 2-bit indexers."""

from __future__ import annotations

import logging
import math
import sys
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Set, Tuple

from isoquant_next.barcodes.molecule_structure import ElementType, MoleculeElement, MoleculeStructure
from isoquant_next.barcodes.common import (
    batch_str_to_2bit,
    detect_exact_positions,
    find_candidate_with_max_score_ssw,
    find_optimal_kmer_size,
    find_polyt,
    reverese_complement,
)
from isoquant_next.barcodes.indexers import ArrayKmerIndexer, Dict2BitKmerIndexer, KmerIndexer
from isoquant_next.barcodes.results import DetectedElement, ExtractionResult
from isoquant_next.core.error_codes import IsoQuantExitCode

logger = logging.getLogger("IsoQuant")


class UniversalSingleMoleculeExtractor:
    """Detect and optionally error-correct elements defined in a MoleculeStructure MDF spec."""

    MIN_SCORE_COEFF = 0.75
    MIN_SCORE_COEFF_TERMMINAL = 0.5
    TERMINAL_MATCH_DELTA = 3
    MAX_LEN_DIFF = 0.25
    DEFAULT_ERROR_RATE = 0.03
    MAX_HITS = 10
    SCORE_DIFF = 1

    def __init__(self, molecule_structure: MoleculeStructure, correct_sequences: bool = True) -> None:
        self.correct_sequences = correct_sequences
        self.molecule_structure = molecule_structure
        self.index_dict: Dict[str, object] = {}
        self.min_scores: Dict[str, int] = {}
        self.const_sequences: Dict[str, str] = {}
        self.has_polyt = False
        self.has_cdna = False
        self.constant_elements_to_detect: Set[str] = set()
        self.elements_to_extract: Set[str] = set()
        self.elements_to_correct: Set[str] = set()

        self._classify_elements()

        if self.correct_sequences:
            self._build_concatenated_indexes()
            self._build_duplicated_indexes()

        if not self.has_cdna:
            logger.critical("Molecule must include a cDNA")
            sys.exit(int(IsoQuantExitCode.INVALID_PARAMETER))

    def _classify_elements(self) -> None:
        for el in self.molecule_structure:
            if el.element_type == ElementType.PolyT:
                if self.has_polyt:
                    logger.critical(
                        "Current version only supports a single polyT, extraction results may be suboptimal"
                    )
                    sys.exit(int(IsoQuantExitCode.INVALID_PARAMETER))
                self.has_polyt = True
            elif el.element_type == ElementType.cDNA:
                if self.has_cdna:
                    logger.critical(
                        "Current version only supports a single cDNA, reads will not be split into molecules"
                    )
                    sys.exit(int(IsoQuantExitCode.INVALID_PARAMETER))
                self.has_cdna = True
            elif el.element_type.is_constant():
                self.constant_elements_to_detect.add(el.element_name)
                self.const_sequences[el.element_name] = el.element_value or ""
                k = max(3, min(el.element_length, max(6, int(el.element_length / 2) - 2)))
                self.index_dict[el.element_name] = KmerIndexer([el.element_value], k)
                self.min_scores[el.element_name] = max(1, int(math.floor(self.MIN_SCORE_COEFF * el.element_length)))
            elif el.element_type.needs_sequence_extraction():
                self.elements_to_extract.add(el.element_name)
            elif el.element_type.needs_correction():
                self._handle_correction_element(el)

    def _handle_correction_element(self, el: MoleculeElement) -> None:
        if el.element_name in self.molecule_structure.elements_to_concatenate:
            return
        if el.element_name in self.molecule_structure.duplicated_elements:
            return
        if self.correct_sequences:
            self.prepare_barcode_index_for_element(el)
            self.elements_to_correct.add(el.element_name)
        else:
            self.elements_to_extract.add(el.element_name)

    def _build_concatenated_indexes(self) -> None:
        for base_name, _count in self.molecule_structure.concatenated_elements_counts.items():
            first_part_el: Optional[MoleculeElement] = None
            total_length = 0
            for el in self.molecule_structure:
                info = self.molecule_structure.elements_to_concatenate.get(el.element_name)
                if info and info[0] == base_name:
                    if first_part_el is None:
                        first_part_el = el
                    total_length += el.element_length
            if first_part_el is not None and isinstance(first_part_el.element_value, list):
                self.prepare_barcode_index(base_name, first_part_el.element_value, total_length)
                self.elements_to_correct.add(base_name)

    def _build_duplicated_indexes(self) -> None:
        for base_name, _count in self.molecule_structure.duplicated_elements_counts.items():
            first_part_el: Optional[MoleculeElement] = None
            for el in self.molecule_structure:
                info = self.molecule_structure.duplicated_elements.get(el.element_name)
                if info and info[0] == base_name:
                    first_part_el = el
                    break
            if first_part_el is not None and isinstance(first_part_el.element_value, list):
                self.prepare_barcode_index(base_name, first_part_el.element_value, first_part_el.element_length)
                self.elements_to_correct.add(base_name)

    def prepare_barcode_index_for_element(self, element: MoleculeElement) -> None:
        assert isinstance(element.element_value, list)
        self.prepare_barcode_index(element.element_name, element.element_value, element.element_length)

    def prepare_barcode_index(self, base_name: str, barcode_list: List[str], barcode_length: int) -> None:
        kmer_size = find_optimal_kmer_size(barcode_length, len(barcode_list))
        k = max(2, min(kmer_size, barcode_length))
        if len(barcode_list) > 100_000 and barcode_length <= 32:
            self.index_dict[base_name] = Dict2BitKmerIndexer(
                batch_str_to_2bit(barcode_list, barcode_length), kmer_size=k, seq_len=barcode_length
            )
        elif k <= 8 and all(set(b.upper()) <= {"A", "C", "G", "T"} for b in barcode_list[:64]):
            self.index_dict[base_name] = ArrayKmerIndexer(barcode_list, kmer_size=k)
        else:
            self.index_dict[base_name] = KmerIndexer(barcode_list, kmer_size=k)
        self.min_scores[base_name] = max(1, int(math.floor(self.MIN_SCORE_COEFF * barcode_length)))

    def result_type(self):
        return ExtractionResult

    def header(self) -> str:
        return ExtractionResult(self.molecule_structure, "", ".", {}).header()

    def detect_const_elements(
        self,
        sequence: str,
        polyt_start: int,
        polyt_end: int,
        detected_elements: Dict[str, DetectedElement],
    ) -> None:
        n = len(sequence)
        for el_name in self.constant_elements_to_detect:
            indexer = self.index_dict[el_name]
            pattern = self.const_sequences[el_name]
            matches = indexer.get_occurrences(sequence, max_hits=self.MAX_HITS)
            if not matches:
                continue
            min_score = self.min_scores[el_name]
            st, en = detect_exact_positions(sequence, 0, n, indexer.k, pattern, matches, min_score=min_score)
            if st is None or en is None:
                # Allow relaxed terminal constant detection near read ends
                if matches[0][1] > 0:
                    st, en = detect_exact_positions(
                        sequence,
                        0,
                        n,
                        indexer.k,
                        pattern,
                        matches,
                        min_score=max(1, int(math.floor(self.MIN_SCORE_COEFF_TERMMINAL * len(pattern)))),
                    )
            if st is not None and en is not None:
                detected_elements[el_name] = DetectedElement(st, en, len(pattern))

    def _detect_potential_start(
        self,
        element: MoleculeElement,
        element_index: int,
        potential_end: int,
        detected_elements: Dict[str, DetectedElement],
    ) -> Optional[int]:
        if element_index > 0:
            prev_el = self.molecule_structure.ordered_elements[element_index - 1]
            if prev_el.element_name in detected_elements:
                det_prev = detected_elements[prev_el.element_name]
                if det_prev.end != -1:
                    return max(0, det_prev.end + 1)
        if potential_end != -1 and element.element_length > 0:
            return max(0, potential_end - element.element_length + 1)
        return None

    def _detect_potential_end(
        self,
        sequence: str,
        element: MoleculeElement,
        element_index: int,
        potential_start: int,
        detected_elements: Dict[str, DetectedElement],
    ) -> Optional[int]:
        ordered = self.molecule_structure.ordered_elements
        if element_index + 1 < len(ordered):
            next_el = ordered[element_index + 1]
            if next_el.element_name in detected_elements:
                det_next = detected_elements[next_el.element_name]
                if det_next.start != -1:
                    return min(len(sequence) - 1, det_next.start - 1)
        if potential_start != -1 and element.element_length > 0:
            return min(len(sequence) - 1, potential_start + element.element_length - 1)
        return None

    def correct_element(
        self,
        element: MoleculeElement,
        sequence: str,
        potential_start: int,
        potential_end: int,
    ) -> DetectedElement:
        if potential_start < 0 or potential_end < potential_start or potential_start >= len(sequence):
            return DetectedElement(-1, -1, -1, ExtractionResult.NOSEQ)
        clamped_end = min(len(sequence) - 1, potential_end)
        # Pad window by +-2 bp to absorb indels near element boundaries
        win_start = max(0, potential_start - 2)
        win_end = min(len(sequence), clamped_end + 3)
        subseq = sequence[win_start:win_end]

        base_name = self._base_name(element.element_name)
        indexer = self.index_dict.get(base_name)
        if indexer is None:
            return DetectedElement(-1, -1, -1, ExtractionResult.NOSEQ)
        min_kmers = 1 if element.element_length < 12 else max(1, element.element_length // max(2, indexer.k) - 2)
        matches = indexer.get_occurrences(subseq, max_hits=self.MAX_HITS, min_kmers=min_kmers)
        if not matches:
            return DetectedElement(-1, -1, -1, ExtractionResult.NOSEQ)

        bc, score, rel_st, rel_en = find_candidate_with_max_score_ssw(
            matches,
            subseq,
            min_score=self.min_scores[base_name],
            score_diff=self.SCORE_DIFF,
        )
        if not bc:
            return DetectedElement(-1, -1, -1, ExtractionResult.NOSEQ)
        abs_st = max(0, win_start + rel_st)
        abs_en = min(len(sequence) - 1, win_start + rel_en - 1)
        return DetectedElement(abs_st, abs_en, score, bc)

    def _base_name(self, element_name: str) -> str:
        if element_name in self.molecule_structure.elements_to_concatenate:
            return self.molecule_structure.elements_to_concatenate[element_name][0]
        if element_name in self.molecule_structure.duplicated_elements:
            return self.molecule_structure.duplicated_elements[element_name][0]
        return element_name

    def _backup_concatenated_element(
        self,
        element: MoleculeElement,
        potential_start: int,
        potential_end: int,
        element_storage: Dict[str, List[Tuple[int, int]]],
    ) -> None:
        base_name, part_idx = self.molecule_structure.elements_to_concatenate[element.element_name]
        parts = element_storage[base_name]
        while len(parts) < part_idx:
            parts.append((-1, -1))
        parts[part_idx - 1] = (potential_start, potential_end)

    def _backup_duplicated_element(
        self,
        element: MoleculeElement,
        potential_start: int,
        potential_end: int,
        element_storage: Dict[str, List[Tuple[int, int]]],
    ) -> None:
        base_name, part_idx = self.molecule_structure.duplicated_elements[element.element_name]
        parts = element_storage[base_name]
        while len(parts) < part_idx:
            parts.append((-1, -1))
        parts[part_idx - 1] = (potential_start, potential_end)

    def _detect_variable_element(
        self,
        el: MoleculeElement,
        sequence: str,
        potential_start: int,
        potential_end: int,
        detected_elements: Dict[str, DetectedElement],
        concat_element_storage: Dict[str, List[Tuple[int, int]]],
        dupl_element_storage: Dict[str, List[Tuple[int, int]]],
        len_diff_threshold: Optional[int] = None,
    ) -> bool:
        if potential_start < 0 or potential_end < potential_start or potential_end >= len(sequence):
            return False
        obs_len = potential_end - potential_start + 1
        tol = (
            len_diff_threshold
            if len_diff_threshold is not None
            else max(2, int(math.ceil(el.element_length * self.MAX_LEN_DIFF)))
        )
        if el.element_type.is_base_separator():
            if el.element_type == ElementType.VAR_ANY_NON_T_SEPARATOR:
                # Trim trailing T homopolymer adjacent to polyT
                while potential_end >= potential_start and sequence[potential_end] == "T":
                    potential_end -= 1
                obs_len = max(0, potential_end - potential_start + 1)
            if obs_len <= 0:
                return False
            detected_elements[el.element_name] = DetectedElement(
                potential_start, potential_end, 0, sequence[potential_start : potential_end + 1]
            )
            return True

        if abs(obs_len - el.element_length) > tol:
            # Clamp to nominal expected length when only one boundary is firm
            if potential_end - potential_start + 1 > el.element_length + tol:
                return False
            return False

        if el.element_name in self.molecule_structure.elements_to_concatenate:
            self._backup_concatenated_element(el, potential_start, potential_end, concat_element_storage)
            detected_elements[el.element_name] = DetectedElement(
                potential_start, potential_end, 0, sequence[potential_start : potential_end + 1]
            )
            return True

        if el.element_name in self.molecule_structure.duplicated_elements:
            self._backup_duplicated_element(el, potential_start, potential_end, dupl_element_storage)
            detected_elements[el.element_name] = DetectedElement(
                potential_start, potential_end, 0, sequence[potential_start : potential_end + 1]
            )
            return True

        if el.element_name in self.elements_to_correct:
            det = self.correct_element(el, sequence, potential_start, potential_end)
            if det.start != -1:
                detected_elements[el.element_name] = det
                return True
            return False

        detected_elements[el.element_name] = DetectedElement(
            potential_start, potential_end, 0, sequence[potential_start : potential_end + 1]
        )
        return True

    def extract_variable_elements3(
        self,
        sequence: str,
        detected_elements: Dict[str, DetectedElement],
        concat_element_storage: Dict[str, List[Tuple[int, int]]],
        dupl_element_storage: Dict[str, List[Tuple[int, int]]],
    ) -> None:
        ordered = self.molecule_structure.ordered_elements
        for idx, el in enumerate(ordered):
            if el.element_name in detected_elements:
                continue
            if not (el.element_type.needs_sequence_extraction() or el.element_type.needs_correction()):
                continue
            pot_start = None
            if idx > 0 and ordered[idx - 1].element_name in detected_elements:
                prev_det = detected_elements[ordered[idx - 1].element_name]
                if prev_det.end != -1:
                    pot_start = prev_det.end + 1
            if pot_start is None:
                continue
            pot_end = self._detect_potential_end(sequence, el, idx, pot_start, detected_elements)
            if pot_end is None:
                pot_end = min(len(sequence) - 1, pot_start + el.element_length - 1)
            self._detect_variable_element(
                el, sequence, pot_start, pot_end, detected_elements, concat_element_storage, dupl_element_storage
            )

    def extract_variable_elements5(
        self,
        sequence: str,
        detected_elements: Dict[str, DetectedElement],
        concat_element_storage: Dict[str, List[Tuple[int, int]]],
        dupl_element_storage: Dict[str, List[Tuple[int, int]]],
    ) -> None:
        ordered = self.molecule_structure.ordered_elements
        for idx in range(len(ordered) - 1, -1, -1):
            el = ordered[idx]
            if el.element_name in detected_elements:
                continue
            if not (el.element_type.needs_sequence_extraction() or el.element_type.needs_correction()):
                continue
            pot_end = None
            if idx + 1 < len(ordered) and ordered[idx + 1].element_name in detected_elements:
                next_det = detected_elements[ordered[idx + 1].element_name]
                if next_det.start != -1:
                    pot_end = next_det.start - 1
            if pot_end is None:
                continue
            pot_start = self._detect_potential_start(el, idx, pot_end, detected_elements)
            if pot_start is None:
                pot_start = max(0, pot_end - el.element_length + 1)
            self._detect_variable_element(
                el, sequence, pot_start, pot_end, detected_elements, concat_element_storage, dupl_element_storage
            )

    def process_concatenated_elements(
        self,
        concatenated_element_storage: Dict[str, List[Tuple[int, int]]],
        detected_results: Dict[str, DetectedElement],
        sequence: str,
    ) -> None:
        for base_name, coords in concatenated_element_storage.items():
            if not coords or any(st == -1 or en == -1 for st, en in coords):
                detected_results[base_name] = DetectedElement(-1, -1, -1, ExtractionResult.NOSEQ)
                continue
            raw_seq = "".join(sequence[st : en + 1] for st, en in coords)
            first_st = coords[0][0]
            last_en = coords[-1][1]
            if not self.correct_sequences:
                detected_results[base_name] = DetectedElement(first_st, last_en, 0, raw_seq)
                continue
            indexer = self.index_dict.get(base_name)
            if indexer is None:
                detected_results[base_name] = DetectedElement(-1, -1, -1, ExtractionResult.NOSEQ)
                continue
            matches = indexer.get_occurrences(raw_seq, max_hits=self.MAX_HITS, min_kmers=1)
            bc, score, _rel_st, _rel_en = find_candidate_with_max_score_ssw(
                matches, raw_seq, min_score=self.min_scores[base_name], score_diff=self.SCORE_DIFF
            )
            if bc:
                detected_results[base_name] = DetectedElement(first_st, last_en, score, bc)
            else:
                detected_results[base_name] = DetectedElement(-1, -1, -1, ExtractionResult.NOSEQ)

    def process_duplicated_elements(
        self,
        duplicated_element_storage: Dict[str, List[Tuple[int, int]]],
        detected_results: Dict[str, DetectedElement],
        sequence: str,
    ) -> None:
        for base_name, coords in duplicated_element_storage.items():
            valid_coords = [(st, en) for st, en in coords if st != -1 and en != -1]
            if not valid_coords:
                detected_results[base_name] = DetectedElement(-1, -1, -1, ExtractionResult.NOSEQ)
                continue

            if not self.correct_sequences:
                raw_seqs = [sequence[st : en + 1] for st, en in valid_coords]
                if len(set(raw_seqs)) == 1:
                    seq_out = raw_seqs[0]
                else:
                    seq_out = ",".join(raw_seqs)
                detected_results[base_name] = DetectedElement(valid_coords[0][0], valid_coords[0][1], 0, seq_out)
                continue

            indexer = self.index_dict.get(base_name)
            if indexer is None:
                detected_results[base_name] = DetectedElement(-1, -1, -1, ExtractionResult.NOSEQ)
                continue
            min_score = self.min_scores[base_name]
            corrected: List[Tuple[str, int, int, int]] = []
            for st, en in valid_coords:
                sub = sequence[st : en + 1]
                matches = indexer.get_occurrences(sub, max_hits=self.MAX_HITS, min_kmers=1)
                bc, score, _rs, _re = find_candidate_with_max_score_ssw(
                    matches, sub, min_score=min_score, score_diff=self.SCORE_DIFF
                )
                if bc:
                    corrected.append((bc, score, st, en))

            if not corrected:
                detected_results[base_name] = DetectedElement(-1, -1, -1, ExtractionResult.NOSEQ)
                continue
            counts = Counter(c[0] for c in corrected).most_common()
            if len(counts) > 1 and counts[0][1] == counts[1][1]:
                # Tied disagreement among duplicated copies cannot be safely resolved
                detected_results[base_name] = DetectedElement(-1, -1, -1, ExtractionResult.NOSEQ)
                continue
            winner_bc = counts[0][0]
            best_copy = max((c for c in corrected if c[0] == winner_bc), key=lambda x: x[1])
            detected_results[base_name] = DetectedElement(best_copy[2], best_copy[3], best_copy[1], winner_bc)

    def _find_patterns_fwd(self, read_id: str, sequence: str) -> ExtractionResult:
        detected_elements: Dict[str, DetectedElement] = {}
        polyt_start, polyt_end = (-1, -1)
        if self.has_polyt:
            polyt_start, polyt_end = find_polyt(sequence)
            if polyt_start != -1:
                detected_elements[ElementType.PolyT.name] = DetectedElement(polyt_start, polyt_end, 0)

        self.detect_const_elements(sequence, polyt_start, polyt_end, detected_elements)

        concat_storage: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
        dupl_storage: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
        # Iterate alternating 5' and 3' extensions until fixed point so chained variable elements resolve
        prev_len = -1
        while len(detected_elements) != prev_len:
            prev_len = len(detected_elements)
            self.extract_variable_elements3(sequence, detected_elements, concat_storage, dupl_storage)
            self.extract_variable_elements5(sequence, detected_elements, concat_element_storage=concat_storage, dupl_element_storage=dupl_storage)

        if concat_storage:
            self.process_concatenated_elements(concat_storage, detected_elements, sequence)
        if dupl_storage:
            self.process_duplicated_elements(dupl_storage, detected_elements, sequence)

        return ExtractionResult(self.molecule_structure, read_id, ".", detected_elements)

    def find_barcode_umi(self, read_id: str, sequence: str) -> ExtractionResult:
        fwd_res = self._find_patterns_fwd(read_id, sequence)
        rev_seq = reverese_complement(sequence)
        rev_res = self._find_patterns_fwd(read_id, rev_seq)
        if rev_res.more_informative_than(fwd_res):
            rev_res.set_strand("-")
            return rev_res
        if fwd_res.detected_results:
            fwd_res.set_strand("+")
        return fwd_res
