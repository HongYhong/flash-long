"""Custom single-cell molecule structure definitions (MDF) for UniversalSingleMoleculeExtractor."""

from __future__ import annotations

import logging
import sys
from collections import defaultdict
from enum import Enum
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

from isoquant_next.barcodes.common import load_barcodes
from isoquant_next.core.error_codes import IsoQuantExitCode

logger = logging.getLogger("IsoQuant")


class ElementType(Enum):
    PolyT = 1
    cDNA = 2
    CONST = 10
    VAR_ANY = 20
    VAR_LIST = 21
    VAR_FILE = 22
    VAR_ANY_SEPARATOR = 31
    VAR_ANY_NON_T_SEPARATOR = 32

    def is_base_separator(self) -> bool:
        return self in (ElementType.VAR_ANY_SEPARATOR, ElementType.VAR_ANY_NON_T_SEPARATOR)

    def is_variable(self) -> bool:
        return self in (ElementType.VAR_ANY, ElementType.VAR_LIST, ElementType.VAR_FILE)

    def needs_sequence_extraction(self) -> bool:
        return self in (ElementType.VAR_ANY, ElementType.VAR_ANY_SEPARATOR, ElementType.VAR_ANY_NON_T_SEPARATOR)

    def needs_correction(self) -> bool:
        return self in (ElementType.VAR_FILE, ElementType.VAR_LIST)

    def is_constant(self) -> bool:
        return self == ElementType.CONST

    def needs_only_coordinates(self) -> bool:
        return self in (ElementType.PolyT, ElementType.cDNA, ElementType.CONST)


class MoleculeElement:
    def __init__(
        self,
        element_name: str,
        element_type: ElementType,
        element_value_1: Optional[str] = None,
        element_value_2: Optional[str] = None,
    ) -> None:
        self.element_name = element_name
        self.element_type = element_type

        if element_type in (ElementType.PolyT, ElementType.cDNA):
            self.element_value = None
            self.element_length = -1
        elif element_type == ElementType.CONST:
            self.element_value = element_value_1 or ""
            self.element_length = len(self.element_value)
        elif element_type == ElementType.VAR_FILE:
            assert element_value_1 is not None
            self.element_value = load_barcodes(element_value_1)
            self.element_length = (
                len(self.element_value[0]) if element_value_2 is None else self._parse_int(element_value_2)
            )
        elif element_type == ElementType.VAR_LIST:
            assert element_value_1 is not None
            self.element_value = element_value_1.split(",")
            self.element_length = (
                len(self.element_value[0]) if element_value_2 is None else self._parse_int(element_value_2)
            )
        elif element_type in (
            ElementType.VAR_ANY,
            ElementType.VAR_ANY_SEPARATOR,
            ElementType.VAR_ANY_NON_T_SEPARATOR,
        ):
            self.element_value = None
            assert element_value_1 is not None
            self.element_length = self._parse_int(element_value_1)
        else:
            logger.critical("Wrong element type %s", element_type)
            sys.exit(int(IsoQuantExitCode.INVALID_FILE_FORMAT))

    def _parse_int(self, val: str) -> int:
        try:
            return int(val)
        except ValueError:
            logger.critical("Incorrectly specified length for element %s: %s", self.element_name, val)
            sys.exit(int(IsoQuantExitCode.INVALID_FILE_FORMAT))


class MoleculeStructure:
    CONCAT_DELIM = "|"
    DUPL_DELIM = "/"

    def __init__(self, str_iterator: Iterator[str]) -> None:
        self.ordered_elements: List[MoleculeElement] = []
        self.elements_to_concatenate: Dict[str, Tuple[str, int]] = {}
        self.concatenated_elements_counts: Dict[str, int] = {}
        self.duplicated_elements: Dict[str, Tuple[str, int]] = {}
        self.duplicated_elements_counts: Dict[str, int] = {}
        self._barcode_elements: List[str] = []
        self._umi_elements: List[str] = []

        elements, element_properties = self._parse_element_properties(str_iterator)
        for element_name in elements:
            element_type = self._make_element(element_name, element_properties)
            if self.CONCAT_DELIM in element_name:
                self.elements_to_concatenate[element_name] = self._parse_linked_element(
                    element_name, element_type, self.CONCAT_DELIM, "concatenated"
                )
            elif self.DUPL_DELIM in element_name:
                self.duplicated_elements[element_name] = self._parse_linked_element(
                    element_name, element_type, self.DUPL_DELIM, "duplicated"
                )

        self._check_linked_elements(self.elements_to_concatenate)
        self.concatenated_elements_counts = self._check_indices(self.elements_to_concatenate)
        self._check_linked_elements(self.duplicated_elements)
        self.duplicated_elements_counts = self._check_indices(self.duplicated_elements)
        self._identify_barcode_umi_elements()

    @staticmethod
    def _parse_element_properties(str_iterator: Iterator[str]):
        line = next(str_iterator)
        elements = [x.strip() for x in line.strip().split(":")]
        element_properties: Dict[str, Tuple[str, Optional[str], Optional[str]]] = {}
        for raw_line in str_iterator:
            v = raw_line.strip().split()
            if not v:
                continue
            element_name = v[0]
            if element_name in element_properties:
                logger.critical("Duplicated element name %s", element_name)
                sys.exit(int(IsoQuantExitCode.INVALID_FILE_FORMAT))
            if len(v) == 3:
                element_properties[element_name] = (v[1], v[2], None)
            elif len(v) == 4:
                element_properties[element_name] = (v[1], v[2], v[3])
            else:
                logger.critical("Incorrect number of properties in line %s", raw_line.strip())
                sys.exit(int(IsoQuantExitCode.INVALID_FILE_FORMAT))
        return elements, element_properties

    def _make_element(self, element_name: str, element_properties: dict) -> ElementType:
        if element_name not in element_properties:
            if element_name not in (ElementType.cDNA.name, ElementType.PolyT.name):
                logger.critical("Molecule element %s was not described in the format file", element_name)
                sys.exit(int(IsoQuantExitCode.INVALID_FILE_FORMAT))
            etype = ElementType[element_name]
            self.ordered_elements.append(MoleculeElement(element_name, etype))
            return etype

        etype_name, v1, v2 = element_properties[element_name]
        if etype_name not in ElementType.__members__:
            logger.critical("Molecule element type %s is not valid", etype_name)
            sys.exit(int(IsoQuantExitCode.INVALID_FILE_FORMAT))
        etype = ElementType[etype_name]
        self.ordered_elements.append(MoleculeElement(element_name, etype, v1, v2))
        return etype

    @staticmethod
    def _parse_linked_element(element_name: str, element_type: ElementType, delim: str, kind: str) -> Tuple[str, int]:
        if not (element_type.needs_sequence_extraction() or element_type.needs_correction()):
            logger.critical("Concatenated elements must be variable (fix element %s)", element_name)
            sys.exit(int(IsoQuantExitCode.INVALID_FILE_FORMAT))
        v = element_name.split(delim)
        if len(v) != 2:
            logger.critical("Incorrect %s element %s", kind, element_name)
            sys.exit(int(IsoQuantExitCode.INVALID_FILE_FORMAT))
        try:
            return (v[0], int(v[1]))
        except ValueError:
            logger.critical("Incorrect index for %s element %s: %s", kind, element_name, v[1])
            sys.exit(int(IsoQuantExitCode.INVALID_FILE_FORMAT))

    def _check_linked_elements(self, element_dict: dict) -> None:
        linked: Dict[str, List[MoleculeElement]] = defaultdict(list)
        for el in self.ordered_elements:
            if el.element_name in element_dict:
                linked[element_dict[el.element_name][0]].append(el)
        for el_list in linked.values():
            if len({el.element_type for el in el_list}) != 1:
                sys.exit(int(IsoQuantExitCode.INVALID_FILE_FORMAT))
            first_val = el_list[0].element_value
            if any(el.element_value != first_val for el in el_list[1:]):
                sys.exit(int(IsoQuantExitCode.INVALID_FILE_FORMAT))

    def _check_indices(self, element_dict: dict) -> Dict[str, int]:
        index_dict: Dict[str, List[int]] = defaultdict(list)
        for el, (base_name, idx) in element_dict.items():
            index_dict[base_name].append(idx)
        for el, idxs in index_dict.items():
            if sorted(idxs) != list(range(1, len(idxs) + 1)):
                sys.exit(int(IsoQuantExitCode.INVALID_FILE_FORMAT))
        return {el: len(val) for el, val in index_dict.items()}

    def _identify_barcode_umi_elements(self) -> None:
        self._barcode_elements = []
        self._umi_elements = []
        seen = set()
        for el in self.ordered_elements:
            if el.element_name in self.elements_to_concatenate:
                name = self.elements_to_concatenate[el.element_name][0]
            elif el.element_name in self.duplicated_elements:
                name = self.duplicated_elements[el.element_name][0]
            else:
                name = el.element_name
            if name in seen:
                continue
            seen.add(name)
            low = name.lower()
            if low.startswith("barcode"):
                self._barcode_elements.append(name)
            elif low.startswith("umi"):
                self._umi_elements.append(name)

    @property
    def barcode_elements(self) -> List[str]:
        return self._barcode_elements

    @property
    def umi_elements(self) -> List[str]:
        return self._umi_elements

    @classmethod
    def from_element_list(cls, element_list: List[MoleculeElement]) -> "MoleculeStructure":
        ms = cls.__new__(cls)
        ms.ordered_elements = element_list
        ms.elements_to_concatenate = {}
        ms.duplicated_elements = {}
        ms.concatenated_elements_counts = {}
        ms.duplicated_elements_counts = {}
        ms._barcode_elements = []
        ms._umi_elements = []
        ms._identify_barcode_umi_elements()
        return ms

    def __iter__(self) -> Iterator[MoleculeElement]:
        return iter(self.ordered_elements)
