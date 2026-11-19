"""Binary serialization utilities using struct packing with signed-int dict fix."""

from __future__ import annotations

import struct
from typing import Any, BinaryIO, Callable, Dict, List, Optional, Sequence, Tuple

ENCODING = "utf-8"
STR_LEN_BYTES = 4
NONE_STR_LEN = 0xFFFFFFFF
SHORT_INT_BYTES = 2
LONG_INT_BYTES = 4
TERMINATION_INT = 0xFFFFFFFF
SHORT_TERMINATION_INT = 0xFFFF
SHORT_FLOAT_MULTIPLIER = 1 << 20
DICT_TYPE_LEN = 1
DICT_INT_TYPE = 9
DICT_INT_PAIR_TYPE = 10
DICT_STR_TYPE = 17

_U8 = struct.Struct(">B")
_U16 = struct.Struct(">H")
_U32 = struct.Struct(">I")
_SIGN_BIT = 1 << 31
_SIGN_MASK = _SIGN_BIT - 1


def write_string(s: str, outf: BinaryIO) -> None:
    encoded = s.encode(ENCODING)
    outf.write(_U32.pack(len(encoded)) + encoded)


def read_string(inf: BinaryIO) -> str:
    (str_len,) = _U32.unpack(inf.read(4))
    return inf.read(str_len).decode(ENCODING) if str_len else ""


def write_string_or_none(s: Optional[str], outf: BinaryIO) -> None:
    if s is None:
        outf.write(_U32.pack(NONE_STR_LEN))
        return
    write_string(s, outf)


def read_string_or_none(inf: BinaryIO) -> Optional[str]:
    (str_len,) = _U32.unpack(inf.read(4))
    if str_len == NONE_STR_LEN:
        return None
    return inf.read(str_len).decode(ENCODING) if str_len else ""


def write_int(val: int, outf: BinaryIO, bytes_len: int = LONG_INT_BYTES) -> None:
    if bytes_len == 4:
        outf.write(_U32.pack(val))
    elif bytes_len == 2:
        outf.write(_U16.pack(val))
    elif bytes_len == 1:
        outf.write(_U8.pack(val))
    else:
        outf.write(int(val).to_bytes(bytes_len, "big"))


def read_int(inf: BinaryIO, bytes_len: int = LONG_INT_BYTES) -> int:
    if bytes_len == 4:
        return _U32.unpack(inf.read(4))[0]
    if bytes_len == 2:
        return _U16.unpack(inf.read(2))[0]
    if bytes_len == 1:
        return _U8.unpack(inf.read(1))[0]
    return int.from_bytes(inf.read(bytes_len), "big")


def write_int_or_none(val: Optional[int], outf: BinaryIO, bytes_len: int = LONG_INT_BYTES) -> None:
    if val is None:
        outf.write(TERMINATION_INT.to_bytes(bytes_len, "big") if bytes_len != 4 else _U32.pack(TERMINATION_INT))
        return
    write_int(val, outf, bytes_len)


def read_int_or_none(inf: BinaryIO, bytes_len: int = LONG_INT_BYTES) -> Optional[int]:
    val = read_int(inf, bytes_len)
    if val == TERMINATION_INT:
        return None
    return val


def write_short_int(val: int, outf: BinaryIO) -> None:
    outf.write(_U16.pack(val))


def read_short_int(inf: BinaryIO) -> int:
    return _U16.unpack(inf.read(2))[0]


def write_int_neg(val: int, outf: BinaryIO) -> None:
    """Encode signed 31-bit integer using bit 31 as sign flag."""
    if val < 0:
        mag = -val
        if mag & _SIGN_BIT:
            raise OverflowError(f"Integer out of 31-bit signed range: {val}")
        outf.write(_U32.pack(mag | _SIGN_BIT))
    else:
        if val & _SIGN_BIT:
            raise OverflowError(f"Integer out of 31-bit signed range: {val}")
        outf.write(_U32.pack(val))


def read_int_neg(inf: BinaryIO) -> int:
    (raw,) = _U32.unpack(inf.read(4))
    if raw & _SIGN_BIT:
        return -(raw & _SIGN_MASK)
    return raw


def write_list(lst: Sequence[Any], outf: BinaryIO, func: Callable[[Any, BinaryIO], None]) -> None:
    outf.write(_U32.pack(len(lst)))
    for val in lst:
        func(val, outf)


def read_list(inf: BinaryIO, func: Callable[[BinaryIO], Any]) -> List[Any]:
    (list_size,) = _U32.unpack(inf.read(4))
    return [func(inf) for _ in range(list_size)]


def write_list_of_pairs(
    lst: Sequence[Tuple[Any, Any]], outf: BinaryIO, func: Callable[[Any, BinaryIO], None]
) -> None:
    outf.write(_U32.pack(len(lst)))
    for a, b in lst:
        func(a, outf)
        func(b, outf)


def read_list_of_pairs(inf: BinaryIO, func: Callable[[BinaryIO], Any]) -> List[Tuple[Any, Any]]:
    (list_size,) = _U32.unpack(inf.read(4))
    return [(func(inf), func(inf)) for _ in range(list_size)]


def write_bool_array(bool_arr: Sequence[bool], outf: BinaryIO) -> None:
    if len(bool_arr) > 8:
        raise ValueError("write_bool_array supports at most 8 booleans per byte")
    byte_val = 0
    for i, val in enumerate(bool_arr):
        if val:
            byte_val |= 1 << i
    outf.write(_U8.pack(byte_val))


def read_bool_array(inf: BinaryIO, arr_size: int = 1) -> List[bool]:
    (byte_val,) = _U8.unpack(inf.read(1))
    return [(byte_val & (1 << i)) != 0 for i in range(arr_size)]


def write_dict(d: Dict[str, Any], outf: BinaryIO) -> None:
    outf.write(_U32.pack(len(d)))
    for k, v in d.items():
        write_string(k, outf)
        if isinstance(v, int):
            outf.write(_U8.pack(DICT_INT_TYPE))
            write_int_neg(v, outf)
        elif isinstance(v, str):
            outf.write(_U8.pack(DICT_STR_TYPE))
            write_string(v, outf)
        elif isinstance(v, tuple) and len(v) == 2:
            outf.write(_U8.pack(DICT_INT_PAIR_TYPE))
            write_int_neg(v[0], outf)
            write_int_neg(v[1], outf)
        else:
            raise ValueError(f"Unsupported value type in dictionary serialization: {type(v).__name__}")


def read_dict(inf: BinaryIO) -> Dict[str, Any]:
    """Read serialized dictionary with fixed signed-integer decoding (`read_int_neg`)."""
    (d_len,) = _U32.unpack(inf.read(4))
    d: Dict[str, Any] = {}
    for _ in range(d_len):
        k = read_string(inf)
        (vtype,) = _U8.unpack(inf.read(1))
        if vtype == DICT_INT_TYPE:
            d[k] = read_int_neg(inf)
        elif vtype == DICT_STR_TYPE:
            d[k] = read_string(inf)
        elif vtype == DICT_INT_PAIR_TYPE:
            d[k] = (read_int_neg(inf), read_int_neg(inf))
        else:
            raise ValueError(f"Serialized dictionary contains unsupported value type {vtype}")
    return d
