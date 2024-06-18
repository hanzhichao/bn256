import os
import sys
from typing import List

from .constants import P, ORDER

BYTE_LEN = 32


def bits_of(k):
    return [int(c) for c in "{0:b}".format(k)]


def sqrt_mod_p(a: int, p: int):
    assert p % 4 == 3
    return pow(a, (p + 1) // 4, p)


def mod_inverse(a: int, p: int):
    # Fermat
    # as golang big.Int.ModInverse()
    return pow(a, p - 2, p)


#
def random_k() -> int:
    rand_elem_bytes = (ORDER.bit_length() + 7) // 8 + 1
    rand_elem_base = 2
    rand_elem_range = ORDER - rand_elem_base
    rand_elem_barrier = (1 << (8 * rand_elem_bytes)) - rand_elem_range

    while True:
        rand_bytes = os.urandom(rand_elem_bytes)
        rand_num = int.from_bytes(rand_bytes, sys.byteorder)
        res = rand_num % rand_elem_range
        if (rand_num - res) <= rand_elem_barrier:
            return res + rand_elem_base


def nums_to_bytes(*args: int) -> bytes:
    return b''.join([i.to_bytes(BYTE_LEN, 'big') for i in args])


def bytes_to_nums(data: bytes, cnt: int) -> List[int]:
    assert isinstance(data, bytes)
    assert len(data) == BYTE_LEN * cnt
    return [int.from_bytes(data[BYTE_LEN * i: BYTE_LEN * (i + 1)], "big")
            for i in range(cnt)]
