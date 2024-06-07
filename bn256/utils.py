import os
import re
import sys

from .constants import P, ORDER


def bits_of(k):
    return [int(c) for c in "{0:b}".format(k)]


def inverse_mod(a, n):
    t = 0
    t2 = 1
    r = n
    r2 = a

    while r2 != 0:
        q = r // r2
        (t, t2) = (t2, t - q * t2)
        (r, r2) = (r2, r - q * r2)

    if r > 1:
        return 0
    if t < 0:
        t += n
    return t


def inv_mod_p(a, p=P):
    # Fermat
    return pow(a, p - 2, p)


def legendre(a, p=P):
    x = pow(a, (p - 1) // 2, p)
    if x == 0 or x == 1:
        return x
    if x == p - 1:
        return -1
    assert False


def rand_elem() -> int:
    """Debiased random element generator"""
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


def split_by_length(data: str, length: int) -> list:
    assert len(data) % length == 0
    return re.findall(r".{%d}" % length, data)
