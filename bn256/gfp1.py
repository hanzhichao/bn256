from .constants import P
from .utils import inv_mod_p, inverse_mod


# Montgomery params
R = pow(2, 256)
R1 = R % P
R2 = (R * R) % P
R3 = (R1 * R2) % P
N = R - inverse_mod(P, R)


class Gfp1(object):
    v: int

    def __init__(self, v: int, redc_needed: bool = True):
        # assert v >= 0
        if redc_needed:
            # assert v < R
            self.v = self._redc(v * R2)
            # assert self.value() == v % p
        else:
            # assert v < p
            self.v = v

    def __repr__(self):
        return "<Gfp1 %d>" % (self.value())

    def __eq__(self, other: "Gfp1"):
        return self.v == other.v

    def __ne__(self, other: "Gfp1"):
        return self.v != other.v

    def __str__(self):
        return "%d" % (self.value())

    def __add__(self, other: "Gfp1") -> "Gfp1":
        x = self.v + other.v
        if (x >> 255) > 0 and x >= P:
            x -= P
        # assert self._redc(x) == (self.value() + other.value()) % p
        return Gfp1(x, False)

    def __sub__(self, other: "Gfp1") -> "Gfp1":
        x = self.v - other.v
        if x < 0:
            x += P
        # assert self._redc(x) == (self.value() - other.value()) % p
        return Gfp1(x, False)

    def __mul__(self, other: "Gfp1") -> "Gfp1":
        return Gfp1(self._redc(self.v * other.v), False)

    def value(self) -> int:
        return self._redc(self.v)

    @staticmethod
    def _redc(T):
        # assert T < (R*p-1)
        m = ((T & (R - 1)) * N) & (R - 1)
        t = (T + m * P) >> 256
        if t >= P:
            t -= P
        return t

    def square(self) -> "Gfp1":
        return self * self

    def double(self) -> "Gfp1":
        return self + self

    def triple(self) -> "Gfp1":
        return self + self + self

    def is_one(self) -> bool:
        return self.v == R1

    def is_zero(self) -> bool:
        return self.v == 0

    def inverse(self) -> "Gfp1":
        # Fermat
        x = Gfp1(self._redc(R3 * inv_mod_p(self.v)), False)
        # assert (self * x).value() == 1
        return x

    def additive_inverse(self):
        x = Gfp1(P, True) - self
        # assert (self + x).value() == 0
        return x

    def to_bytes(self):
        p_bytes = (P.bit_length() // 8) + (1 if (P.bit_length() % 8) > 0 else 0)
        return self.value().to_bytes(p_bytes, "big")


GFP1_ZERO = Gfp1(0)
GFP1_ONE = Gfp1(1)
CURVE_B = Gfp1(3)
