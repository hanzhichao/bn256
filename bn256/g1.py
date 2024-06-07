import hashlib
from typing import List

from .point import CurvePoint, CURVE_G
from .gfp1 import P, GFP1_ONE, Gfp1, CURVE_B
from .utils import inv_mod_p, legendre, rand_elem, split_by_length


def sqrt_mod_p(a):
    assert P % 4 == 3
    return pow(a, (P + 1) // 4, P)


class G1:
    def __init__(self, p: CurvePoint):
        self.p = p

    def __repr__(self):
        return "<G1 p=%s>" % self.p

    def __add__(self, other: "G1") -> "G1":
        return G1(self.p + other.p)

    def add(self, other: "G1") -> "G1":
        return self.__add__(other)

    def marshal(self) -> bytes:
        self.p.force_affine()
        return b"".join(
            map(
                lambda x: x.to_bytes(32, byteorder="big"),
                [self.p.x.value(), self.p.y.value()],
            ),
        )

    @classmethod
    def unmarshal(cls, data: bytes) -> "G1":
        data_hex = data.hex()
        assert len(data_hex) == 128
        x_hex, y_hex = split_by_length(data_hex, 64)
        x, y = int(x_hex, 16), int(y_hex, 16)
        return G1(CurvePoint(Gfp1(x), Gfp1(y)))

    @classmethod
    def random_g1(cls) -> (int, "G1"):
        k = rand_elem()
        return k, cls.scalar_base_mult(k)

    @classmethod
    def hash_to_point(cls, msg: bytes) -> "G1":
        # From "Indifferentiable Hashing to Barreto-Naehrig Curves"
        # https://www.di.ens.fr/~fouque/pub/latincrypt12.pdf

        # constants
        sqrt_neg_3 = sqrt_mod_p(P - 3)
        inv_2 = inv_mod_p(2)
        b = CURVE_B.value()

        # compute t in F_q
        sha = hashlib.sha512()
        sha.update(msg)
        t = int(sha.hexdigest(), 16) % P

        if t == 0:
            # TODO handle this case as described in paper
            assert False

        t2 = (t * t) % P

        chi_t = legendre(t)

        w = sqrt_neg_3 * t * inv_mod_p(1 + b + t2)

        def g(x):
            return (x * x * x + b) % P

        x1 = ((sqrt_neg_3 - 1) * inv_2 - t * w) % P
        g_x1 = g(x1)
        if legendre(g_x1) == 1:
            x1_sqrt = sqrt_mod_p(g_x1)
            return cls(CurvePoint(Gfp1(x1), Gfp1(chi_t * x1_sqrt)))

        x2 = (-1 - x1) % P
        g_x2 = g(x2)

        if legendre(g_x2) == 1:
            x2_sqrt = sqrt_mod_p(g_x2)
            return cls(CurvePoint(Gfp1(x2), Gfp1(chi_t * x2_sqrt)))

        x3 = 1 + inv_mod_p(w * w)
        g_x3 = g(x3)

        assert legendre(g_x3) == 1
        x3_sqrt = sqrt_mod_p(g_x3)
        p = CurvePoint(Gfp1(x3), Gfp1(chi_t * x3_sqrt))
        return cls(p)

    @classmethod
    def scalar_base_mult(cls, k: int) -> "G1":
        p = CURVE_G.scalar_mul(k)
        return cls(p)

    def scalar_mult(self, k: int) -> "G1":
        return G1(self.p.scalar_mul(k))

    def compress(self) -> (int, int):
        self.p.force_affine()
        x = self.p.x.value()
        y = self.p.y.value()
        return (x, y & 1)

    @classmethod
    def uncompress(cls, g: List[int]) -> "G1":
        x = g[0]
        y_sign = g[1]

        assert y_sign == 0 or y_sign == 1
        assert 0 <= x < P

        xxx = (x * x * x + CURVE_B.value()) % P

        y = sqrt_mod_p(xxx)

        if y_sign != y & 1:
            y = P - y

        p = CurvePoint(Gfp1(x), Gfp1(y))
        return cls(p)
