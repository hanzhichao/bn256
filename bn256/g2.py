from .gfp2 import Gfp2
from .twist import TwistPoint, TWIST_G
from .utils import random_k, bytes_to_nums


class G2:

    def __init__(self, p: TwistPoint):
        self.p: TwistPoint = p

    def __repr__(self):
        return "<G2 p=%s>" % self.p

    def __add__(self, other: "G2") -> "G2":
        return self.add(other)

    def __copy__(self) -> "G2":
        return self.copy()

    def __bytes__(self) -> bytes:
        return self.marshal()

    def __eq__(self, other: "G2") -> bool:
        assert isinstance(other, G2)
        return self.p == other.p

    def __ne__(self, other: "G2") -> bool:
        assert isinstance(other, G2)
        return self.p != other.p

    @staticmethod
    def base() -> "G2":
        return G2(TWIST_G)

    @classmethod
    def scalar_base_mult(cls, k: int) -> "G2":  # ✅
        assert isinstance(k, int)
        p = cls.base().p.mul_scalar(k)
        return cls(p)

    @classmethod
    def random_g2(cls) -> (int, "G2"):  # ✅
        k = random_k()
        return k, cls.scalar_base_mult(k)

    def add(self, other: "G2") -> "G2":  # ✅
        assert isinstance(other, G2)
        # ?? a + b != b + a
        return G2(self.p + other.p)

    def scalar_mult(self, k: int) -> "G2":
        assert isinstance(k, int)
        p = self.p.mul_scalar(k)
        return G2(p)

    def marshal(self) -> bytes:
        return self.p.bytes()

    @classmethod
    def unmarshal(cls, data: bytes) -> "G2":
        xx, xy, yx, yy = bytes_to_nums(data, cnt=4)
        p = TwistPoint(Gfp2(xx, xy), Gfp2(yx, yy))
        return G2(p)

    def copy(self) -> "G2":
        p = self.p.copy()
        return G2(p)
