from .curve import CurvePoint, CURVE_G
from .utils import random_k, bytes_to_nums


class G1:
    def __init__(self, p: CurvePoint):
        self.p: CurvePoint = p

    def __repr__(self):
        return "<G1 p=%s>" % self.p

    def __add__(self, other: "G1") -> "G1":
        return self.add(other)

    def __neg__(self) -> "G1":
        return self.neg()

    def __mul__(self, other: int) -> "G1":
        return self.scalar_mult(other)

    def __copy__(self) -> "G1":
        return self.copy()

    def __bytes__(self) -> bytes:
        return self.marshal()

    def __eq__(self, other: "G1") -> bool:
        assert isinstance(other, G1)
        return self.p == other.p

    def __ne__(self, other: "G1") -> bool:
        assert isinstance(other, G1)
        return self.p != other.p

    @staticmethod
    def base() -> "G1":  # ✅
        return G1(CURVE_G)

    @classmethod
    def random_g1(cls) -> (int, "G1"):  # ✅
        """RandomG1 返回 x 和 g₁x，其中 x 是从 r 读取的随机非零数。"""
        k = random_k()
        return k, cls.scalar_base_mult(k)

    @classmethod
    def unmarshal(cls, data: bytes) -> "G1":  # ✅
        x, y = bytes_to_nums(data, cnt=2)
        p = CurvePoint(x, y)
        return G1(p)

    @classmethod
    def scalar_base_mult(cls, k: int) -> "G1":  # ✅
        p = cls.base().p.mul_scalar(k)
        return cls(p)

    def add(self, other: "G1") -> "G1":  # ✅
        assert isinstance(other, G1)
        p = self.p + other.p
        return G1(p)

    def marshal(self) -> bytes:  # ✅
        return self.p.bytes()

    def scalar_mult(self, k: int) -> "G1":  # ✅
        assert isinstance(k, int)
        p = self.p.mul_scalar(k)
        return G1(p)

    def neg(self) -> "G1":  # ✅
        p = self.p.negative()
        return G1(p)

    def copy(self) -> "G1":  # ✅
        p = self.p.copy()
        return G1(p)
