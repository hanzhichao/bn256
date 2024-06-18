from .g1 import G1
from .g2 import G2
from .gfp12 import Gfp12, Gfp6, Gfp2
from .optate import optimal_ate
from .utils import nums_to_bytes, bytes_to_nums


class GT(object):

    def __init__(self, p: Gfp12):
        self.p: Gfp12 = p

    def __repr__(self):
        return '<Gt p=%s>' % self.p

    def __add__(self, other: "GT") -> "GT":  # Fixme ❌
        return self.add(other)

    def __neg__(self) -> "GT":
        return self.neg()

    def __eq__(self, other: "GT") -> bool:
        assert isinstance(other, GT)
        return self.p == other.p

    def __ne__(self, other: "GT") -> bool:
        assert isinstance(other, GT)
        return self.p != other.p

    @classmethod
    def base(cls) -> "GT":
        return GT.pair(G1.base(), G2.base())

    @classmethod
    def random_gt(cls) -> (int, int, "GT"):
        k1, g1 = G1.random_g1()
        k2, g2 = G2.random_g2()
        gt = GT.pair(g1, g2)
        return k1, k2, gt

    def scalar_mult(self, k: int) -> "GT":  # ✅
        assert isinstance(k, int)
        return GT(self.p.exp(k))

    def add(self, other: "GT") -> "GT":
        assert isinstance(other, GT)
        return GT(self.p.mul(other.p))

    def neg(self) -> "GT":
        return GT(-self.p)

    def marshal(self) -> bytes:  # ✅
        self.p.minimal()
        return nums_to_bytes(self.p.x.x.x, self.p.x.x.y, self.p.x.y.x, self.p.x.y.y, self.p.x.z.x, self.p.x.z.y,
                             self.p.y.x.x, self.p.y.x.y, self.p.y.y.x, self.p.y.y.y, self.p.y.z.x, self.p.y.z.y)

    @classmethod
    def pair(cls, g1: G1, g2: G2) -> "GT":
        assert isinstance(g1, G1)
        assert isinstance(g2, G2)
        p = optimal_ate(g2.p, g1.p)
        return cls(p)

    @classmethod
    def unmarshal(cls, data: bytes) -> "GT":
        p0_x, p0_y, p1_x, p1_y, p2_x, p2_y, p3_x, p3_y, p4_x, p4_y, p5_x, p5_y = bytes_to_nums(data, 12)
        p = Gfp12(
            Gfp6(Gfp2(p0_x, p0_y), Gfp2(p1_x, p1_y), Gfp2(p2_x, p2_y)),
            Gfp6(Gfp2(p3_x, p3_y), Gfp2(p4_x, p4_y), Gfp2(p5_x, p5_y)),
        )
        return cls(p)
