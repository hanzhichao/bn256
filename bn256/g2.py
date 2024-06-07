from .gfp1 import Gfp1
from .gfp2 import Gfp2, GFP2_ZERO, GFP2_ONE
from .point import TwistPoint, TWIST_G
from .utils import rand_elem, split_by_length


class G2:
    def __init__(self, p: TwistPoint):
        self.p = p

    def __repr__(self):
        return "<G2 p=%s>" % self.p

    @classmethod
    def scalar_base_mult(cls, k: int) -> "G2":
        return cls(TWIST_G.scalar_mul(k))

    @classmethod
    def random_g2(cls) -> (int, "G2"):
        k = rand_elem()
        return k, cls.scalar_base_mult(k)

    def add(self, other: "G2") -> "G2":
        return G2(self.p.add(other.p))

    def scalar_mult(self, k: int) -> "G2":
        return G2(self.p.scalar_mul(k))

    def marshal(self) -> bytes:
        self.p.force_affine()
        return b"".join(
            map(
                lambda x: x.to_bytes(32, byteorder="big"),
                [self.p.x.x.value(), self.p.x.y.value(), self.p.y.x.value(), self.p.y.y.value()],
            ),
        )

    @classmethod
    def unmarshal(cls, data: bytes) -> "G2":
        if len(data) == 0:
            return G2(TwistPoint(GFP2_ZERO, GFP2_ONE, GFP2_ZERO))

        num_bytes = 32
        x_x = int.from_bytes(data[0:num_bytes], "big")
        x_y = int.from_bytes(data[num_bytes : num_bytes * 2], "big")
        y_x = int.from_bytes(data[num_bytes * 2 : num_bytes * 3], "big")
        y_y = int.from_bytes(data[num_bytes * 3 : num_bytes * 4], "big")
        return G2(TwistPoint(Gfp2(x_x, x_y), Gfp2(y_x, y_y)))
