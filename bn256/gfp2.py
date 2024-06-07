from typing import Union

from .gfp1 import P, Gfp1, CURVE_B
from .utils import bits_of


class Gfp2(object):
    x: Gfp1
    y: Gfp1

    def __init__(self, x: Union[Gfp1, int], y: Union[Gfp1, int]):
        """
        Represented as i*x + y
        """
        if isinstance(x, Gfp1):
            self.x = x
            self.y = y
        else:
            # Assumed to be integers
            self.x = Gfp1(x)
            self.y = Gfp1(y)

    def __str__(self):
        return "(%d,%d)" % (self.x.value(), self.y.value())

    def __repr__(self):
        return "<Gfp2 (%d,%d)>" % (self.x.value(), self.y.value())

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def __ne__(self, other):
        return self.x != other.x or self.y != other.y

    def is_zero(self) -> bool:
        return self.x.is_zero() and self.y.is_zero()

    def is_one(self) -> bool:
        return self.x.is_zero() and self.y.is_one()

    def conjugate_of(self) -> "Gfp2":
        """
        For gamma = A + iB \in gfp2
        gamma^p = A - iB
        """
        return Gfp2(self.x.additive_inverse(), self.y)

    def negative_of(self) -> "Gfp2":
        return Gfp2(self.x.additive_inverse(), self.y.additive_inverse())

    def add(self, other: "Gfp2") -> "Gfp2":
        assert isinstance(other, Gfp2)
        return Gfp2((self.x + other.x), (self.y + other.y))

    def sub(self, other: "Gfp2") -> "Gfp2":
        assert isinstance(other, Gfp2)
        return Gfp2((self.x - other.x), (self.y - other.y))

    def double(self) -> "Gfp2":
        return Gfp2((self.x.double()), (self.y.double()))

    def mul(self, other: "Gfp2") -> "Gfp2":
        assert isinstance(other, Gfp2)
        # Karatsuba
        vy = self.y * other.y
        vx = self.x * other.x
        c0 = vy - vx
        c1 = (self.x + self.y) * (other.x + other.y) - vy - vx

        return Gfp2(c1, c0)

    def __mul__(self, other: "Gfp2") -> "Gfp2":
        return self.mul(other)

    def __sub__(self, other: "Gfp2") -> "Gfp2":
        return self.sub(other)

    def __add__(self, other: "Gfp2") -> "Gfp2":
        return self.add(other)

    def mul_scalar(self, k) -> "Gfp2":
        return Gfp2((self.x * k), (self.y * k))

    # Multiply by i+3
    def mul_xi(self) -> "Gfp2":
        # (xi + y)(3 + i) = 3xi + 3y - x + yi = (3x + y)i + (3y - x)
        tx = (self.x.triple()) + self.y
        ty = (self.y.triple()) - self.x
        return Gfp2(tx, ty)

    def square(self) -> "Gfp2":
        assert isinstance(self.x, Gfp1) and isinstance(self.y, Gfp1)
        # Complex squaring
        t1 = self.y - self.x
        t2 = self.y + self.x
        ty = t1 * t2
        # ty = self.y*self.y - self.x*self.x
        tx = self.x * self.y
        tx = tx.double()
        return Gfp2(tx, ty)

    def inverse(self) -> "Gfp2":
        # Algorithm 8 from http://eprint.iacr.org/2010/354.pdf
        t = self.x.square() + self.y.square()
        inv = t.inverse()
        c_x = self.x.additive_inverse() * inv
        c_y = self.y * inv

        return Gfp2(c_x, c_y)

    def exp(self, k: int) -> "Gfp2":
        assert isinstance(k, int) and isinstance(self, Gfp2)
        R = [Gfp2(Gfp1(0), Gfp1(1)), self]
        for kb in bits_of(k):
            R[kb ^ 1] = R[kb].mul(R[kb ^ 1])
            assert type(R[kb]) == Gfp2
            R[kb] = R[kb].square()  # FIXME
        return R[0]


GFP2_ZERO = Gfp2(0, 0)
GFP2_ONE = Gfp2(0, 1)


XI = Gfp2(1, 3)  # i + 3

XI1 = [
    XI.exp(1 * (P - 1) // 6),
    XI.exp(2 * (P - 1) // 6),
    XI.exp(3 * (P - 1) // 6),
    XI.exp(4 * (P - 1) // 6),
    XI.exp(5 * (P - 1) // 6),
]

XI2 = [(x * x.conjugate_of()) for x in XI1]
