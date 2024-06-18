
from .gfp2 import (Gfp2, XI_TO_2P_MINUS_2_OVER_3, XI_TO_P_MINUS_1_OVER_3, XI_TO_2P_SQUARED_MINUS_2_OVER_3,
                   XI_TO_P_SQUARED_MINUS_1_OVER_3)


class Gfp6(object):

    def __init__(self, x: Gfp2, y: Gfp2, z: Gfp2):
        assert isinstance(x, Gfp2) and isinstance(y, Gfp2) and isinstance(z, Gfp2)
        self.x: Gfp2 = x
        self.y: Gfp2 = y
        self.z: Gfp2 = z

    @classmethod
    def zero(cls) -> "Gfp6":
        return cls(Gfp2.zero(), Gfp2.zero(), Gfp2.zero())

    @classmethod
    def one(cls) -> "Gfp6":
        return cls(Gfp2.zero(), Gfp2.zero(), Gfp2.one())

    def __str__(self) -> str:  # ✅
        return self.string()

    def __repr__(self):  # ✅
        return "<Gfp6 (%s,%s,%s)>" % (self.x, self.y, self.z)

    def __eq__(self, other: "Gfp6") -> bool:  # ✅
        assert isinstance(other, Gfp6)
        return self.x == other.x and self.y == other.y and self.z == other.z

    def __ne__(self, other: "Gfp6") -> bool:  # ✅
        assert isinstance(other, Gfp6)
        return self.x != other.x or self.y != other.y or self.z != other.z

    def __add__(self, other: "Gfp6") -> "Gfp6":  # ✅
        return self.add(other)

    def __sub__(self, other: "Gfp6") -> "Gfp6":  # ✅
        return self.sub(other)

    def __mul__(self, other: "Gfp6") -> "Gfp6":  # ✅
        return self.mul(other)

    def __neg__(self) -> "Gfp6":  # ✅
        return self.negative()

    def __invert__(self) -> "Gfp6":  # ✅
        return self.invert()

    def __copy__(self) -> "Gfp6":
        return Gfp6(self.x.copy(), self.y.__copy__(), self.z.copy())

    def set(self, xx: int, xy: int, yx: int, yy: int, zx: int, zy: int) -> "Gfp6":  # ✅
        assert isinstance(xx, int)
        assert isinstance(xy, int)
        assert isinstance(yx, int)
        assert isinstance(yy, int)
        assert isinstance(zx, int)
        assert isinstance(zy, int)

        self.x.set(xx, xy)
        self.y.set(yx, yy)
        self.z.set(zx, zy)
        return self

    def set_zero(self) -> "Gfp6":  # ✅
        self.x.set_zero()
        self.y.set_zero()
        self.z.set_zero()
        return self

    def set_one(self) -> "Gfp6":  # ✅
        self.x.set_zero()
        self.y.set_zero()
        self.z.set_one()
        return self

    def is_zero(self):  # ✅
        return self.x.is_zero() and self.y.is_zero() and self.z.is_zero()

    def is_one(self):  # ✅
        return self.x.is_zero() and self.y.is_zero() and self.z.is_one()

    def negative(self) -> "Gfp6":  # ✅
        return Gfp6(-self.x, -self.y, -self.z)

    def add(self, other: "Gfp6") -> "Gfp6":  # ✅
        return Gfp6(self.x.add(other.x), self.y.add(other.y), self.z.add(other.z))

    def sub(self, other: "Gfp6") -> "Gfp6":  # ✅
        return Gfp6(self.x.sub(other.x), self.y.sub(other.y), self.z.sub(other.z))

    def mul(self, other: "Gfp6") -> "Gfp6":  # ✅
        v2 = self.x * other.x
        v1 = self.y * other.y
        v0 = self.z * other.z

        x = (self.x + self.z) * (other.x + other.z) - v0 + v1 - v2
        y = (self.y + self.z) * (other.y + other.z) - v0 - v1 + v2.mul_xi()
        z = ((self.x + self.y) * (other.x + other.y) - v1 - v2).mul_xi() + v0

        return Gfp6(x, y, z)

    def mul_scalar(self, k: Gfp2) -> "Gfp6":  # ✅
        assert isinstance(k, Gfp2)

        return Gfp6(self.x.mul(k), self.y.mul(k), self.z.mul(k))

    def mul_gfp(self, k: int) -> "Gfp6":  # ✅
        assert isinstance(k, int)
        x = self.x.mul_scalar(k)
        y = self.y.mul_scalar(k)
        z = self.z.mul_scalar(k)
        return Gfp6(x, y, z)

    def mul_tau(self) -> "Gfp6":
        x = self.y
        y = self.z
        z = self.x.mul_xi()
        return Gfp6(x, y, z)

    def double(self) -> "Gfp6":  # ✅
        return Gfp6(self.x.double(), self.y.double(), self.z.double())

    def square(self) -> "Gfp6":  # ✅
        # Algorithm 16 from http://eprint.iacr.org/2010/354.pdf
        y2 = self.y.double()
        zy2 = self.z * y2
        xx = self.x.square()
        y = xx.mul_xi() + zy2

        zz = self.z.square()
        t1 = self.x + self.z - self.y
        t1 = t1.square()
        t2 = y2 * self.x
        z = t2.mul_xi() + zz

        x = zy2 - xx + t1 + t2 - zz
        return Gfp6(x, y, z)

    def invert(self) -> "Gfp6":  # ✅
        xx = self.x.square()
        yy = self.y.square()
        zz = self.z.square()

        xy = self.x * self.y
        xz = self.x * self.z
        yz = self.y * self.z

        a = zz - xy.mul_xi()
        b = xx.mul_xi() - yz
        c = yy - xz

        f = ~((c * self.y).mul_xi() + a * self.z + (b * self.x).mul_xi())

        x = c * f
        y = b * f
        z = a * f
        return Gfp6(x, y, z)

    def string(self) -> str:  # ✅
        return "(%s,%s,%s)" % (self.x, self.y, self.z)

    def minimal(self) -> "Gfp6":  # ✅
        self.x.minimal()
        self.y.minimal()
        self.z.minimal()
        return self

    def frobenius(self) -> "Gfp6":  # ✅
        x = self.x.conjugate()
        y = self.y.conjugate()
        z = self.z.conjugate()

        x = x * XI_TO_2P_MINUS_2_OVER_3
        y = y * XI_TO_P_MINUS_1_OVER_3
        return Gfp6(x, y, z)

    def frobenius_p2(self) -> "Gfp6":  # ✅
        x = self.x.mul_scalar(XI_TO_2P_SQUARED_MINUS_2_OVER_3)
        y = self.y.mul_scalar(XI_TO_P_SQUARED_MINUS_1_OVER_3)
        z = self.z
        return Gfp6(x, y, z)

    def copy(self) -> "Gfp6":
        return self.__copy__()

    def bytes(self): # TODO
        return [self.x.bytes() ,self.y.bytes() , self.z.bytes()]
