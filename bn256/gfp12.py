from .gfp6 import Gfp6
from .gfp2 import Gfp2
from .utils import bits_of

XI_TO_P_MINUS_1_OVER_6 = Gfp2(
    16469823323077808223889137241176536799009286646108169935659301613961712198316,
    8376118865763821496583973867626364092589906065868298776909617916018768340080,
)

XI_TO_P_SQUARED_MINUS_1_OVER_6 = 21888242871839275220042445260109153167277707414472061641714758635765020556617


class Gfp12(object):
    def __init__(self, x: Gfp6, y: Gfp6 = None):
        assert isinstance(x, Gfp6) and isinstance(y, Gfp6)
        self.x: Gfp6 = x
        self.y: Gfp6 = y

    def __eq__(self, other: "Gfp12") -> bool:  # ✅
        return self.x == other.x and self.y == other.y

    def __ne__(self, other: "Gfp12") -> bool:  # ✅
        return self.x != other.x or self.y == other.y

    def __str__(self):  # ✅
        return self.string()

    def __repr__(self):  # ✅
        return "<Gfp12 %s>" % self

    def __add__(self, other: "Gfp12") -> "Gfp12":  # ✅
        return self.add(other)

    def __sub__(self, other: "Gfp12") -> "Gfp12":  # ✅
        return self.sub(other)

    def __mul__(self, other: "Gfp12") -> "Gfp12":  # ✅
        return self.mul(other)

    def __neg__(self) -> "Gfp12":  # ✅
        return self.negative()

    def __pow__(self, k: int) -> "Gfp12":  # ✅
        return self.exp(k)

    def __copy__(self) -> "Gfp12":  # ✅
        return self.copy()

    @classmethod
    def zero(cls) -> "Gfp12":  # ✅
        return cls(Gfp6.zero(), Gfp6.zero())

    @classmethod
    def one(cls) -> "Gfp12":  # ✅
        return cls(Gfp6.zero(), Gfp6.one())

    def is_zero(self) -> bool:  # ✅
        return self.x.is_zero() and self.y.is_zero()

    def is_one(self) -> bool:  # ✅
        return self.x.is_zero() and self.y.is_one()

    def conjugate(self) -> "Gfp12":  # ✅
        return Gfp12(-self.x, self.y)

    def negative(self) -> "Gfp12":  # ✅
        return Gfp12(-self.x, -self.y)

    def frobenius(self) -> "Gfp12":  # ✅
        x = self.x.frobenius().mul_scalar(XI_TO_P_MINUS_1_OVER_6)
        y = self.y.frobenius()
        return Gfp12(x, y)

    def frobenius_p2(self) -> "Gfp12":  # ✅
        x = self.x.frobenius_p2().mul_gfp(XI_TO_P_SQUARED_MINUS_1_OVER_6)
        y = self.y.frobenius_p2()
        return Gfp12(x, y)

    def add(self, other: "Gfp12") -> "Gfp12":  # ✅
        assert isinstance(other, Gfp12)
        x = self.x + other.x
        y = self.y + other.y
        return Gfp12(x, y)

    def sub(self, other: "Gfp12") -> "Gfp12":  # ✅
        assert isinstance(other, Gfp12)
        x = self.x - other.x
        y = self.y - other.y
        return Gfp12(x, y)

    def mul(self, other: "Gfp12") -> "Gfp12":  #
        assert isinstance(other, Gfp12)
        x = (self.x * other.y) + (other.x * self.y)
        y = (self.y * other.y) + (self.x * other.x).mul_tau()
        return Gfp12(x, y)

    def mul_scalar(self, k: Gfp6) -> "Gfp12":  # ✅
        assert isinstance(k, Gfp6)
        x = self.x.mul(k)
        y = self.y.mul(k)
        return Gfp12(x, y)

    def exp(self, k: int) -> "Gfp12":  #  ✅
        assert isinstance(k, int)
        r = Gfp12.one()
        for b in bits_of(k):
            r = r.square() * self if b != 0 else r.square()
        return r

    def square(self) -> "Gfp12":  # ✅
        t = self.x * self.y
        x = t.double()
        y = (self.x + self.y) * (self.y + self.x.mul_tau()) - t - t.mul_tau()

        return Gfp12(x, y)

    def __invert__(self) -> "Gfp12":  # ✅
        t = ~(self.y.square() - self.x.square().mul_tau())
        inv = Gfp12(-self.x, self.y).mul_scalar(t)
        return inv

    def invert(self) -> "Gfp12":  # ✅
        return self.__invert__()

    def set_zero(self) -> "Gfp12":  # ✅
        self.x.set_zero()
        self.y.set_zero()
        return self

    def set_one(self) -> "Gfp12":  # ✅
        self.x.set_zero()
        self.y.set_one()
        return self

    def set(self, axx: int, axy: int, ayx: int, ayy: int, azx: int, azy: int,
            bxx: int, bxy: int, byx: int, byy: int, bzx: int, bzy: int) -> "Gfp12":  # ✅
        self.x.set(axx, axy, ayx, ayy, azx, azy)
        self.y.set(bxx, bxy, byx, byy, bzx, bzy)
        return self

    def string(self) -> str:  # ✅
        return "(%s,%s)" % (self.x, self.y)

    def copy(self) -> "Gfp12":  # ✅
        return Gfp12(self.x.copy(), self.y.copy())

    def minimal(self) -> "Gfp12":  # ✅
        self.x.minimal()
        self.y.minimal()
        return self
