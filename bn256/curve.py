from .constants import P
from .utils import bits_of, mod_inverse, nums_to_bytes

CURVE_B = 3


class CurvePoint(object):
    def __init__(self, x: int, y: int, z: int = 1, t: int = 0):
        assert isinstance(x, int)
        assert isinstance(y, int)
        assert isinstance(z, int)
        assert isinstance(t, int)

        self.x: int = x
        self.y: int = y
        self.z: int = z
        self.t: int = t

    def __str__(self):
        return self.string()

    def __repr__(self):  # ✅
        return "<CurvePoint (%d,%d,%d)>" % (self.x, self.y, self.z)

    def __eq__(self, other: "CurvePoint") -> bool:  # ✅
        assert isinstance(other, CurvePoint)
        a = self.copy().make_affine()
        b = other.copy().make_affine()
        return a.x == b.x and a.y == b.y

    def __ne__(self, other: "CurvePoint") -> bool:  # ✅
        assert isinstance(other, CurvePoint)
        a = self.copy().make_affine()
        b = other.copy().make_affine()
        return a.x != b.x or a.y != b.y

    def __add__(self, other: "CurvePoint") -> "CurvePoint":  # ✅
        return self.add(other)

    def __neg__(self) -> "CurvePoint":  # ✅
        return self.negative()

    def __mul__(self, other: int) -> "CurvePoint":
        return self.mul_scalar(other)

    def __copy__(self) -> "CurvePoint":  # ✅
        return self.copy()

    def __bytes__(self):
        return self.bytes()

    @classmethod
    def zero(cls) -> "CurvePoint":
        return cls(0, 0, 0, 0)

    @classmethod
    def one(cls) -> "CurvePoint":
        return cls(0, 1, 0, 0)

    def set(self, x: int, y: int, z: int, t: int = 0) -> "CurvePoint":  # ✅
        self.x = x
        self.y = y
        self.z = z
        self.t = t
        return self

    def set_zero(self) -> "CurvePoint":
        self.x = 0
        self.y = 0
        self.z = 0
        self.t = 0
        return self

    def set_one(self) -> "CurvePoint":
        self.x = 0
        self.y = 1
        self.z = 0
        self.t = 0
        return self

    def set_infinity(self) -> "CurvePoint":  # ✅
        self.z = 0
        return self

    def is_infinity(self) -> bool:  # ✅
        return self.z == 0

    def is_on_curve(self) -> bool:  # ✅
        r = (self.y * self.y) - self.x * self.x * self.x - CURVE_B
        if r < 0 or r >= P:
            r = r % P
        return r == 0

    def add(self, other: "CurvePoint") -> "CurvePoint":  # ✅
        if self.is_infinity():
            return other

        if other.is_infinity():
            return self

        # z1² mod P
        z1z1 = (self.z * self.z) % P
        # z2² mod P
        z2z2 = (other.z * other.z) % P

        # u1 = x1 * z2² mod P
        u1 = (self.x * z2z2) % P
        # u1 = x2 * z1² mod P
        u2 = (other.x * z1z1) % P
        # s1 = y1 * z2³ mod P
        s1 = (self.y * (other.z * z2z2 % P)) % P
        # s2 = y2 * z1³ mod P
        s2 = (other.y * (self.z * z1z1 % P)) % P

        # h = x2 * z1² - x1 * z2²
        h = u2 - u1
        # r = y2 * z1³ - y1 * z2³
        r = s2 - s1
        if h == 0 and r == 0:
            return self.double()

        # Set x = (2h)²(s²-u1-u2)
        _4h2 = (h + h) ** 2 % P
        _4h3 = (h * _4h2) % P
        _2r = r + r
        v = (u1 * _4h2) % P
        x = (_2r * _2r) % P - _4h3 - (v + v)

        # Set y = -(2h)³(s1 + s*(x/4h²-u1))
        t = (s1 * _4h3) % P
        y = (_2r * (v - x)) % P - (t + t)

        #  Set z = 2(u2-u1)·z1·z2 = 2h·z1·z2
        t = self.z + other.z
        z = ((t * t) % P - z1z1 - z2z2) * h % P

        return CurvePoint(x, y, z)

    def double(self) -> "CurvePoint":
        # A = x²
        x2 = (self.x * self.x) % P
        # B = y²
        y2 = (self.y * self.y) % P
        # C = y⁴
        y4 = (y2 * y2) % P

        # t = x + y²
        t = self.x + y2
        # (x + y²)² - x² - y⁴ = 2xy²
        _2xy2 = (t * t) % P - x2 - y4
        _4xy2 = _2xy2 + _2xy2

        _3x2 = x2 + x2 + x2
        _9x4 = (_3x2 * _3x2) % P

        _8xy2 = _4xy2 + _4xy2
        x = _9x4 - _8xy2

        _2y4 = y4 + y4
        _4y4 = _2y4 + _2y4
        _8y4 = _4y4 + _4y4
        y = (_3x2 * (_4xy2 - x)) % P - _8y4

        # z = 2yz
        yz = (self.y * self.z) % P
        z = yz + yz

        return CurvePoint(x, y, z)

    def mul_scalar(self, k: int) -> "CurvePoint":  # ✅
        assert isinstance(k, int)
        if int(k) == 0:
            return CurvePoint.zero()
        if int(k) == 1:
            return self.copy()
        r = CurvePoint.zero()
        for b in [0] + bits_of(k):
            r = r.double() + self if b != 0 else r.double()
        return r

    def negative(self) -> "CurvePoint":  # ✅
        return CurvePoint(self.x, -self.y, self.z)

    def make_affine(self) -> "CurvePoint":  # ✅
        if self.z == 1:
            return self

        if self.is_infinity():  # self.z == 0
            return self.set_one()

        z_inv = mod_inverse(self.z, P)
        z_inv2 = (z_inv * z_inv) % P

        self.x = (self.x * z_inv2) % P
        self.y = ((self.y * z_inv) % P * z_inv2) % P
        self.z = 1
        self.t = 1
        return self

    def string(self) -> str:  # ✅
        copy = self.copy().make_affine()
        return "(%d,%d)" % (copy.x, copy.y)

    def copy(self) -> "CurvePoint":  # ✅
        return CurvePoint(self.x, self.y, self.z, self.t)

    def bytes(self) -> bytes:
        copy = self.copy()
        copy.make_affine()
        return nums_to_bytes(copy.x, copy.y)


CURVE_G = CurvePoint(1, 2, 1, 1)
