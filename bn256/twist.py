from .gfp2 import Gfp2
from .utils import bits_of, nums_to_bytes

TWIST_B = Gfp2(
    266929791119991161246907387137283842545076965332900288569378510910307636690,
    19485874751759354771024239261021720505790618469301721065564631296452457478373,
)


class TwistPoint(object):
    """
    TwistPoint implements the elliptic curve y²=x³+3/ξ over GF(p²). Points are
    kept in Jacobian form and t=z² when valid. The group G₂ is the set of
    n-torsion points of this curve over GF(p²) (where n = Order)
    """

    def __init__(self, x: Gfp2, y: Gfp2, z: Gfp2 = Gfp2.one(), t: Gfp2 = Gfp2.zero()):
        assert isinstance(x, Gfp2)
        assert isinstance(y, Gfp2)
        assert isinstance(z, Gfp2)
        assert isinstance(t, Gfp2)
        self.x: Gfp2 = x
        self.y: Gfp2 = y
        self.z: Gfp2 = z
        self.t: Gfp2 = t

    def __str__(self):
        return self.string()

    def __repr__(self):  # ✅
        return "<TwistPoint (%s,%s,%s)>" % (self.x, self.y, self.z)

    def __add__(self, other: "TwistPoint") -> "TwistPoint":
        return self.add(other)

    def __mul__(self, other: int) -> "TwistPoint":
        return self.mul_scalar(other)

    def __neg__(self) -> "TwistPoint":  # ✅
        return self.negative()

    def __copy__(self) -> "TwistPoint":
        return self.copy()

    def __bytes__(self) -> bytes:
        return self.bytes()

    def __eq__(self, other: "TwistPoint") -> bool:
        a = self.copy().make_affine()
        b = other.copy().make_affine()
        return a.x == b.x and a.y == b.y and a.z == b.z

    def __ne__(self, other: "TwistPoint") -> bool:
        a = self.copy().make_affine()
        b = other.copy().make_affine()
        return a.x != b.x or a.y != b.y or a.z != b.z

    @classmethod
    def zero(cls) -> "TwistPoint":  # ✅
        return cls(Gfp2.zero(), Gfp2.zero(), Gfp2.zero(), Gfp2.zero())

    @classmethod
    def one(cls) -> "TwistPoint":  # ✅
        return cls(Gfp2.zero(), Gfp2.one(), Gfp2.zero(), Gfp2.zero())

    def set_zero(self) -> "TwistPoint":
        self.x = Gfp2.zero()
        self.y = Gfp2.zero()
        self.z = Gfp2.zero()
        self.t = Gfp2.zero()
        return self

    def set_one(self) -> "TwistPoint":
        self.x = Gfp2.zero()
        self.y = Gfp2.one()
        self.z = Gfp2.zero()
        self.t = Gfp2.zero()
        return self

    def set_infinity(self) -> "TwistPoint":  # ✅
        self.z.set_zero()
        return self

    def is_infinity(self) -> bool:  # ✅
        return self.z.is_zero()

    def is_on_curve(self) -> bool:  # ✅
        y2 = self.y.square()
        x2 = self.x.square()
        x3 = x2 * self.x

        r = y2 - x3 - TWIST_B
        r.minimal()
        return r.is_zero()

    def add(self, other: "TwistPoint") -> "TwistPoint":  # ✅
        if self.is_infinity():
            return other

        if other.is_infinity():
            return self

        z1z1 = self.z.square()
        z2z2 = other.z.square()

        u1 = self.x * z2z2
        u2 = other.x * z1z1

        s1 = self.y * (other.z * z2z2)
        s2 = (other.y * (self.z * z1z1))

        h = u2 - u1
        r = s2 - s1
        if h.is_zero() and r.is_zero():
            return self.double()

        # Set x = (2h)²(s²-u1-u2)
        _4h2 = (h + h).square()
        _4h3 = (h * _4h2)
        _2r = r + r
        v = (u1 * _4h2)
        x = _2r.square() - _4h3 - (v + v)

        # Set y = -(2h)³(s1 + s*(x/4h²-u1))
        t = s1 * _4h3
        y = (_2r * (v - x)) - (t + t)

        #  Set z = 2(u2-u1)·z1·z2 = 2h·z1·z2
        t = self.z + other.z
        z = (t.square() - z1z1 - z2z2) * h

        return TwistPoint(x, y, z)

    def mul_scalar(self, k: int) -> "TwistPoint":
        assert isinstance(k, int)
        if int(k) == 0:
            return TwistPoint.zero()
        if int(k) == 1:
            return self.copy()
        r = TwistPoint.zero()
        for b in [0] + bits_of(k):
            r = r.double() + self if b != 0 else r.double()
        return r

    def double(self) -> "TwistPoint":  # ✅
        x2 = self.x.square()
        y2 = self.y.square()
        y4 = y2.square()

        t = self.x + y2
        _2xy2 = t.square() - x2 - y4
        _4xy2 = _2xy2 + _2xy2

        _3x2 = x2 + x2 + x2
        _9x4 = _3x2 * _3x2

        _8xy2 = _4xy2 + _4xy2
        x = _9x4 - _8xy2

        _2y4 = y4 + y4
        _4y4 = _2y4 + _2y4
        _8y4 = _4y4 + _4y4
        y = _3x2 * (_4xy2 - x) - _8y4

        # z = 2yz
        yz = self.y * self.z
        z = yz + yz

        return TwistPoint(x, y, z)

    def negative(self) -> "TwistPoint":  # ✅
        return TwistPoint(self.x, -self.y, self.z, self.t)

    def make_affine(self) -> "TwistPoint":
        if self.z.is_one():
            return self

        if self.z.is_zero():
            return self.set_one()

        z_inv = self.z.invert()
        z_inv2 = z_inv.square()
        z_inv3 = z_inv2 * z_inv

        x = self.x * z_inv2
        y = self.y * z_inv3
        z = Gfp2.one()
        return TwistPoint(x, y, z)

    def string(self) -> str:
        copy = self.copy().make_affine()
        return "(%s,%s)" % (copy.x, copy.y)

    def copy(self) -> "TwistPoint":
        return TwistPoint(self.x.copy(), self.y.copy(), self.z.copy(), self.t.copy())

    def bytes(self) -> bytes:
        copy = self.copy().make_affine()
        return nums_to_bytes(copy.x.x, copy.x.y, copy.y.x, copy.y.y)


TWIST_G = TwistPoint(
    Gfp2(
        11559732032986387107991004021392285783925812861821192530917403151452391805634,
        10857046999023057135944570762232829481370756359578518086990519993285655852781,
    ),
    Gfp2(
        4082367875863433681332203403145435568316851327593401208105741076214120093531,
        8495653923123431417604973247489272438418190587263600148770280649306958101930,
    ),
    Gfp2(0, 1),
    Gfp2(0, 1),
)
