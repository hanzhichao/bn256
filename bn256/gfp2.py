from .constants import P

from .utils import bits_of, mod_inverse


class Gfp2(object):
    def __init__(self, x: int, y: int):
        self.x: int = x
        self.y: int = y

    @classmethod
    def zero(cls) -> "Gfp2":
        return cls(0, 0)

    @classmethod
    def one(cls) -> "Gfp2":
        return cls(0, 1)

    def __str__(self):  # ✅
        return self.string()

    def __repr__(self):  # ✅
        return "<Gfp2 (%d,%d)>" % (self.x % P, self.y % P)

    def __eq__(self, other: "Gfp2") -> bool:  # ✅
        assert isinstance(other, Gfp2)
        ax = self.x % P
        ay = self.y % P
        bx = other.x % P
        by = other.y % P
        return ax == bx and ay == by

    def __ne__(self, other: "Gfp2") -> bool:  # ✅
        assert isinstance(other, Gfp2)
        ax = self.x % P
        ay = self.y % P
        bx = other.x % P
        by = other.y % P
        return ax != bx or ay != by

    def __neg__(self) -> "Gfp2":  # ✅
        return self.negative()

    def __add__(self, other: "Gfp2") -> "Gfp2":  # ✅
        return self.add(other)

    def __sub__(self, other: "Gfp2") -> "Gfp2":  # ✅
        return self.sub(other)

    def __mul__(self, other: "Gfp2") -> "Gfp2":  # ✅
        return self.mul(other)

    def __pow__(self, k: int) -> "Gfp2":  # ✅
        return self.exp(k)

    def __invert__(self) -> "Gfp2":  # ✅
        return self.invert()

    def __copy__(self) -> "Gfp2":  # ✅
        return self.copy()

    def set(self, x: int, y: int) -> "Gfp2":
        assert isinstance(x, int) and isinstance(y, int)
        self.x = x
        self.y = y
        return self

    def set_zero(self) -> "Gfp2":
        self.x = 0
        self.y = 0
        return self

    def set_one(self) -> "Gfp2":
        self.x = 0
        self.y = 1
        return self

    def is_zero(self) -> bool:  # ✅
        return self.x % P == 0 and self.y % P == 0

    def is_one(self) -> bool:  # ✅
        return self.x % P == 0 and self.y % P == 1

    def add(self, other: "Gfp2") -> "Gfp2":  # ✅
        assert isinstance(other, Gfp2)
        x = self.x + other.x
        y = self.y + other.y
        return Gfp2(x, y)

    def sub(self, other: "Gfp2") -> "Gfp2":  # ✅
        assert isinstance(other, Gfp2)
        x = (self.x - other.x) % P
        y = (self.y - other.y) % P
        return Gfp2((self.x - other.x), (self.y - other.y))

    def mul(self, other: "Gfp2") -> "Gfp2":  # ✅
        assert isinstance(other, Gfp2)
        x = (self.x * other.y + other.x * self.y) % P
        y = (self.y * other.y - self.x * other.x) % P
        return Gfp2(x, y)

    def double(self) -> "Gfp2":  # ✅
        x = self.x + self.x
        y = self.y + self.y
        return Gfp2(x, y)

    def mul_scalar(self, k: int) -> "Gfp2":  # ✅
        x = self.x * k
        y = self.y * k
        return Gfp2(x, y)

    # Multiply by i+3
    def mul_xi(self) -> "Gfp2":  # ✅
        # (xi+y)(i+3) = (9x+y)i+(9y-x)
        x = (self.x << 3) + self.x + self.y
        y = (self.y << 3) + self.y - self.x
        return Gfp2(x, y)

    def square(self) -> "Gfp2":  # ✅
        x = self.x * self.y
        x = (x + x) % P
        y = (self.y - self.x) * (self.y + self.x) % P
        return Gfp2(x, y)

    def invert(self) -> "Gfp2":  # ✅
        t = self.y * self.y + self.x * self.x
        inv = mod_inverse(t, P)
        x = (-self.x * inv) % P
        y = (self.y * inv) % P
        return Gfp2(x, y)

    def exp(self, k: int) -> "Gfp2":  # ✅
        assert isinstance(k, int)
        r = Gfp2.one()
        for b in bits_of(k):
            r = r * r if b == 0 else r * r * self
        return r

    def negative(self) -> "Gfp2":  # ✅
        return Gfp2(-self.x, -self.y)

    def conjugate(self) -> "Gfp2":  # ✅
        return Gfp2(-self.x, self.y)

    def string(self) -> str:  # ✅
        return "(%d,%d)" % (self.x % P, self.y % P)

    def minimal(self) -> "Gfp2":
        if self.x < 0 or self.x >= P:
            self.x = self.x % P
        if self.y < 0 or self.y >= P:
            self.y = self.y % P
        return self

    def real(self) -> int:
        return self.x

    def imag(self) -> int:
        return self.y

    def copy(self) -> "Gfp2":
        return Gfp2(self.x, self.y)


XI_TO_P_MINUS_1_OVER_6 = Gfp2(
    16469823323077808223889137241176536799009286646108169935659301613961712198316,
    8376118865763821496583973867626364092589906065868298776909617916018768340080,
)

XI_TO_P_MINUS_1_OVER_3 = Gfp2(
    10307601595873709700152284273816112264069230130616436755625194854815875713954,
    21575463638280843010398324269430826099269044274347216827212613867836435027261,
)

XI_TO_P_MINUS_1_OVER_2 = Gfp2(
    3505843767911556378687030309984248845540243509899259641013678093033130930403,
    2821565182194536844548159561693502659359617185244120367078079554186484126554,
)

XI_TO_P_SQUARED_MINUS_1_OVER_3 = 21888242871839275220042445260109153167277707414472061641714758635765020556616

XI_TO_2P_MINUS_2_OVER_3 = Gfp2(
    19937756971775647987995932169929341994314640652964949448313374472400716661030,
    2581911344467009335267311115468803099551665605076196740867805258568234346338,
)

XI_TO_2P_SQUARED_MINUS_2_OVER_3 = 2203960485148121921418603742825762020974279258880205651966
