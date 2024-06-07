from .gfp6 import Gfp6, GFP6_ZERO, GFP6_ONE
from .gfp2 import XI1, XI2
from .utils import bits_of


class Gfp12(object):
    x: Gfp6
    y: Gfp6

    def __init__(self, x: Gfp6, y: Gfp6 = None):
        assert isinstance(x, Gfp6) and isinstance(y, Gfp6)
        self.x = x
        self.y = y

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def __repr__(self):
        return "(%s,%s)" % (self.x, self.y)

    def is_zero(self) -> bool:
        return self.x.is_zero() and self.y.is_zero()

    def is_one(self) -> bool:
        return self.x.is_zero() and self.y.is_one()

    def conjugate_of(self) -> "Gfp12":
        return Gfp12(self.x.negative_of(), self.y)

    def negative_of(self) -> "Gfp12":
        return Gfp12(self.x.negative_of(), self.y.negative_of())

    def frobenius(self) -> "Gfp12":
        e1_x = self.x.x.conjugate_of().mul(XI1[4])
        e1_y = self.x.y.conjugate_of().mul(XI1[2])
        e1_z = self.x.z.conjugate_of().mul(XI1[0])

        e2_x = self.y.x.conjugate_of().mul(XI1[3])
        e2_y = self.y.y.conjugate_of().mul(XI1[1])
        e2_z = self.y.z.conjugate_of()

        return Gfp12(Gfp6(e1_x, e1_y, e1_z), Gfp6(e2_x, e2_y, e2_z))

    def frobenius_p2(self) -> "Gfp12":
        e1_x = self.x.x.mul(XI2[4])
        e1_y = self.x.y.mul(XI2[2])
        e1_z = self.x.z.mul(XI2[0])

        e2_x = self.y.x.mul(XI2[3])
        e2_y = self.y.y.mul(XI2[1])
        e2_z = self.y.z

        return Gfp12(Gfp6(e1_x, e1_y, e1_z), Gfp6(e2_x, e2_y, e2_z))

    def sub(self, other: "Gfp12") -> "Gfp12":
        return Gfp12(self.x - other.x, self.y - other.y)

    def mul(self, other: "Gfp12") -> "Gfp12":
        # TODO Karatsuba (algo 20)
        AXBX = self.x * other.x
        AXBY = self.x * other.y
        AYBX = self.y * other.x
        AYBY = self.y * other.y
        return Gfp12(AXBY + AYBX, AYBY + AXBX.mul_tau())

    def mul_scalar(self, k: Gfp6):
        assert isinstance(k, Gfp6)
        return Gfp12(self.x.mul(k), self.y.mul(k))

    def exp(self, k: int) -> "Gfp12":
        assert isinstance(k, int)

        R = [Gfp12(GFP6_ZERO, GFP6_ONE), self]

        for kb in bits_of(k):
            R[kb ^ 1] = R[kb].mul(R[kb ^ 1])
            R[kb] = R[kb].square()

        return R[0]

    def square(self) -> "Gfp12":
        v0 = self.x * self.y
        t = self.x.mul_tau()
        t += self.y
        ty = self.x + self.y
        ty *= t
        ty -= v0
        t = v0.mul_tau()
        ty -= t

        c_x = v0.double()
        c_y = ty

        return Gfp12(c_x, c_y)

    def inverse(self) -> "Gfp12":
        e = Gfp12(self.x.negative_of(), self.y)

        t1 = self.x.square()
        t2 = self.y.square()
        t1 = t1.mul_tau()
        t1 = t2 - t1
        t2 = t1.inverse()

        e = e.mul_scalar(t2)
        return e
