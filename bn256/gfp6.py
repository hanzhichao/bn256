from .gfp2 import Gfp2, GFP2_ZERO, GFP2_ONE


class Gfp6(object):
    x: Gfp2
    y: Gfp2
    z: Gfp2

    def __init__(self, x: Gfp2, y: Gfp2, z: Gfp2):
        assert isinstance(x, Gfp2) and isinstance(y, Gfp2) and isinstance(z, Gfp2)
        self.x = x
        self.y = y
        self.z = z

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y and self.z == other.z

    def __repr__(self):
        return "(%s,%s,%s)" % (self.x, self.y, self.z)

    def is_zero(self):
        return self.x.is_zero() and self.y.is_zero() and self.z.is_zero()

    def is_one(self):
        return self.x.is_zero() and self.y.is_zero() and self.z.is_one()

    def negative_of(self):
        return Gfp6(self.x.negative_of(), self.y.negative_of(), self.z.negative_of())

    def add(self, other: "Gfp6"):
        return Gfp6(self.x.add(other.x), self.y.add(other.y), self.z.add(other.z))

    def sub(self, other: "Gfp6"):
        return Gfp6(self.x.sub(other.x), self.y.sub(other.y), self.z.sub(other.z))

    def double(self):
        return Gfp6(self.x.double(), self.y.double(), self.z.double())

    def mul(self, other: "Gfp6") -> "Gfp6":
        # Algorithm 13 from http://eprint.iacr.org/2010/354.pdf
        # plus some short-circuits

        if self.x.is_zero():
            if self.y.is_zero():
                return other.mul_scalar(self.z)

            t0 = other.z * self.z
            t1 = other.y * self.y

            tz = (other.x + other.y) * (self.y)
            tz -= t1
            tz = tz.mul_xi()
            tz += t0

            ty = (other.y + other.z) * (self.y + self.z)
            ty -= t0
            ty -= t1

            tx = (other.x) * (self.z)
            tx += t1

            return Gfp6(tx, ty, tz)

        if other.x.is_zero():
            if other.y.is_zero():
                return self.mul_scalar(other.z)

            t0 = self.z * other.z
            t1 = self.y * other.y

            tz = (self.x + self.y) * (other.y)
            tz -= t1
            tz = tz.mul_xi()
            tz += t0

            ty = (self.y + self.z) * (other.y + other.z)
            ty -= t0
            ty -= t1

            tx = (self.x) * (other.z)
            tx += t1

            return Gfp6(tx, ty, tz)

        t0 = self.z * other.z
        t1 = self.y * other.y
        t2 = self.x * other.x

        tz = (self.x + self.y) * (other.x + other.y)
        tz -= t1
        tz -= t2
        tz = tz.mul_xi()
        tz += t0

        ty = (self.y + self.z) * (other.y + other.z)
        ty -= t0
        ty -= t1
        ty += t2.mul_xi()

        tx = (self.x + self.z) * (other.x + other.z)
        tx -= t0
        tx += t1
        tx -= t2

        return Gfp6(tx, ty, tz)

    def __mul__(self, other: "Gfp6") -> "Gfp6":
        return self.mul(other)

    def __add__(self, other: "Gfp6") -> "Gfp6":
        return self.add(other)

    def __sub__(self, other: "Gfp6") -> "Gfp6":
        return self.sub(other)

    def mul_scalar(self, k: Gfp2) -> "Gfp6":
        assert isinstance(k, Gfp2)

        return Gfp6(self.x.mul(k), self.y.mul(k), self.z.mul(k))

    def mul_tau(self) -> "Gfp6":
        tx = self.y
        ty = self.z
        tz = self.x.mul_xi()
        return Gfp6(tx, ty, tz)

    def square(self) -> "Gfp6":
        # Algorithm 16 from http://eprint.iacr.org/2010/354.pdf
        ay2 = self.y.double()
        c4 = self.z * ay2
        c5 = self.x.square()
        c1 = c5.mul_xi() + c4
        c2 = c4 - c5
        c3 = self.z.square()
        c4 = self.x + self.z - self.y
        c5 = ay2 * self.x
        c4 = c4.square()
        c0 = c5.mul_xi() + c3
        c2 = c2 + c4 + c5 - c3
        n = Gfp6(c2, c1, c0)
        return n

    def inverse(self) -> "Gfp6":
        # Algorithm 17
        XX = self.x.square()
        YY = self.y.square()
        ZZ = self.z.square()

        XY = self.x * self.y
        XZ = self.x * self.z
        YZ = self.y * self.z

        A = ZZ - XY.mul_xi()
        B = XX.mul_xi() - YZ
        # There is an error in the paper for this line
        C = YY - XZ

        F = (C * self.y).mul_xi()
        F += A * self.z
        F += (B * self.x).mul_xi()

        F = F.inverse()

        c_x = C * F
        c_y = B * F
        c_z = A * F
        return Gfp6(c_x, c_y, c_z)


GFP6_ZERO = Gfp6(GFP2_ZERO, GFP2_ZERO, GFP2_ZERO)
GFP6_ONE = Gfp6(GFP2_ZERO, GFP2_ZERO, GFP2_ONE)
