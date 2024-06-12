from typing import Union

from .constants import P
from .utils import bits_of
from .gfp1 import Gfp1, GFP1_ZERO, GFP1_ONE, CURVE_B

from .gfp2 import Gfp2, GFP2_ZERO, GFP2_ONE, XI

# origin
TWIST_B = XI.inverse().mul(Gfp2(0, CURVE_B.value()))

# chainmaker
# TWIST_B = Gfp2(
#     266929791119991161246907387137283842545076965332900288569378510910307636690,
#     19485874751759354771024239261021720505790618469301721065564631296452457478373,
# )


class Point(object):
    x: Union[Gfp1, Gfp2]
    y: Union[Gfp1, Gfp2]
    z: Union[Gfp1, Gfp2]

    def __repr__(self):
        self.force_affine()
        return "(%s, %s)" % (self.x, self.y)

    #     self.z = z
    def is_infinite(self):
        return self.z.is_zero()

    def add(self, other):
        return self.__add__(other)

    def __add__(self, other):
        if self.is_infinite():
            return other
        if other.is_infinite():
            return self

        """
        http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#addition-add-2007-bl
          Z1Z1 = self.z^2
          Z2Z2 = other.z^2
          U1 = self.x*Z2Z2
          U2 = other.x*Z1Z1
          S1 = self.y*other.z*Z2Z2
          S2 = other.y*self.z*Z1Z1
          H = U2-U1
          I = (2*H)^2
          J = H*I
          r = 2*(S2-S1)
          V = U1*I
          X3 = r^2-J-2*V
          Y3 = r*(V-X3)-2*S1*J
          Z3 = ((self.z+other.z)^2-Z1Z1-Z2Z2)*H
            """

        z1z1 = self.z.square()
        z2z2 = other.z.square()
        u1 = z2z2 * self.x
        u2 = z1z1 * other.x
        h = u2 - u1

        s1 = self.y * other.z * z2z2
        s2 = other.y * self.z * z1z1
        r = s2 - s1

        if h.is_zero() and r.is_zero():
            return self.double()  # FIXME

        r = r.double()
        i = h.square()
        i = i.double().double()
        j = h * i

        V = u1 * i

        c_x = r.square() - j - V.double()
        c_y = r * (V - c_x) - s1 * j.double()

        c_z = self.z + other.z
        c_z = c_z.square()
        c_z -= z1z1
        c_z -= z2z2
        c_z *= h

        return self.__class__(c_x, c_y, c_z)

    def double(self):
        # http://hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l
        """
        compute A = X1^2
            compute B = Y1^2
            compute C = B^2
            compute D = 2 ((X1 + B)^2 - A - C)
            compute E = 3 A
            compute F = E^2
            compute X3 = F - 2 D
            compute Y3 = E (D - X3) - 8 C
            compute Z3 = 2 Y1 Z1
        """
        A = self.x.square()
        B = self.y.square()
        C = B.square()

        t = self.x + B
        t = t.square()

        D = t - A - C
        D = D.double()

        E = A.double() + A
        F = E.square()

        C8 = C.double().double().double()

        c_x = F - D.double()
        c_y = E * (D - c_x) - C8
        c_z = (self.y * self.z).double()

        return self.__class__(c_x, c_y, c_z)

    @staticmethod
    def one_element() -> Union[Gfp1, Gfp2]:
        pass

    @staticmethod
    def zero_element() -> Union[Gfp1, Gfp2]:
        pass

    def _is_on_curve(self, curve) -> bool:
        self.force_affine()
        yy = self.y.square()
        xxx = self.x.square() * self.x
        yy -= xxx
        yy -= curve
        return yy.is_zero()

    def force_affine(self) -> None:
        if self.z.is_one():
            return

        if self.z.is_zero():
            return

        zinv = self.z.inverse()
        zinv2 = zinv * zinv
        zinv3 = zinv2 * zinv

        self.x = self.x * zinv2
        self.y = self.y * zinv3
        self.z = self.one_element()

    def scalar_mul(self, k: int):
        assert isinstance(k, int)

        if int(k) == 0:
            return self.__class__(self.one_element(), self.one_element(), self.zero_element())

        R = [self.__class__(self.zero_element(), self.zero_element(), self.zero_element()),
             self]

        for kb in bits_of(k):
            R[kb ^ 1] = R[kb].add(R[kb ^ 1])
            R[kb] = R[kb].double()
        return R[0]

    def negate(self):
        self.y = self.y.negative_of()


class CurvePoint(Point):
    x: Gfp1
    y: Gfp1
    z: Gfp1

    def __init__(self, x: Gfp1, y: Gfp1, z: Gfp1 = GFP1_ONE):
        assert isinstance(x, Gfp1) and isinstance(y, Gfp1) and isinstance(z, Gfp1)

        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        self.force_affine()
        return "<CurvePoint(%s, %s)>" % (self.x, self.y)

    @staticmethod
    def zero_element():
        return GFP1_ZERO

    @staticmethod
    def one_element():
        return GFP1_ONE

    def is_on_curve(self):
        return self._is_on_curve(CURVE_B)


class TwistPoint(Point):
    x: Gfp2
    y: Gfp2
    z: Gfp2

    def __init__(self, x: Gfp2, y: Gfp2, z: Gfp2 = GFP2_ONE):
        assert isinstance(x, Gfp2) and isinstance(y, Gfp2) and isinstance(z, Gfp2)
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        self.force_affine()
        return "<TwistPoint(%s, %s)>" % (self.x, self.y)

    @staticmethod
    def one_element():
        return GFP2_ONE

    @staticmethod
    def zero_element():
        return GFP2_ZERO

    def is_on_curve(self):
        return self._is_on_curve(TWIST_B)

    def negate(self):
        self.y = self.y.negative_of()


CURVE_G = CurvePoint(GFP1_ONE, Gfp1(P - 2))

# twistGen is the generator of group G₂.
# chainmaker
# TWIST_G = TwistPoint(
#     Gfp2(
#         11559732032986387107991004021392285783925812861821192530917403151452391805634,
#         10857046999023057135944570762232829481370756359578518086990519993285655852781,
#     ),
#     Gfp2(
#         4082367875863433681332203403145435568316851327593401208105741076214120093531,
#         8495653923123431417604973247489272438418190587263600148770280649306958101930,
#     ),
#     GFP2_ONE,
# )

# origin
TWIST_G = TwistPoint(
    Gfp2(
        21167961636542580255011770066570541300993051739349375019639421053990175267184,
        64746500191241794695844075326670126197795977525365406531717464316923369116492,
    ),
    Gfp2(
        20666913350058776956210519119118544732556678129809273996262322366050359951122,
        17778617556404439934652658462602675281523610326338642107814333856843981424549,
    ),
    GFP2_ONE,
)
