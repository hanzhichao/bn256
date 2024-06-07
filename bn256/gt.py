import copy
import hashlib

from .constants import U
from .g1 import G1
from .g2 import G2
from .gfp1 import Gfp1
from .gfp12 import Gfp12
from .gfp2 import Gfp2, GFP2_ZERO, GFP2_ONE, XI1, XI2
from .gfp6 import Gfp6, GFP6_ZERO, GFP6_ONE
from .point import TwistPoint, CurvePoint


def _to_naf(x):
    z = []
    while x > 0:
        if x % 2 == 0:
            z.append(0)
        else:
            zi = 2 - (x % 4)
            x -= zi
            z.append(zi)
        x = x // 2
    return z


# 6u+2 in NAF
NAF_6UP2 = list(reversed(_to_naf(6 * U + 2)))[1:]


def line_func_add(r: TwistPoint, p: TwistPoint, q: CurvePoint, r2: Gfp2):
    assert (
        isinstance(r, TwistPoint) and isinstance(p, TwistPoint) and isinstance(q, CurvePoint) and isinstance(r2, Gfp2)
    )

    r_t = r.z.square()
    B = p.x * r_t
    D = p.y + r.z
    D = D.square()
    D -= r2
    D -= r_t
    D *= r_t

    H = B - r.x
    I = H.square()

    E = I.double().double()

    J = H * E
    L1 = D - r.y
    L1 -= r.y

    V = r.x * E

    r_x = L1.square()
    r_x -= J
    r_x -= V.double()

    r_z = r.z + H
    r_z = r_z.square()
    r_z -= r_t
    r_z -= I

    t = V - r_x
    t *= L1
    t2 = r.y * J
    t2 = t2.double()
    r_y = t - t2

    r_out = TwistPoint(r_x, r_y, r_z)

    t = p.y + r_z
    t = t.square()
    t = t - r2
    t = t - (r_z.square())

    t2 = L1 * p.x
    t2 = t2.double()
    a = t2 - t

    c = r_z.mul_scalar(q.y).double()

    b = L1.negative_of()
    b = b.mul_scalar(q.x).double()

    return (a, b, c, r_out)


def line_func_double(r: TwistPoint, q: CurvePoint):
    assert isinstance(r, TwistPoint) and isinstance(q, CurvePoint)

    # cache this?
    r_t = r.z.square()

    A = r.x.square()
    B = r.y.square()
    C = B.square()

    D = r.x + B
    D = D.square()
    D -= A
    D -= C
    D = D.double()

    E = A.double() + A
    F = E.square()

    C8 = C.double().double().double()  # C*8

    r_x = F - D.double()
    r_y = E * (D - r_x) - C8

    # (y+z)*(y+z) - (y*y) - (z*z) = 2*y*z
    r_z = (r.y + r.z).square() - B - r_t

    assert r_z == r.y * r.z.double()

    r_out = TwistPoint(r_x, r_y, r_z)
    # assert r_out.is_on_curve()

    a = r.x + E
    a = a.square()
    a -= A + F + B.double().double()

    t = E * r_t
    t = t.double()
    b = t.negative_of()
    b = b.mul_scalar(q.x)

    c = r_z * r_t
    c = c.double().mul_scalar(q.y)

    return (a, b, c, r_out)


def mul_line(r: Gfp12, a: Gfp2, b: Gfp2, c: Gfp2):
    assert isinstance(r, Gfp12) and isinstance(a, Gfp2) and isinstance(b, Gfp2) and isinstance(c, Gfp2)

    # See function fp12e_mul_line in dclxvi

    t1 = Gfp6(GFP2_ZERO, a, b)
    t2 = Gfp6(GFP2_ZERO, a, b + c)

    t1 = t1 * r.x
    t3 = r.y.mul_scalar(c)
    r.x += r.y
    r.y = t3
    r.x *= t2
    r.x -= t1
    r.x -= r.y
    r.y += t1.mul_tau()


def miller(q: TwistPoint, p: CurvePoint):
    assert isinstance(q, TwistPoint) and isinstance(p, CurvePoint)

    Q = copy.deepcopy(q)
    Q.force_affine()

    P = copy.deepcopy(p)
    P.force_affine()

    mQ = copy.deepcopy(Q)
    mQ.negate()

    f = Gfp12(GFP6_ZERO, GFP6_ONE)
    T = Q

    Qp = Q.y.square()

    for naf_i in NAF_6UP2:
        # Skip on first iteration?
        f = f.square()

        a, b, c, T = line_func_double(T, P)
        mul_line(f, a, b, c)

        if naf_i == 1:
            a, b, c, T = line_func_add(T, Q, P, Qp)
            mul_line(f, a, b, c)
        elif naf_i == -1:
            a, b, c, T = line_func_add(T, mQ, P, Qp)
            mul_line(f, a, b, c)

    # Q1 = pi(Q)
    Q1 = TwistPoint(Q.x.conjugate_of().mul(XI1[1]), Q.y.conjugate_of().mul(XI1[2]), GFP2_ONE)

    # Q2 = pi2(Q)
    Q2 = TwistPoint(Q.x.mul_scalar(XI2[1].y), Q.y, GFP2_ONE)

    Qp = Q1.y.square()
    a, b, c, T = line_func_add(T, Q1, P, Qp)
    mul_line(f, a, b, c)

    Qp = Q2.y.square()
    a, b, c, T = line_func_add(T, Q2, P, Qp)
    mul_line(f, a, b, c)

    return f


def final_exp(inp: Gfp12, u=U):
    assert isinstance(inp, Gfp12)

    # Algorithm 31 from https://eprint.iacr.org/2010/354.pdf

    t1 = inp.conjugate_of()
    inv = inp.inverse()

    t1 = t1.mul(inv)
    # Now t1 = inp^(p**6-1)

    t2 = t1.frobenius_p2()
    t1 = t1.mul(t2)

    fp1 = t1.frobenius()
    fp2 = t1.frobenius_p2()
    fp3 = fp2.frobenius()

    fu1 = t1.exp(u)
    fu2 = fu1.exp(u)
    fu3 = fu2.exp(u)

    y3 = fu1.frobenius()
    fu2p = fu2.frobenius()
    fu3p = fu3.frobenius()
    y2 = fu2.frobenius_p2()

    y0 = fp1.mul(fp2)
    y0 = y0.mul(fp3)

    y1 = t1.conjugate_of()
    y5 = fu2.conjugate_of()
    y3 = y3.conjugate_of()
    y4 = fu1.mul(fu2p)
    y4 = y4.conjugate_of()

    y6 = fu3.mul(fu3p)
    y6 = y6.conjugate_of()

    t0 = y6.square()
    t0 = t0.mul(y4)
    t0 = t0.mul(y5)

    t1 = y3.mul(y5)
    t1 = t1.mul(t0)
    t0 = t0.mul(y2)
    t1 = t1.square()
    t1 = t1.mul(t0)
    t1 = t1.square()
    t0 = t1.mul(y1)
    t1 = t1.mul(y0)
    t0 = t0.square()
    t0 = t0.mul(t1)

    return t0


def optimal_ate(a: TwistPoint, b: CurvePoint) -> Gfp12:
    assert isinstance(a, TwistPoint) and isinstance(b, CurvePoint)

    e = miller(a, b)
    ret = final_exp(e)

    if a.is_infinite() or b.is_infinite():
        return Gfp12(GFP6_ZERO, GFP6_ONE)

    return ret


class Gt(object):
    p: Gfp12

    def __init__(self, p: Gfp12):
        self.p = p

    def scalar_mult(self, k: int) -> "Gt":
        return Gt(self.p.exp(k))

    def add(self, other: "Gt") -> "Gt":
        return Gt(self.p.mul(other.p))

    def marshal(self) -> bytes:
        return b"".join(
            map(
                lambda x: x.to_bytes(32, byteorder="big"),
                [
                    self.p.x.x.x.value(),
                    self.p.x.x.y.value(),
                    self.p.x.y.x.value(),
                    self.p.x.y.y.value(),
                    self.p.x.z.x.value(),
                    self.p.x.z.y.value(),
                    self.p.y.x.x.value(),
                    self.p.y.x.y.value(),
                    self.p.y.y.x.value(),
                    self.p.y.y.y.value(),
                    self.p.y.z.x.value(),
                    self.p.y.z.y.value(),
                ],
            ),
        )

    @classmethod
    def from_pair(cls, g1: G1, g2: G2) -> "Gt":
        p = optimal_ate(g2.p, g1.p)
        return cls(p)

    @classmethod
    def unmarshal(cls, data: bytes) -> "Gt":
        num_bytes = 32
        p0_x = int.from_bytes(data[0:num_bytes], "big")
        p0_y = int.from_bytes(data[num_bytes : num_bytes * 2], "big")
        p1_x = int.from_bytes(data[num_bytes * 2 : num_bytes * 3], "big")
        p1_y = int.from_bytes(data[num_bytes * 3 : num_bytes * 4], "big")
        p2_x = int.from_bytes(data[num_bytes * 4 : num_bytes * 5], "big")
        p2_y = int.from_bytes(data[num_bytes * 5 : num_bytes * 6], "big")
        p3_x = int.from_bytes(data[num_bytes * 6 : num_bytes * 7], "big")
        p3_y = int.from_bytes(data[num_bytes * 7 : num_bytes * 8], "big")
        p4_x = int.from_bytes(data[num_bytes * 8 : num_bytes * 9], "big")
        p4_y = int.from_bytes(data[num_bytes * 9 : num_bytes * 10], "big")
        p5_x = int.from_bytes(data[num_bytes * 10 : num_bytes * 11], "big")
        p5_y = int.from_bytes(data[num_bytes * 11 : num_bytes * 12], "big")

        p = Gfp12(
            Gfp6(Gfp2(p0_x, p0_y), Gfp2(p1_x, p1_y), Gfp2(p2_x, p2_y)),
            Gfp6(Gfp2(p3_x, p3_y), Gfp2(p4_x, p4_y), Gfp2(p5_x, p5_y)),
        )
        return cls(p)

    def hash(self) -> bytes:
        sha = hashlib.sha512()

        for parts in self.marshal():
            parts = parts.to_bytes()
            sha.update(parts)

        return sha.digest()
