import copy

from bn256.constants import U, P
from bn256.curve import CurvePoint
from bn256.gfp12 import Gfp12
from bn256.gfp2 import (Gfp2, XI_TO_P_MINUS_1_OVER_2, XI_TO_P_MINUS_1_OVER_3, XI_TO_P_SQUARED_MINUS_1_OVER_3)
from bn256.gfp6 import Gfp6
from bn256.twist import TwistPoint

XI = Gfp2(1, 3)  # i + 3

XI1 = [
    XI.exp(1 * (P - 1) // 6),
    XI.exp(2 * (P - 1) // 6),
    XI.exp(3 * (P - 1) // 6),
    XI.exp(4 * (P - 1) // 6),
    XI.exp(5 * (P - 1) // 6),
]

XI2 = [(x * x.conjugate()) for x in XI1]

# 6u+2 in NAF
NAF_6UP2 = [1, 1, 0, 1, 0, 0, -1, 0, 1, 1, 0, 0, 0, -1, 0, 0, 1, 1, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0, 0,
            1, 1, 1, 0, 0, 0, 0, -1, 0, 1, 0, 0, -1, 0, 1, 1, 0, 0, 1, 0, 0, -1, 1, 0, 0, -1, 0, 1, 0, 1, 0, 0, 0]


def line_func_add(r: TwistPoint, p: TwistPoint, q: CurvePoint, r2: Gfp2):  # ✅
    assert isinstance(r, TwistPoint)
    assert isinstance(p, TwistPoint)
    assert isinstance(q, CurvePoint)
    assert isinstance(r2, Gfp2)

    B = p.x * r.t
    D = ((p.y + r.z).square() - r2 - r.t) * r.t
    H = B - r.x
    I = H.square()
    E = I.double().double()
    J = H * E
    L1 = D - r.y - r.y
    V = r.x * E

    r_out_x = L1.square() - J - V.double()
    r_out_y = (V - r_out_x) * L1 - (r.y * J).double()
    r_out_z = (r.z + H).square() - r.t - I
    r_out_t = r_out_z.square()
    r_out = TwistPoint(r_out_x, r_out_y, r_out_z, r_out_t)

    a = (L1 * p.x).double() - ((p.y + r_out_z).square() - r2 - r_out_t)
    b = L1.negative().mul_scalar(q.x).double()
    c = r_out_z.mul_scalar(q.y).double()

    return a, b, c, r_out


def line_func_double(r: TwistPoint, q: CurvePoint):
    assert isinstance(r, TwistPoint)
    assert isinstance(q, CurvePoint)
    # A = x²
    A = r.x.square()
    # B = y²
    B = r.y.square()
    # C = y⁴
    C = B.square()
    # D = 2((x + y²)² - x² - y⁴)
    D = ((r.x + B).square() - A - C).double()
    # E = 3yx²
    E = A + A + A
    # G = 9y²x⁴
    G = E.square()

    # R.x = 9y²x⁴ - 4((x + y²)² - x² - y⁴)
    r_out_x = G - D.double()
    # R.y = (2((x + y²)² - x² - y⁴) - (9y²x⁴ - 4((x + y²)² - x² - y⁴))) * 3yx² - 8y²
    r_out_y = (D - r_out_x) * E - C.double().double().double()
    # R.z = (y + z)² - y² - t
    r_out_z = (r.y + r.z).square() - B - r.t
    # R.t = ((y + z)² - y² - t)²
    r_out_t = r_out_z.square()
    r_out = TwistPoint(r_out_x, r_out_y, r_out_z, r_out_t)

    # a = (x + 3yx²)² - x² - 9y²x⁴ - 4y²
    a = (r.x + E).square() - A - G - B.double().double()
    # b = (0 - 6tyx²) * k1

    b = (Gfp2.zero() - (E * r.t).double()).mul_scalar(q.x)
    # c = 2t((y + z)² - y² - t) * k2
    c = (r_out_z * r.t).double().mul_scalar(q.y)

    return a, b, c, r_out


def mul_line(r: Gfp12, a: Gfp2, b: Gfp2, c: Gfp2):  # ✅
    assert isinstance(r, Gfp12)
    assert isinstance(a, Gfp2)
    assert isinstance(b, Gfp2)
    assert isinstance(c, Gfp2)

    a2 = Gfp6(Gfp2.zero(), a, b) * r.x
    t2 = Gfp6(Gfp2.zero(), a, b + c)
    yc = r.y.mul_scalar(c)
    r.x = (r.x + r.y) * t2 - a2 - yc
    r.y = yc + a2.mul_tau()


def miller(q: TwistPoint, p: CurvePoint):  # ❌
    assert isinstance(q, TwistPoint)
    assert isinstance(p, CurvePoint)

    ret = Gfp12.one()

    a_affine = q.copy().make_affine()
    b_affine = p.copy().make_affine()
    minus_a = -a_affine

    r = a_affine.copy()
    r2 = a_affine.y.square()

    for index, naf_i in enumerate(NAF_6UP2[:-1]):
        a, b, c, new_r = line_func_double(r, b_affine)
        if index != 0:
            ret = ret.square()
        mul_line(ret, a, b, c)

        r = new_r
        next_naf_i = NAF_6UP2[index + 1]
        if next_naf_i == 1:
            a, b, c, new_r = line_func_add(r, a_affine, b_affine, r2)
        elif next_naf_i == -1:
            a, b, c, new_r = line_func_add(r, minus_a, b_affine, r2)
        else:
            continue
        mul_line(ret, a, b, c)
        r = new_r

    # In order to calculate Q1 we have to convert q from the sextic twist
    # to the full GF(p^12) group, apply the Frobenius there, and convert
    # back.
    #
    # The twist isomorphism is (x', y') -> (xω², yω³). If we consider just
    # x for a moment, then after applying the Frobenius, we have x̄ω^(2p)
    # where x̄ is the conjugate of x. If we are going to apply the inverse
    # isomorphism we need a value with a single coefficient of ω² so we
    # rewrite this as x̄ω^(2p-2)ω². ξ⁶ = ω and, due to the construction of
    # p, 2p-2 is a multiple of six. Therefore we can rewrite as
    # x̄ξ^((p-1)/3)ω² and applying the inverse isomorphism eliminates the
    # ω².
    #
    # A similar argument can be made for the y value.

    q1 = TwistPoint(
        a_affine.x.conjugate().mul(XI_TO_P_MINUS_1_OVER_3),
        a_affine.y.conjugate().mul(XI_TO_P_MINUS_1_OVER_2),
        Gfp2.one(),
        Gfp2.one()
    )  # ✅

    # For Q2 we are applying the p² Frobenius. The two conjugations cancel
    # out and we are left only with the factors from the isomorphism. In
    # the case of x, we end up with a pure number which is why
    # xiToPSquaredMinus1Over3 is ∈ GF(p). With y we get a factor of -1. We
    # ignore this to end up with -Q2.

    minus_q2 = TwistPoint(
        a_affine.x.mul_scalar(XI_TO_P_SQUARED_MINUS_1_OVER_3),
        a_affine.y,
        Gfp2.one(),
        Gfp2.one()
    )  # ✅

    r2 = q1.y.square()
    a, b, c, new_r = line_func_add(r, q1, b_affine, r2)  # ❌ ret
    mul_line(ret, a, b, c)
    r = new_r

    r2 = minus_q2.y.square()
    a, b, c, new_r = line_func_add(r, minus_q2, b_affine, r2)
    mul_line(ret, a, b, c)
    r = new_r
    return ret


def final_exp(inp: Gfp12, u: int = U):  # ✅
    """
    finalExponentiation computes the (p¹²-1)/Order-th power of an element of
    GF(p¹²) to obtain an element of GT (steps 13-15 of algorithm 1 from
    http://cryptojedi.org/papers/dclxvi-20100714.pdf)
    """
    assert isinstance(inp, Gfp12)

    inv = inp.invert()

    t1 = inp.conjugate() * inv
    t2 = t1.frobenius_p2()
    t1 = t1 * t2

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

    y0 = fp1 * fp2 * fp3
    y1 = t1.conjugate()
    y3 = y3.conjugate()
    y4 = fu1 * fu2p
    y4 = y4.conjugate()
    y5 = fu2.conjugate()
    y6 = (fu3 * fu3p).conjugate()

    t0 = y6.square() * y4 * y5
    t1 = ((y3 * y5 * t0).square() * t0 * y2).square()

    ret = (t1 * y1).square() * (t1 * y0)
    return ret


def optimal_ate(q: TwistPoint, p: CurvePoint) -> Gfp12:  # ✅
    assert isinstance(q, TwistPoint)
    assert isinstance(p, CurvePoint)

    e = miller(q, p)
    ret = final_exp(e)
    if q.is_infinity() or p.is_infinity():
        ret.set_one()
    return ret
