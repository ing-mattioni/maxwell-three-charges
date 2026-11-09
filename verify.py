"""Interval enclosures of the quantities behind conditions (J), (I), (T), (W).

Chart (WLOG by cyclic relabelling, theta3 is the largest angle):
    y1 = theta3 in [2pi/3, pi],  y2 = tau in [-1, 1]  (see thetas_from_y),
    lambda1, lambda2 >= 0, lambda3 = 1 - lambda1 - lambda2 >= 0.
enclose_y returns enclosures of the rescaled quantities of Section 6 of the paper:
    zeta2 = |zeta|^2,  zeta = sum lambda_i v_i^2
    jac   = sigma3^2 J                         (J): jac < 0 on { |zeta| <= kappa }
    Ih    = inflection function, eq. (6.9)     (I): Ih > 0 on { |zeta| = kappa }
    turn  = nu = d arg t / d omega             (T): turn < 1/2 at the cusps
    czre, czim = real and imaginary parts of conj(t)^2 zeta / sigma3^2
                 (a cusp needs czim = 0 and czre < 0)
    Wq    = rho_1 Im t                         (W): |Wq| < 1/2 on { |zeta| = kappa }
with kappa = s/(s+2).  This module only evaluates; the branch and bound, the mean-value
forms and the decision to close a box are in certify.py.
"""
import numpy as np
from interval import I, C, sin_i, cos_i, PI_LO, PI_HI

TWO_PI = I(2 * PI_LO, 2 * PI_HI)


def thetas_from_y(y1, y2):
    """rectangle coordinates: y1 = theta3 in [2pi/3, pi], y2 = tau in [-1, 1];
    theta1,2 = (2pi - theta3)/2 +- tau (3 theta3 - 2pi)/2   (so theta3 = max, sum = 2pi exactly)"""
    base = (TWO_PI - y1) * 0.5
    half = (3.0 * y1 - TWO_PI) * 0.5
    phi = y2 * half
    return base + phi, base - phi, y1


def enclose_y(y1, y2, l1, l2, s, need_bdry=True):
    th1, th2, th3 = thetas_from_y(y1, y2)
    return enclose(th1, th2, l1, l2, s, need_bdry, th3=th3)


def enclose(th1, th2, l1, l2, s, need_bdry=True, th3=None):
    if th3 is None:
        th3 = TWO_PI - th1 - th2
    l3 = 1.0 - l1 - l2
    s1, s2, s3 = sin_i(th1), sin_i(th2), sin_i(th3)
    c1, c2, c3 = cos_i(th1), cos_i(th2), cos_i(th3)
    n = th1.lo.shape[-1]
    zero = I(np.zeros(n))
    v = [C(I(np.ones(n)), zero), C(c3, s3), C(c2, -s2)]
    lam = [l1, l2, l3]
    rho = [l1 / s1, l2 / s2, None]
    rhoh = [s3 * rho[0], s3 * rho[1], l3]
    v2 = [vi.sq() for vi in v]
    zeta = v2[0] * l1 + v2[1] * l2 + v2[2] * l3
    zeta2 = zeta.abs2()

    def s3rr(a, b):
        if a == 2:
            return l3 * rho[b]
        if b == 2:
            return l3 * rho[a]
        return s3 * rho[a] * rho[b]

    muh = []
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        w = v[j].conj() * s3rr(i, k) - v[k].conj() * s3rr(i, j)
        muh.append(C(-0.5 * w.im, 0.5 * w.re))  # (i/2) * w
    Zz = None
    Zzb = None
    for i in range(3):
        a = muh[i] * v2[i] - v[i] * (lam[i] * rhoh[i])
        b = muh[i].conj() * v2[i] + (v2[i] * v[i]) * (lam[i] * rhoh[i])
        Zz = a if Zz is None else Zz + a
        Zzb = b if Zzb is None else Zzb + b
    jac = Zz.abs2() - Zzb.abs2()
    out = dict(zeta2=zeta2, jac=jac)
    if not need_bdry:
        return out
    Z = zeta
    iZ = C(-Z.im, Z.re)
    num = Zz.conj() * iZ - Zzb * iZ.conj()
    t = num / jac
    absz = zeta2.sqrt()
    Kb = C(-Z.re, Z.im) / absz          # -conj(zeta)/|zeta|
    dot = lambda a, b: a.re * b.re + a.im * b.im

    def Q(a, b):
        ab = a * b * Kb
        return 0.5 * dot(a, b) + 0.5 * ab.re

    v1, v3 = v[0], v[2]
    Q13 = Q(v1, v3)
    Ih = 0.5 * (s1 * s2) + (s + 2) * Q13 * (l1 * s3 * dot(v1, t) - s1 * l3 * dot(v3, t)) \
        - l1 * s3 * Q(v3, t) + s1 * l3 * Q(v1, t)
    Zzz = None; Zzzb = None; Zzbzb = None
    for i in range(3):
        a = (muh[i] * v[i]) * (-2.0 * rhoh[i])
        b = (muh[i] * (v2[i] * v[i])) * rhoh[i] - (muh[i].conj() * v[i]) * rhoh[i] - v2[i] * (lam[i] * rhoh[i].sq())
        c = (muh[i].conj() * (v2[i] * v[i])) * (2.0 * rhoh[i]) + v2[i].sq() * (2.0 * lam[i] * rhoh[i].sq())
        Zzz = a if Zzz is None else Zzz + a
        Zzzb = b if Zzzb is None else Zzzb + b
        Zzbzb = c if Zzbzb is None else Zzbzb + c
    rhs = -Z - (Zzz * t.sq() + Zzzb * (2.0 * t.abs2()) + Zzbzb * t.conj().sq())
    tp = (Zz.conj() * rhs - Zzb * rhs.conj()) / jac
    tpt = tp * t.conj()
    turn = tpt.im / t.abs2()
    # cusp locus: t parallel to k  <=>  conj(t)^2 * zeta is real and negative
    cz = t.conj().sq() * Z
    # (W): |d arg u_1 / d omega| = sigma3 * |Im(th * conj(v1))| * rho_1, with v1 = 1
    Wq = s3 * t.im * rho[0]
    out.update(Ih=Ih, turn=turn, czre=cz.re, czim=cz.im, Wq=Wq)
    return out
