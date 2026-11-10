"""Second, independent interval evaluation of the quantities of Section 6 (s = 1).

Written from the manuscript only, during the review of 9 November 2026, by a different AI
system (Claude Fable 5.2) from the one that wrote the main verifier.  Nothing here imports
interval.py, ad.py, verify.py or certify.py.  Differences from that implementation:
  * own interval class and complex helpers (pairs of intervals);
  * no automatic differentiation and no mean-value forms: plain interval evaluation;
  * own sin/cos enclosures: Taylor coefficients 1/n! built by interval division, degree
    23/22, remainders |r|^25/25! and |r|^24/24!; arguments outside [-1.5, 7.8], which lies
    inside (-pi/2, 5pi/2), are refused;
  * Z is evaluated with the double angles 2*theta_3 and 2*theta_2;
  * |Z|^2 is compared with an interval enclosure of 1/9, without the 1e-12 margins.
Shared with the main verifier: IEEE 754 binary64 round-to-nearest for + - * / sqrt in NumPy
(every result is widened by one ulp on each side), and the formulas of the manuscript.
That assumption is tested by environment_check() below, independently of interval.py, when
this module is imported; the import fails if it does not hold.
"""
import struct
import numpy as np

INF = np.inf


def dn(x):
    return np.nextafter(x, -INF)


def up(x):
    return np.nextafter(x, INF)


# ------------------------------------------------------------------ the arithmetic environment
def _word(x):
    return struct.unpack('<Q', struct.pack('<d', float(x)))[0]


# operation, operands, correctly rounded result, all as hexadecimal literals; checked from
# exact rationals when this file was written.  Rounding up, down or toward zero, flushing
# subnormal results (FTZ) or reading subnormal operands as zero (DAZ) each break some row.
ENVIRONMENT_CASES = [
    ('mul', '0x1.0000002p+0', '0x1.0000002p+0', '0x1.0000004p+0'),                # 1+2^-26+2^-54
    ('add', '0x1p+53', '0x1p+0', '0x1p+53'),                                      # tie, to even
    ('add', '0x1p+53', '0x1.8p+1', '0x1.0000000000002p+53'),                      # tie, to even
    ('sub', '0x1p+0', '0x1p-54', '0x1p+0'),                                       # tie below 1
    ('sub', '-0x1p+0', '0x1p-54', '-0x1p+0'),                                     # same, negative
    ('div', '0x1p+0', '0x1.4p+3', '0x1.999999999999ap-4'),                       # 1/10, rounds up
    ('div', '0x1.cp+2', '0x1.8p+1', '0x1.2aaaaaaaaaaabp+1'),                     # 7/3, rounds up
    ('div', '-0x1.cp+2', '0x1.8p+1', '-0x1.2aaaaaaaaaaabp+1'),                   # -7/3
    ('sqrt', '0x1.8p+1', None, '0x1.bb67ae8584caap+0'),                          # sqrt(3)
    ('sqrt', '0x1p-1', None, '0x1.6a09e667f3bcdp-1'),                            # sqrt(1/2)
    ('mul', '0x0.0000000000001p-1022', '0x1p+1', '0x0.0000000000002p-1022'),     # DAZ
    ('div', '0x1p-1022', '0x1p+1', '0x0.8p-1022'),                               # FTZ
    ('sub', '0x1.0000000000001p-1022', '0x1p-1022', '0x0.0000000000001p-1022'),  # FTZ
    ('up', '0x0p+0', None, '0x0.0000000000001p-1022'),
    ('dn', '0x1p+0', None, '0x1.fffffffffffffp-1'),
]
_OPS = {'add': np.add, 'sub': np.subtract, 'mul': np.multiply, 'div': np.divide,
        'sqrt': np.sqrt, 'up': up, 'dn': dn}


def environment_check(length=997):
    """Raise ArithmeticError unless every row of ENVIRONMENT_CASES holds bit for bit, on
    scalars and on arrays (whose NumPy loops may differ from the scalar ones)."""
    bad = []
    for op, a, b, want in ENVIRONMENT_CASES:
        args = [float.fromhex(a)] + ([] if b is None else [float.fromhex(b)])
        w = _word(float.fromhex(want))
        scalar = _OPS[op](*[np.float64(x) for x in args])
        array = _OPS[op](*[np.full(length, x) for x in args])
        if _word(scalar) != w or any(_word(v) != w for v in array):
            bad.append(f'{op}({", ".join(x for x in (a, b) if x)})')
    if bad:
        raise ArithmeticError('evaluator.py: the arithmetic is not IEEE binary64 with round to nearest '
                              'and gradual underflow; failing cases: ' + ', '.join(bad))
    return True


environment_check()


class Iv:
    __slots__ = ('a', 'b')

    def __init__(self, a, b=None):
        self.a = np.asarray(a, dtype=np.float64)
        self.b = self.a if b is None else np.asarray(b, dtype=np.float64)

    def __add__(self, o):
        o = iv(o)
        return Iv(dn(self.a + o.a), up(self.b + o.b))
    __radd__ = __add__

    def __neg__(self):
        return Iv(-self.b, -self.a)

    def __sub__(self, o):
        o = iv(o)
        return Iv(dn(self.a - o.b), up(self.b - o.a))

    def __rsub__(self, o):
        return iv(o) - self

    def __mul__(self, o):
        o = iv(o)
        p1, p2, p3, p4 = self.a * o.a, self.a * o.b, self.b * o.a, self.b * o.b
        return Iv(dn(np.minimum(np.minimum(p1, p2), np.minimum(p3, p4))), up(np.maximum(np.maximum(p1, p2), np.maximum(p3, p4))))
    __rmul__ = __mul__

    def sqr(self):
        x, y = self.a * self.a, self.b * self.b
        lo = np.where((self.a <= 0) & (self.b >= 0), 0.0, np.minimum(x, y))   # NaN propagates through minimum
        return Iv(np.maximum(dn(lo), 0.0), up(np.maximum(x, y)))

    def inv(self):
        ok = (self.a > 0) | (self.b < 0)          # False for NaN and for intervals containing 0
        lo = np.where(ok, dn(1.0 / np.where(ok, self.b, 1.0)), -INF)
        hi = np.where(ok, up(1.0 / np.where(ok, self.a, 1.0)), INF)
        return Iv(lo, hi)

    def __truediv__(self, o):
        return self * iv(o).inv()

    def __rtruediv__(self, o):
        return iv(o) * self.inv()

    def sqrt(self):
        return Iv(np.maximum(dn(np.sqrt(np.maximum(self.a, 0.0))), 0.0), up(np.sqrt(self.b)))


def iv(x):
    return x if isinstance(x, Iv) else Iv(x)


# ------------------------------------------------------------------ complex helpers (re, im)
def cadd(z, w):
    return (z[0] + w[0], z[1] + w[1])


def csub(z, w):
    return (z[0] - w[0], z[1] - w[1])


def cmul(z, w):
    return (z[0] * w[0] - z[1] * w[1], z[0] * w[1] + z[1] * w[0])


def cconj(z):
    return (z[0], -z[1])


def cscale(z, r):
    return (z[0] * r, z[1] * r)


def cabs2(z):
    return z[0].sqr() + z[1].sqr()


def csqr(z):
    return (z[0].sqr() - z[1].sqr(), (z[0] * z[1]) * 2.0)


def cneg(z):
    return (-z[0], -z[1])


def rdot(z, w):
    return z[0] * w[0] + z[1] * w[1]


# ------------------------------------------------------------------ trigonometry
PI = Iv(3.141592653589793, float(np.nextafter(3.141592653589793, 4.0)))   # float(pi) < pi < next
HALF_PI = PI * 0.5
_FACT_INV = [Iv(1.0)]
for _n in range(1, 26):
    _FACT_INV.append(_FACT_INV[-1] / float(_n))        # rigorous enclosures of 1/n!


def sincos_points(x):
    x = np.asarray(x, dtype=np.float64)
    k = np.rint(x / 1.5707963267948966)
    r = Iv(x) - HALF_PI * Iv(k)
    r2 = r.sqr()
    zero = np.zeros_like(x)
    S = Iv(zero); Cc = Iv(zero)
    for j in range(11, -1, -1):
        sgn = 1.0 if j % 2 == 0 else -1.0
        S = S * r2 + _FACT_INV[2 * j + 1] * sgn
        Cc = Cc * r2 + _FACT_INV[2 * j] * sgn
    S = S * r
    m = Iv(np.maximum(np.abs(r.a), np.abs(r.b)))
    m2 = m.sqr(); m4 = m2.sqr(); m8 = m4.sqr(); m16 = m8.sqr()
    m24 = m16 * m8
    rs = ((m24 * m) * _FACT_INV[25]).b
    rc = (m24 * _FACT_INV[24]).b
    S = S + Iv(-rs, rs); Cc = Cc + Iv(-rc, rc)
    q = np.mod(k, 4.0)
    s_lo = np.select([q == 0, q == 1, q == 2, q == 3], [S.a, Cc.a, -S.b, -Cc.b])
    s_hi = np.select([q == 0, q == 1, q == 2, q == 3], [S.b, Cc.b, -S.a, -Cc.a])
    c_lo = np.select([q == 0, q == 1, q == 2, q == 3], [Cc.a, -S.b, -Cc.b, S.a])
    c_hi = np.select([q == 0, q == 1, q == 2, q == 3], [Cc.b, -S.a, -Cc.a, S.b])
    return Iv(s_lo, s_hi), Iv(c_lo, c_hi)


_ARG_LO, _ARG_HI = -1.5, 7.8       # inside (-pi/2, 5pi/2): the only extrema are those handled below
_PAD = 1e-9


def sincos(x):
    """enclosures of sin and cos over interval arguments inside [-1.5, 7.8] (asserted)."""
    a, b = x.a, x.b
    if not (np.all(a >= _ARG_LO) and np.all(b <= _ARG_HI) and np.all(a <= b)):
        raise ArithmeticError('trig argument outside the admissible range')
    sa, ca = sincos_points(a)
    sb, cb = sincos_points(b)
    s_lo = np.minimum(sa.a, sb.a); s_hi = np.maximum(sa.b, sb.b)
    c_lo = np.minimum(ca.a, cb.a); c_hi = np.maximum(ca.b, cb.b)
    pl, ph = PI.a, PI.b

    def hits(c_lo_, c_hi_):
        return (a <= c_hi_ + _PAD) & (b >= c_lo_ - _PAD)
    s_hi = np.where(hits(0.5 * pl, 0.5 * ph), 1.0, s_hi)        # max of sin at pi/2
    s_lo = np.where(hits(1.5 * pl, 1.5 * ph), -1.0, s_lo)       # min of sin at 3pi/2
    c_hi = np.where(hits(0.0, 0.0) | hits(2 * pl, 2 * ph), 1.0, c_hi)   # max of cos at 0, 2pi
    c_lo = np.where(hits(pl, ph), -1.0, c_lo)                   # min of cos at pi
    return Iv(np.maximum(s_lo, -1.0), np.minimum(s_hi, 1.0)), Iv(np.maximum(c_lo, -1.0), np.minimum(c_hi, 1.0))


# ------------------------------------------------------------------ the chart
TWO_PI = PI * 2.0
NINTH = Iv(1.0) / 9.0          # enclosure of kappa^2 = 1/9  (s = 1)


def angles(y1, y2):
    base = (TWO_PI - y1) * 0.5
    half = (y1 * 3.0 - TWO_PI) * 0.5
    ph = y2 * half
    return base + ph, base - ph, y1


def zeta(y1, y2, l1, l2):
    th1, th2, th3 = angles(y1, y2)
    l3 = 1.0 - l1 - l2
    s23, c23 = sincos(th3 * 2.0)        # v2^2 = e^{2 i theta3}
    s22, c22 = sincos(th2 * 2.0)        # v3^2 = e^{-2 i theta2}
    Z = (l1 + l2 * c23 + l3 * c22, l2 * s23 - l3 * s22)
    return Z, (th1, th2, th3, l3, (c23, s23), (c22, -s22))


def quantities(y1, y2, l1, l2):
    """all rescaled quantities of Section 6 (s = 1) by naive interval evaluation."""
    Z, (th1, th2, th3, l3, w2, w3) = zeta(y1, y2, l1, l2)
    s1, c1 = sincos(th1); s2, c2 = sincos(th2); s3, c3 = sincos(th3)
    n = y1.a.shape
    one = Iv(np.ones(n)); zero = Iv(np.zeros(n))
    v = [(one, zero), (c3, s3), (c2, -s2)]
    v2 = [(one, zero), w2, w3]
    v3 = [(one, zero), cmul(w2, v[1]), cmul(w3, v[2])]
    v4 = [(one, zero), csqr(w2), csqr(w3)]
    lam = [l1, l2, l3]
    r1 = l1 / s1; r2 = l2 / s2
    R = [s3 * r1, s3 * r2, l3]                       # sigma3 * rho_i
    P12 = s3 * r1 * r2; P13 = l3 * r1; P23 = l3 * r2   # sigma3 * rho_a * rho_b
    half_i = lambda z: (z[1] * (-0.5), z[0] * 0.5)   # (i/2) z
    M = [half_i(csub(cscale(cconj(v[1]), P13), cscale(cconj(v[2]), P12))),
         half_i(csub(cscale(cconj(v[2]), P12), cscale(cconj(v[0]), P23))),
         half_i(csub(cscale(cconj(v[0]), P23), cscale(cconj(v[1]), P13)))]
    A = None; B = None; Zzz = None; Zzb = None; Zbb = None
    for i in range(3):
        lr = lam[i] * R[i]
        lrr = lam[i] * R[i].sqr()
        a = csub(cmul(M[i], v2[i]), cscale(v[i], lr))
        b = cadd(cmul(cconj(M[i]), v2[i]), cscale(v3[i], lr))
        zz = cscale(cmul(M[i], v[i]), R[i] * (-2.0))
        zb = csub(csub(cscale(cmul(M[i], v3[i]), R[i]), cscale(cmul(cconj(M[i]), v[i]), R[i])), cscale(v2[i], lrr))
        bb = cadd(cscale(cmul(cconj(M[i]), v3[i]), R[i] * 2.0), cscale(v4[i], lrr * 2.0))
        A = a if A is None else cadd(A, a); B = b if B is None else cadd(B, b)
        Zzz = zz if Zzz is None else cadd(Zzz, zz); Zzb = zb if Zzb is None else cadd(Zzb, zb); Zbb = bb if Zbb is None else cadd(Zbb, bb)
    jac = cabs2(A) - cabs2(B)
    z2 = cabs2(Z)
    jinv = jac.inv()
    iZ = (-Z[1], Z[0])
    t = cscale(csub(cmul(cconj(A), iZ), cmul(B, cconj(iZ))), jinv)          # t / sigma3
    t2 = csqr(t); tt = cabs2(t)
    Rr = cneg(cadd(Z, cadd(cadd(cmul(Zzz, t2), cscale(Zzb, tt * 2.0)), cmul(Zbb, cconj(t2)))))
    tp = cscale(csub(cmul(cconj(A), Rr), cmul(B, cconj(Rr))), jinv)
    nu = cmul(tp, cconj(t))[1] / tt
    cz = cmul(cconj(t2), Z)
    W = r1 * s3 * t[1]
    kb = cscale((-Z[0], Z[1]), z2.sqrt().inv())                              # conj(k)^2 = -conj(Z)/|Z|

    def Q(p, q):
        return rdot(p, q) * 0.5 + cmul(cmul(p, q), kb)[0] * 0.5
    Ich = s1 * s2 * 0.5 + Q(v[0], v[2]) * 3.0 * (l1 * s3 * rdot(v[0], t) - s1 * l3 * rdot(v[2], t)) \
        - l1 * s3 * Q(v[2], t) + s1 * l3 * Q(v[0], t)
    return dict(z2=z2, jac=jac, Ich=Ich, nu=nu, czre=cz[0], czim=cz[1], W=W)
