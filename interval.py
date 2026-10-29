"""Minimal vectorised interval arithmetic with outward rounding (numpy float64).

Every arithmetic result is widened by one ulp in each direction (np.nextafter), which encloses
the exact result of the IEEE round-to-nearest operation.  sin/cos are evaluated at the
endpoints and widened by an absolute 1e-13 (far larger than the libm error), with interior
extrema handled explicitly.  Complex intervals are rectangles (pairs of real intervals).
"""
import numpy as np

_dn = lambda x: np.nextafter(x, -np.inf)
_up = lambda x: np.nextafter(x, np.inf)

PI_LO = 3.141592653589793      # float(pi) < pi
PI_HI = np.nextafter(PI_LO, 4)  # > pi
TRIG_PAD = 0.0  # endpoint values are now rigorous (sincos_point); no libm assumption


class I:
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        lo = np.asarray(lo, dtype=float)
        self.lo = lo
        self.hi = lo if hi is None else np.asarray(hi, dtype=float)

    # ---- helpers
    @staticmethod
    def c(x):
        return x if isinstance(x, I) else I(x, x)

    def _defer(self, o):
        return type(o).__name__ == 'DI'

    def __repr__(self):
        return f'I({self.lo}, {self.hi})'

    @property
    def mid(self):
        return 0.5 * (self.lo + self.hi)

    @property
    def wid(self):
        return self.hi - self.lo

    # ---- arithmetic
    def __add__(self, o):
        if self._defer(o):
            return NotImplemented
        o = I.c(o)
        return I(_dn(self.lo + o.lo), _up(self.hi + o.hi))
    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, o):
        if self._defer(o):
            return NotImplemented
        o = I.c(o)
        return I(_dn(self.lo - o.hi), _up(self.hi - o.lo))

    def __rsub__(self, o):
        return I.c(o) - self

    def __mul__(self, o):
        if self._defer(o):
            return NotImplemented
        if not isinstance(o, I):
            o = np.asarray(o, dtype=float)
            a, b = self.lo * o, self.hi * o
            return I(_dn(np.minimum(a, b)), _up(np.maximum(a, b)))
        p1, p2, p3, p4 = self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi
        lo = np.minimum(np.minimum(p1, p2), np.minimum(p3, p4))
        hi = np.maximum(np.maximum(p1, p2), np.maximum(p3, p4))
        return I(_dn(lo), _up(hi))
    __rmul__ = __mul__

    def sq(self):
        a, b = self.lo * self.lo, self.hi * self.hi
        lo = np.where((self.lo <= 0) & (self.hi >= 0), 0.0, np.minimum(a, b))
        return I(np.maximum(_dn(lo), 0.0), _up(np.maximum(a, b)))

    def recip(self):
        bad = (self.lo <= 0) & (self.hi >= 0)
        lo = np.where(bad, -np.inf, _dn(1.0 / np.where(bad, 1.0, self.hi)))
        hi = np.where(bad, np.inf, _up(1.0 / np.where(bad, 1.0, self.lo)))
        return I(lo, hi)

    def __truediv__(self, o):
        if self._defer(o):
            return NotImplemented
        if not isinstance(o, I):
            return self * I(1.0 / np.asarray(o, float)) if False else self * I.c(o).recip()
        return self * o.recip()

    def __rtruediv__(self, o):
        return I.c(o) * self.recip()

    def sqrt(self):
        return I(np.maximum(_dn(np.sqrt(np.maximum(self.lo, 0.0))), 0.0), _up(np.sqrt(self.hi)))

    def __getitem__(self, k):
        return I(self.lo[k], self.hi[k])


def cos_i(x):
    """cos on intervals contained in [-0.1, 2*pi+0.1]"""
    if type(x).__name__ == 'DI':
        return x.cos()
    a, b = x.lo, x.hi
    _, cA = sincos_point(a)
    _, cB = sincos_point(b)
    lo = np.minimum(cA.lo, cB.lo)
    hi = np.maximum(cA.hi, cB.hi)
    # minimum -1 at pi; maximum 1 at 0 and 2pi
    eps = 1e-9  # generous: including an extremum that is not there only widens the enclosure
    has_pi = (a <= PI_HI + eps) & (b >= PI_LO - eps)
    has_0 = (a <= eps) & (b >= -eps)
    has_2pi = (a <= 2 * PI_HI + eps) & (b >= 2 * PI_LO - eps)
    lo = np.where(has_pi, -1.0, lo - TRIG_PAD)
    hi = np.where(has_0 | has_2pi, 1.0, hi + TRIG_PAD)
    return I(np.maximum(lo, -1.0), np.minimum(hi, 1.0))


def sin_i(x):
    """sin on intervals contained in [-0.1, 2*pi+0.1]"""
    if type(x).__name__ == 'DI':
        return x.sin()
    a, b = x.lo, x.hi
    sA, _ = sincos_point(a)
    sB, _ = sincos_point(b)
    lo = np.minimum(sA.lo, sB.lo)
    hi = np.maximum(sA.hi, sB.hi)
    eps = 1e-9
    has_p2 = (a <= PI_HI / 2 + eps) & (b >= PI_LO / 2 - eps)
    has_3p2 = (a <= 1.5 * PI_HI + eps) & (b >= 1.5 * PI_LO - eps)
    hi = np.where(has_p2, 1.0, hi + TRIG_PAD)
    lo = np.where(has_3p2, -1.0, lo - TRIG_PAD)
    return I(np.maximum(lo, -1.0), np.minimum(hi, 1.0))


class C:
    """complex interval (rectangle)"""
    __slots__ = ('re', 'im')

    def __init__(self, re, im=None):
        self.re = re if type(re).__name__ == 'DI' else I.c(re)
        if im is None:
            self.im = 0.0 * self.re
        else:
            self.im = im if type(im).__name__ == 'DI' else I.c(im)

    @staticmethod
    def c(x):
        if isinstance(x, C):
            return x
        if type(x).__name__ == 'DI':
            return C(x, 0.0 * x)
        if isinstance(x, I):
            return C(x, I(np.zeros_like(x.lo)))
        x = complex(x)
        return C(I(x.real), I(x.imag))

    def __add__(self, o):
        o = C.c(o)
        return C(self.re + o.re, self.im + o.im)
    __radd__ = __add__

    def __sub__(self, o):
        o = C.c(o)
        return C(self.re - o.re, self.im - o.im)

    def __rsub__(self, o):
        return C.c(o) - self

    def __neg__(self):
        return C(-self.re, -self.im)

    def __mul__(self, o):
        if isinstance(o, I) or type(o).__name__ == 'DI':
            return C(self.re * o, self.im * o)
        if not isinstance(o, C):
            o = complex(o)
            if o.imag == 0:
                return C(self.re * o.real, self.im * o.real)
            o = C.c(o)
        return C(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    __rmul__ = __mul__

    def conj(self):
        return C(self.re, -self.im)

    def abs2(self):
        return self.re.sq() + self.im.sq()

    def sq(self):
        return C(self.re.sq() - self.im.sq(), 2.0 * (self.re * self.im))

    def __truediv__(self, o):
        if isinstance(o, I) or type(o).__name__ == 'DI':
            r = o.recip()
            return C(self.re * r, self.im * r)
        o = C.c(o)
        d = o.abs2().recip()
        return C(self.re, self.im) * o.conj() * d


# ---------------------------------------------------------------------------
# Rigorous point enclosures of sin / cos using only IEEE +,-,* with outward rounding
# (no reliance on the accuracy of the platform libm).
_HALF_PI = None


def _poly_sin_cos(y):
    """y: interval with |y| <= 0.8. Returns (sin y, cos y) enclosures via Taylor series
    to degree 21/20 plus a rigorous Lagrange remainder bound."""
    y2 = y.sq()
    # sin y = y * sum_{k=0}^{10} (-1)^k y^{2k}/(2k+1)!,  remainder <= |y|^23/23!
    # cos y =     sum_{k=0}^{10} (-1)^k y^{2k}/(2k)!,    remainder <= |y|^22/22!
    import math
    ss = I(np.zeros_like(y.lo)); cc = I(np.zeros_like(y.lo))
    for k in range(10, -1, -1):
        cs = (-1) ** k / math.factorial(2 * k + 1)
        cc_ = (-1) ** k / math.factorial(2 * k)
        # coefficients are not exactly representable: enclose them
        csI = I(np.nextafter(cs, -np.inf), np.nextafter(cs, np.inf))
        ccI = I(np.nextafter(cc_, -np.inf), np.nextafter(cc_, np.inf))
        ss = ss * y2 + csI
        cc = cc * y2 + ccI
    ss = ss * y
    m = np.maximum(np.abs(y.lo), np.abs(y.hi))
    rs = _up(_up(m ** 23) / 2.5e22)    # 23! = 2.585e22 > 2.5e22
    rc = _up(_up(m ** 22) / 1.1e21)    # 22! = 1.124e21 > 1.1e21
    return ss + I(-rs, rs), cc + I(-rc, rc)


def sincos_point(x):
    """rigorous enclosures of sin(x), cos(x) for float arrays x in [-1, 2*pi+1]"""
    x = np.asarray(x, dtype=float)
    hp = I(PI_LO / 2, PI_HI / 2)
    k = np.round(x / (PI_LO / 2)).astype(int)
    y = I(x, x) - hp * k.astype(float)
    sy, cy = _poly_sin_cos(y)
    km = np.mod(k, 4)
    # sin(x) = sin(y + k pi/2), cos(x) = cos(y + k pi/2)
    s_lo = np.select([km == 0, km == 1, km == 2, km == 3], [sy.lo, cy.lo, -sy.hi, -cy.hi])
    s_hi = np.select([km == 0, km == 1, km == 2, km == 3], [sy.hi, cy.hi, -sy.lo, -cy.lo])
    c_lo = np.select([km == 0, km == 1, km == 2, km == 3], [cy.lo, -sy.hi, -cy.hi, sy.lo])
    c_hi = np.select([km == 0, km == 1, km == 2, km == 3], [cy.hi, -sy.lo, -cy.lo, sy.hi])
    return I(s_lo, s_hi), I(c_lo, c_hi)
