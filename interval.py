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
TRIG_PAD = 1e-13


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
    ca, cb = np.cos(a), np.cos(b)
    lo = np.minimum(ca, cb)
    hi = np.maximum(ca, cb)
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
    sa, sb = np.sin(a), np.sin(b)
    lo = np.minimum(sa, sb)
    hi = np.maximum(sa, sb)
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
