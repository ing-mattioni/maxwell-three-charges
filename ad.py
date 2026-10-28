"""Forward-mode automatic differentiation over intervals (value + 4 partial derivatives).

DI.lo / DI.hi have shape (5, N): row 0 = value enclosure, rows 1..4 = enclosures of the
partial derivatives with respect to the 4 box coordinates.  All operations reuse the
outward-rounded interval primitives of interval.py, so the enclosures are rigorous.
"""
import numpy as np
import interval as iv
from interval import I

ND = 4


class DI:
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi):
        self.lo = lo
        self.hi = hi

    @staticmethod
    def var(x, k):
        """independent variable number k (0..3) on the interval x"""
        n = x.lo.shape[0]
        lo = np.zeros((ND + 1, n)); hi = np.zeros((ND + 1, n))
        lo[0], hi[0] = x.lo, x.hi
        lo[k + 1] = 1.0; hi[k + 1] = 1.0
        return DI(lo, hi)

    @staticmethod
    def const(x):
        if isinstance(x, DI):
            return x
        if isinstance(x, I):
            n = x.lo.shape
            lo = np.zeros((ND + 1,) + n); hi = np.zeros((ND + 1,) + n)
            lo[0], hi[0] = x.lo, x.hi
            return DI(lo, hi)
        return None

    def val(self):
        return I(self.lo[0], self.hi[0])

    def d(self, k):
        return I(self.lo[k + 1], self.hi[k + 1])

    def _all(self):
        return I(self.lo, self.hi)

    def _coerce(self, o):
        if isinstance(o, DI):
            return o
        if isinstance(o, I):
            n = self.lo.shape[1]
            lo = np.zeros_like(self.lo); hi = np.zeros_like(self.hi)
            lo[0] = np.broadcast_to(o.lo, (n,)); hi[0] = np.broadcast_to(o.hi, (n,))
            return DI(lo, hi)
        o = float(o)
        lo = np.zeros_like(self.lo); lo[0] = o
        return DI(lo, lo.copy())

    def __add__(self, o):
        o = self._coerce(o)
        r = self._all() + o._all()
        return DI(r.lo, r.hi)
    __radd__ = __add__

    def __neg__(self):
        return DI(-self.hi, -self.lo)

    def __sub__(self, o):
        o = self._coerce(o)
        r = self._all() - o._all()
        return DI(r.lo, r.hi)

    def __rsub__(self, o):
        return self._coerce(o) - self

    def __mul__(self, o):
        if not isinstance(o, (DI, I)):
            r = self._all() * float(o)
            return DI(r.lo, r.hi)
        o = self._coerce(o)
        a0, b0 = self.val(), o.val()
        v = a0 * b0
        # d = a0*b' + a'*b0  (broadcast over rows 1..4)
        A0 = I(np.broadcast_to(a0.lo, self.lo[1:].shape), np.broadcast_to(a0.hi, self.lo[1:].shape))
        B0 = I(np.broadcast_to(b0.lo, self.lo[1:].shape), np.broadcast_to(b0.hi, self.lo[1:].shape))
        dd = A0 * I(o.lo[1:], o.hi[1:]) + I(self.lo[1:], self.hi[1:]) * B0
        return DI(np.concatenate([v.lo[None], dd.lo]), np.concatenate([v.hi[None], dd.hi]))
    __rmul__ = __mul__

    def _chain(self, f0, fp):
        """value f0 (I), derivative factor fp (I): result' = fp * self'"""
        FP = I(np.broadcast_to(fp.lo, self.lo[1:].shape), np.broadcast_to(fp.hi, self.lo[1:].shape))
        dd = FP * I(self.lo[1:], self.hi[1:])
        return DI(np.concatenate([f0.lo[None], dd.lo]), np.concatenate([f0.hi[None], dd.hi]))

    def sq(self):
        a0 = self.val()
        return self._chain(a0.sq(), 2.0 * a0)

    def recip(self):
        a0 = self.val()
        r = a0.recip()
        return self._chain(r, -(r.sq()))

    def __truediv__(self, o):
        if not isinstance(o, (DI, I)):
            return self * (1.0 / float(o)) if False else self * self._coerce(o).recip()
        return self * self._coerce(o).recip()

    def __rtruediv__(self, o):
        return self._coerce(o) * self.recip()

    def sqrt(self):
        a0 = self.val()
        r = a0.sqrt()
        return self._chain(r, 0.5 * r.recip())

    def sin(self):
        a0 = self.val()
        return self._chain(iv.sin_i(a0), iv.cos_i(a0))

    def cos(self):
        a0 = self.val()
        return self._chain(iv.cos_i(a0), -iv.sin_i(a0))


def mvf(f_box, f_center, box_lo, box_hi):
    """mean value form enclosure intersected with the naive one.
    f_box: DI evaluated on the box; f_center: I evaluated at the centre point."""
    c = 0.5 * (box_lo + box_hi)
    acc = f_center
    for k in range(ND):
        dk = f_box.d(k)
        acc = acc + dk * I(box_lo[:, k] - c[:, k], box_hi[:, k] - c[:, k])
    nv = f_box.val()
    return I(np.maximum(acc.lo, nv.lo), np.minimum(acc.hi, nv.hi))
