"""Validate the second evaluator (evaluator.py) against 40-digit reference values.

Usage:  python independent/validate.py [BOXES]          (needs mpmath)

1. Point enclosures of sin and cos, and their interval extensions, against mpmath.
2. On random boxes of the chart (generic, on theta3 = pi, on lambda_i = 0), the enclosures of
   all seven quantities must contain the exact values at the centre, at two corners and at a
   random point.  The exact values come from first_principles.paper_formulas, a transcription
   of the manuscript's formulas that first_principles.py compares with the definitions.
Exit status 1 on any violation.
"""
import sys
import os
import math
import numpy as np
import mpmath as mp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import evaluator as ev
from evaluator import Iv
from first_principles import paper_formulas

mp.mp.dps = 40
rng = np.random.default_rng(99)
violations = 0

# ---- trigonometry: points
xs = np.concatenate([rng.uniform(-1.5, 7.8, 20000), np.array([k * math.pi / 4 for k in range(-1, 10)]),
                     np.nextafter(np.array([k * math.pi / 2 for k in range(0, 5)]), 10), np.nextafter(np.array([k * math.pi / 2 for k in range(0, 5)]), -10)])
S, C = ev.sincos_points(xs)
bad = 0; wmax = 0
for i, x in enumerate(xs):
    m = mp.mpf(float(x))
    bad += not (mp.mpf(float(S.a[i])) <= mp.sin(m) <= mp.mpf(float(S.b[i])))
    bad += not (mp.mpf(float(C.a[i])) <= mp.cos(m) <= mp.mpf(float(C.b[i])))
    wmax = max(wmax, float(S.b[i] - S.a[i]), float(C.b[i] - C.a[i]))
print('sin/cos at', len(xs), 'points: violations', bad, ' largest enclosure width', wmax)
violations += bad
# ---- trigonometry: intervals
lo = rng.uniform(-1.5, 7.7, 5000); hi = np.minimum(lo + 10.0 ** rng.uniform(-16, 0.5, 5000), 7.8)
Si, Ci = ev.sincos(Iv(lo, hi))
bad = 0
for i in range(len(lo)):
    a, b = mp.mpf(float(lo[i])), mp.mpf(float(hi[i]))
    pts = [a, b] + [e for e in (mp.pi / 2, 3 * mp.pi / 2, 0, mp.pi, 2 * mp.pi) if a <= e <= b] + [a + (b - a) * j / 9 for j in range(1, 9)]
    sv = [mp.sin(p) for p in pts]; cv = [mp.cos(p) for p in pts]
    bad += not (mp.mpf(float(Si.a[i])) <= min(sv) and max(sv) <= mp.mpf(float(Si.b[i])))
    bad += not (mp.mpf(float(Ci.a[i])) <= min(cv) and max(cv) <= mp.mpf(float(Ci.b[i])))
print('sin/cos on', len(lo), 'intervals: violations', bad)
violations += bad
try:
    ev.sincos(Iv(np.array([6.0]), np.array([8.0])))
    print('ERROR: out-of-range argument accepted'); violations += 1
except ArithmeticError:
    print('out-of-range argument refused, as intended')

# ---- the seven quantities
KEYS = ('z2', 'jac', 'Ich', 'nu', 'czre', 'czim', 'W')


def exact(p):
    y1, y2, l1, l2 = (mp.mpf(float(v)) for v in p)
    th1 = (2 * mp.pi - y1) / 2 + y2 * (3 * y1 - 2 * mp.pi) / 2
    th2 = (2 * mp.pi - y1) / 2 - y2 * (3 * y1 - 2 * mp.pi) / 2
    P = paper_formulas((th1, th2, y1), [l1, l2, 1 - l1 - l2], mp.mpf(1))
    s3 = P['sig'][2]; tc = P['t'] / s3; cz = mp.conj(tc) ** 2 * P['Z']
    return dict(z2=abs(P['Z']) ** 2, jac=s3 ** 2 * P['J'], Ich=P['Ich'], nu=P['nu'], czre=cz.real, czim=cz.imag, W=P['rho'][0] * P['t'].imag)


n = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
a0 = 2 * math.pi / 3
c = np.column_stack([rng.uniform(a0, math.pi, n), rng.uniform(-1, 1, n), rng.uniform(0, 1, n), rng.uniform(0, 1, n)])
sw = c[:, 2] + c[:, 3] > 1
c[sw, 2], c[sw, 3] = 1 - c[sw, 2], 1 - c[sw, 3]
edge = rng.integers(0, 6, n)
c[edge == 0, 0] = np.nextafter(math.pi, 4); c[edge == 1, 2] = 0.0; c[edge == 2, 3] = 0.0
w = 10.0 ** rng.uniform(-10, -2, (n, 4))
f = rng.uniform(0, 1, (n, 4))
blo = c - f * w; bhi = c + (1 - f) * w
blo[:, 0] = np.maximum(blo[:, 0], np.nextafter(a0, 0)); bhi[:, 0] = np.minimum(bhi[:, 0], np.nextafter(math.pi, 4))
blo[:, 1] = np.maximum(blo[:, 1], -1); bhi[:, 1] = np.minimum(bhi[:, 1], 1)
blo[:, 2:] = np.maximum(blo[:, 2:], 0); bhi[:, 2:] = np.minimum(bhi[:, 2:], 1)
bhi = np.maximum(bhi, blo)
with np.errstate(all='ignore'):
    q = ev.quantities(*(Iv(blo[:, j], bhi[:, j]) for j in range(4)))
fin = {k: int((np.isfinite(q[k].a) & np.isfinite(q[k].b)).sum()) for k in KEYS}
print('boxes', n, ' finite enclosures per quantity:', fin)
viol = {k: 0 for k in KEYS}; chk = {k: 0 for k in KEYS}
for i in range(n):
    for p in (0.5 * (blo[i] + bhi[i]), blo[i], bhi[i], blo[i] + (bhi[i] - blo[i]) * rng.uniform(size=4)):
        p = np.minimum(np.maximum(p, blo[i]), bhi[i])
        try:
            ex = exact(p)
        except ZeroDivisionError:
            continue
        for k in KEYS:
            lo_, hi_ = float(q[k].a[i]), float(q[k].b[i])
            if math.isfinite(lo_) and math.isfinite(hi_):
                chk[k] += 1
                if not (mp.mpf(lo_) <= ex[k] <= mp.mpf(hi_)):
                    viol[k] += 1
                    if sum(viol.values()) < 6:
                        print('   violation', k, i, p.tolist(), lo_, hi_, mp.nstr(ex[k], 20))
print('point checks per quantity:', chk)
print('violations per quantity  :', viol)
violations += sum(viol.values())
if violations:
    raise SystemExit(1)
print('all validation checks passed')
