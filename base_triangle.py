"""Rigorous check for one base triangle Delta0 = (0, 1, 1/2 + i/8) (counter-clockwise):

    for every x in Delta0:   |Z(x)| > kappa   or   theta3(x) > max(theta1(x), theta2(x)).

So on dM(Delta0) the largest angle is always theta3, the certificate frame of every boundary point
uses the same labelling as Delta0, and condition (W) holds for the fixed vertex p1 along all of
dM(Delta0).  This pins the total turning of the caustic at Delta0 (see notes/theory.md, Sec. 6).

Boxes live in barycentric coordinates (lambda1, lambda2); x = sum lambda_i p_i.  Boxes with some
lambda_i > 2/3 are discarded analytically, since |Z| >= 2 lambda_i - 1 > 1/3 there.
"""
import sys
import json
import time
import numpy as np
from interval import I, C
from domains import simplex_boxes
from provenance import metadata

s = 1.0
kappa = s / (s + 2)
kappa_hi = float((I(s) / (I(s) + 2.0)).hi)
P = [0.0 + 0j, 1.0 + 0j, 0.5 + 0.125j]


def check(blo, bhi):
    l1, l2 = I(blo[:, 0], bhi[:, 0]), I(blo[:, 1], bhi[:, 1])
    l3 = 1.0 - l1 - l2
    lam = [l1, l2, l3]
    # analytic discard near vertices: |Z| >= 2 lambda_i - 1 > kappa
    near_vertex = np.zeros(len(blo), bool)
    for li in lam:
        near_vertex |= (2.0 * li - 1.0).lo > kappa_hi
    # x = sum lambda_i p_i
    x = C(l1 * P[0].real + l2 * P[1].real + l3 * P[2].real, l1 * P[0].imag + l2 * P[1].imag + l3 * P[2].imag)
    w = [C.c(complex(p)) - x for p in P]                 # p_i - x
    r2 = [wi.abs2() for wi in w]
    v2 = [wi.sq() / ri for wi, ri in zip(w, r2)]         # v_i^2 = w_i^2 / |w_i|^2
    Z = v2[0] * l1 + v2[1] * l2 + v2[2] * l3
    outside = Z.abs2().lo > kappa * kappa * (1 + 1e-12)
    r = [ri.sqrt() for ri in r2]
    dot = lambda a, b: a.re * b.re + a.im * b.im
    d12, d23, d13 = dot(w[0], w[1]), dot(w[1], w[2]), dot(w[0], w[2])
    # theta3 > theta1  <=>  cos th3 < cos th1  <=>  d12 r3 < d23 r1 ;  theta3 > theta2 <=> d12 r3 < d13 r2
    ok3 = ((d12 * r[2] - d23 * r[0]).hi < 0) & ((d12 * r[2] - d13 * r[1]).hi < 0)
    return near_vertex | outside | ok3


def main(n0=64, maxit=40):
    lo, hi = simplex_boxes(n0)
    initial = len(lo)
    t0 = time.time()
    total = 0
    for it in range(maxit):
        if len(lo) == 0:
            break
        total += len(lo)
        with np.errstate(all='ignore'):
            done = check(lo, hi)
        lo, hi = lo[~done], hi[~done]
        w = hi - lo
        dim = np.argmax(w, axis=1)
        ar = np.arange(len(lo))
        mid = 0.5 * (lo[ar, dim] + hi[ar, dim])
        h1 = hi.copy(); h1[ar, dim] = mid
        l2 = lo.copy(); l2[ar, dim] = mid
        lo, hi = np.concatenate([lo, l2]), np.concatenate([h1, hi])
        if len(lo) and (hi - lo).max() < 1e-12:
            break
    res = dict(triangle=[str(p) for p in P], boxes=int(total), unresolved=int(len(lo)), seconds=time.time() - t0)
    res['initial_boxes'] = initial
    res.update(metadata(['base_triangle.py', 'interval.py', 'domains.py'], n0=n0, maxit=maxit))
    print(json.dumps(res))
    json.dump(res, open('cert_base_triangle.json', 'w'), indent=1)
    if len(lo):
        raise SystemExit(1)


if __name__ == '__main__':
    main(*(int(a) for a in sys.argv[1:]))
