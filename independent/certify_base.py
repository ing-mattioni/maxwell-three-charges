"""Independent check of item 3 of Theorem 5.1 on the base triangle (0, 1, 1/2 + i/8).

Usage:  python independent/certify_base.py
Writes independent_base.json in the current directory.  Exit status 1 if a box is left
unresolved.  Only NumPy is required; nothing is imported from the main verifier.

Claim: at every point x of the open triangle with |Z(x)| <= 1/3, the angle subtended by
the side [0, 1] is strictly the largest of the three angles at x.
Boxes live in barycentric coordinates (lambda1, lambda2); the starting box is the whole unit
square.  A box is closed when one of these holds on the whole box:
    lambda3 < 0                          (outside the triangle);
    some lambda_i > 2/3                  (then |Z| >= 2 lambda_i - 1 > 1/3);
    |Z|^2 > 1/9;
    (w1.w2) r3 < (w2.w3) r1  and  (w1.w2) r3 < (w1.w3) r2,   with w_i = p_i - x, r_i = |w_i|,
which says cos(theta3) < cos(theta1) and cos(theta3) < cos(theta2).
"""
import sys
import os
import json
import time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import evaluator as ev
from evaluator import Iv
from certify_domain import provenance

P0 = [(0.0, 0.0), (1.0, 0.0), (0.5, 0.125)]     # exactly representable
MIN_W = 1e-12


def step(lo, hi):
    l1, l2 = Iv(lo[:, 0], hi[:, 0]), Iv(lo[:, 1], hi[:, 1])
    l3 = 1.0 - l1 - l2
    lam = [l1, l2, l3]
    outdom = l3.b < 0
    third = Iv(1.0) / 3.0
    nearv = np.zeros(len(lo), bool)
    for li in lam:
        nearv |= (li * 2.0 - 1.0).a > third.b
    x = (l1 * P0[0][0] + l2 * P0[1][0] + l3 * P0[2][0], l1 * P0[0][1] + l2 * P0[1][1] + l3 * P0[2][1])
    w = [(Iv(p[0]) - x[0], Iv(p[1]) - x[1]) for p in P0]
    r2 = [ev.cabs2(wi) for wi in w]
    Z = None
    for wi, ri, li in zip(w, r2, lam):
        term = ev.cscale(ev.csqr(wi), li / ri)          # lambda_i w_i^2 / |w_i|^2 = lambda_i u_i^2
        Z = term if Z is None else ev.cadd(Z, term)
    outside = ev.cabs2(Z).a > ev.NINTH.b
    r = [ri.sqrt() for ri in r2]
    d12, d23, d13 = ev.rdot(w[0], w[1]), ev.rdot(w[1], w[2]), ev.rdot(w[0], w[2])
    largest = ((d12 * r[2] - d23 * r[0]).b < 0) & ((d12 * r[2] - d13 * r[1]).b < 0)
    return outdom | nearv | outside | largest


if __name__ == '__main__':
    t0 = time.time()
    lo = np.array([[0.0, 0.0]]); hi = np.array([[1.0, 1.0]])
    total = 0; unresolved = 0; depth = 0
    while len(lo) and depth < 80:
        total += len(lo)
        with np.errstate(all='ignore'):
            done = step(lo, hi)
        lo, hi = lo[~done], hi[~done]
        w_ = hi - lo
        tiny = w_.max(axis=1) < MIN_W
        unresolved += int(tiny.sum())
        lo, hi, w_ = lo[~tiny], hi[~tiny], w_[~tiny]
        dim = np.argmax(w_, axis=1); ar = np.arange(len(lo))
        mid = 0.5 * (lo[ar, dim] + hi[ar, dim])
        h1 = hi.copy(); h1[ar, dim] = mid
        l2_ = lo.copy(); l2_[ar, dim] = mid
        lo, hi = np.concatenate([lo, l2_]), np.concatenate([h1, hi])
        depth += 1
    unresolved += len(lo)
    summary = dict(schema_version=1,
                   description='independent check that theta3 is strictly the largest angle on the closed minimum region of the base triangle',
                   triangle=['0', '1', '1/2 + i/8'], initial_boxes=1, boxes=int(total), unresolved=int(unresolved),
                   seconds=time.time() - t0, parameters=dict(min_width=MIN_W))
    summary.update(provenance(['evaluator.py', 'certify_base.py', 'certify_domain.py']))
    json.dump(summary, open('independent_base.json', 'w'), indent=1)
    print(json.dumps(summary, indent=1))
    if unresolved:
        raise SystemExit(1)
