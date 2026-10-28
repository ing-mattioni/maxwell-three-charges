"""Precise tracing of dM by bisection along rays from argmax det T, caustic, winding numbers."""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from core import triangle, barycentric, charges_at, detT, critical_points


def argmax_detT(P):
    f = lambda z: -detT(P, z[None])[0] if barycentric(P, z[None]).min() > 0 else 1.0
    best = None
    for lam0 in [(1/3, 1/3, 1/3), (0.6, 0.2, 0.2), (0.2, 0.6, 0.2), (0.2, 0.2, 0.6)]:
        r = minimize(f, np.array(lam0) @ P, method='Nelder-Mead', options=dict(xatol=1e-13, fatol=1e-15, maxiter=5000))
        if best is None or r.fun < best.fun:
            best = r
    return best.x, -best.fun


def trace_boundary(P, s, m=4000):
    c, Dmax = argmax_detT(P)
    thr = (s + 1) / (s + 2) ** 2
    ang = np.linspace(0, 2 * np.pi, m, endpoint=False)
    pts = []
    for t in ang:
        d = np.array([np.cos(t), np.sin(t)])
        # distance to triangle boundary along ray
        hi = 1e-9
        while barycentric(P, (c + 2 * hi * d)[None]).min() > 0 and hi < 10:
            hi *= 2
        lo_, hi_ = hi, 2 * hi
        for _ in range(60):
            mid = 0.5 * (lo_ + hi_)
            if barycentric(P, (c + mid * d)[None]).min() > 0:
                lo_ = mid
            else:
                hi_ = mid
        R = lo_ * (1 - 1e-12)
        # check monotonicity of D along ray coarse (star-shapedness)
        rr = np.linspace(0, R, 400)[1:]
        Dv = detT(P, c + rr[:, None] * d)
        sgn = Dv > thr
        changes = np.count_nonzero(sgn[1:] != sgn[:-1])
        if changes != 1:
            print('ray', t, 'has', changes, 'crossings of threshold')
        k = np.argmax(~sgn)
        a_, b_ = rr[k - 1] if k > 0 else 0.0, rr[k]
        for _ in range(60):
            mid = 0.5 * (a_ + b_)
            if detT(P, (c + mid * d)[None])[0] > thr:
                a_ = mid
            else:
                b_ = mid
        pts.append(c + 0.5 * (a_ + b_) * d)
    return c, np.array(pts)


def winding(Y, q):
    d = Y - q
    ang = np.arctan2(d[:, 1], d[:, 0])
    da = np.diff(np.concatenate([ang, ang[:1]]))
    da = (da + np.pi) % (2 * np.pi) - np.pi
    return np.round(da.sum() / (2 * np.pi)).astype(int)


if __name__ == '__main__':
    a, b, s = map(float, sys.argv[1:4])
    P = triangle(a, b)
    c, B = trace_boundary(P, s)
    q = charges_at(P, B, s)
    Y = np.log(q[:, :2] / q[:, 2:3])
    from stress import self_intersections
    print('caustic self-crossings (precise):', self_intersections(Y))
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    ax[0].plot(*np.vstack([P, P[:1]]).T, 'k-'); ax[0].plot(B[:, 0], B[:, 1], '.-', ms=1); ax[0].set_aspect('equal')
    ax[1].plot(Y[:, 0], Y[:, 1], '.-', ms=1)
    plt.savefig(sys.argv[4] if len(sys.argv) > 4 else 'figs/inspect.png', dpi=90)
