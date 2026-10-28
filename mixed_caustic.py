"""Mixed signs q = (+,+,-).  Minima live in R_A = {lam1, lam2 > 0 > lam3} (beyond edge p1p2) where |Z| < kappa;
maxima are impossible.  Trace dM_A, map it to charge space, compute winding numbers (= number of minima),
and the projective invariants of the caustic (cusps, inflections) as in the positive case.
"""
import sys
import numpy as np
from scipy.optimize import minimize
from fast import tri, bary, Zfun


def in_RA(p, z):
    l = bary(p, z)
    return (l[..., 0] > 0) & (l[..., 1] > 0) & (l[..., 2] < 0)


def centers_RA(p, s, ngrid=400, R=8.0):
    """local minima of |Z| inside R_A (candidates for isotropic points)"""
    kappa = s / (s + 2)
    c = 0.5 * (p[0] + p[1])
    n_out = (p[0] + p[1]) / 2 - p[2]
    n_out /= abs(n_out)
    xs = np.linspace(-R, R, ngrid)
    X, Y = np.meshgrid(xs, xs)
    z = c + (X + 1j * Y).ravel()
    m = in_RA(p, z)
    z = z[m]
    a = np.abs(Zfun(p, z))
    idx = np.argsort(a)[:40]
    found = []
    for i in idx:
        f = lambda v: abs(Zfun(p, np.array([v[0] + 1j * v[1]]))[0]) if in_RA(p, np.array([v[0] + 1j * v[1]]))[0] else 9.0
        r = minimize(f, [z[i].real, z[i].imag], method='Nelder-Mead', options=dict(xatol=1e-12, fatol=1e-14, maxiter=4000))
        w = r.x[0] + 1j * r.x[1]
        if r.fun < kappa and all(abs(w - u) > 1e-4 for u, _ in found):
            found.append((w, r.fun))
    return found


def trace_RA(p, s, c, n=4096, Rmax=50.0):
    kappa = s / (s + 2)
    th = 2 * np.pi * np.arange(n) / n
    e = np.exp(1j * th)
    # march outwards until leaving R_A or |Z| >= kappa
    ts = np.concatenate([np.linspace(0, 1, 200)[1:], np.linspace(1, Rmax, 400)[1:]])
    P = c + np.outer(ts, e)  # (nt, n)
    inside = in_RA(p, P) & (np.abs(Zfun(p, P)) < kappa)
    first = np.argmax(~inside, axis=0)
    nonstar = np.count_nonzero(np.count_nonzero(inside[1:] != inside[:-1], axis=0) != 1)
    lo = np.where(first > 0, ts[np.maximum(first - 1, 0)], 0.0)
    hi = ts[first]
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        zz = c + mid * e
        ins = in_RA(p, zz) & (np.abs(Zfun(p, zz)) < kappa)
        lo = np.where(ins, mid, lo); hi = np.where(ins, hi, mid)
    return c + 0.5 * (lo + hi) * e, nonstar


def charges_mixed(p, z, s):
    lam = bary(p, z)
    r = np.abs(z[..., None] - p)
    q = lam * r ** (s + 2)
    return q / q[..., :2].sum(axis=-1, keepdims=True)  # class (+,+,-), normalised q1+q2 = 1


def winding(Y, P):
    d = Y - P
    ang = np.arctan2(d[:, 1], d[:, 0])
    da = np.diff(np.concatenate([ang, ang[:1]]))
    da = (da + np.pi) % (2 * np.pi) - np.pi
    return int(np.round(da.sum() / 2 / np.pi))


def analyse(p, s, n=4096, grid=120):
    cs = centers_RA(p, s)
    if not cs:
        return dict(ncomp=0, maxwind=0)
    curves = []
    nonstar = 0
    for c, _ in cs:
        zb, ns = trace_RA(p, s, c, n)
        nonstar += ns
        q = charges_mixed(p, zb, s)
        Y = np.stack([np.log(q[:, 0]), np.log(-q[:, 2])], -1)  # log q1, log|q3| with q1+q2=1
        curves.append(Y)
    allY = np.concatenate(curves)
    xs = np.linspace(allY[:, 0].min(), allY[:, 0].max(), grid)
    ys = np.linspace(allY[:, 1].min(), allY[:, 1].max(), grid)
    mw = 0
    for x in xs:
        for y in ys:
            w = sum(abs(winding(Y, np.array([x, y]))) for Y in curves)
            mw = max(mw, w)
    return dict(ncomp=len(cs), maxwind=mw, nonstar=nonstar)


if __name__ == '__main__':
    s = 1.0
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    mx = 0
    for k in range(N):
        a, b = rng.uniform(-1.0, 2.0), 10 ** rng.uniform(-1.5, 0.5)
        r = analyse(tri(a, b), s)
        mx = max(mx, r['maxwind'])
        print(f'a={a:+.3f} b={b:.3f}: {r}', flush=True)
    print('max total |winding| (= max #minima):', mx)
