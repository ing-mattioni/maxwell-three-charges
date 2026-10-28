"""Vectorised tracing of dM = {|Z| = s/(s+2)} and its dual curve.

Z(x) = sum_i lam_i(x) u_i(x)^2  (complex),  M = {|Z| < kappa},  kappa = s/(s+2).
"""
import numpy as np
from scipy.optimize import minimize


def tri(a, b):
    return np.array([0.0, 1.0, a + 1j * b])


def bary(p, z):
    p1, p2, p3 = p
    den = ((p2 - p3).conjugate() * (p1 - p3)).imag
    l1 = ((p2 - p3).conjugate() * (z - p3)).imag / den
    l2 = ((z - p3).conjugate() * (p1 - p3)).imag / den
    return np.stack([l1, l2, 1 - l1 - l2], axis=-1)


def Zfun(p, z):
    lam = bary(p, z)
    d = z[..., None] - p
    u2 = d / d.conjugate()
    return np.sum(lam * u2, axis=-1)


def center(p):
    """the isotropic point Z = 0 (minimiser of |Z|)"""
    f = lambda v: abs(Zfun(p, np.array([v[0] + 1j * v[1]]))[0]) ** 2 if bary(p, np.array([v[0] + 1j * v[1]])).min() > 0 else 10.0
    g = p.mean()
    r = minimize(f, [g.real, g.imag], method='Nelder-Mead', options=dict(xatol=1e-14, fatol=1e-30, maxiter=20000))
    return r.x[0] + 1j * r.x[1], np.sqrt(r.fun)


def ray_exit(p, c, e):
    """distance from c along unit complex e to triangle boundary"""
    lam_c = bary(p, np.array([c]))[0]
    lam_e = bary(p, c + e) - lam_c  # linear part
    t = np.full(e.shape, np.inf)
    for i in range(3):
        li = lam_e[..., i]
        with np.errstate(divide='ignore'):
            ti = np.where(li < 0, -lam_c[i] / li, np.inf)
        t = np.minimum(t, ti)
    return t


def trace(p, s, n=4096, c=None, iters=70):
    kappa = s / (s + 2)
    if c is None:
        c, _ = center(p)
    th = 2 * np.pi * np.arange(n) / n
    e = np.exp(1j * th)
    R = ray_exit(p, c, e)
    lo = np.zeros(n)
    hi = R * (1 - 1e-13)
    # sanity: |Z| increasing along ray? find first crossing by coarse scan
    ts = np.linspace(0, 1, 257)[1:]
    Zs = np.abs(Zfun(p, c + np.outer(hi, ts) * e[:, None]))  # (n, 256)
    above = Zs >= kappa
    multi = np.count_nonzero(above[:, 1:] != above[:, :-1], axis=1)
    first = np.argmax(above, axis=1)
    lo = np.where(first > 0, hi * ts[np.maximum(first - 1, 0)], 0.0)
    hi2 = hi * ts[first]
    for _ in range(iters):
        mid = 0.5 * (lo + hi2)
        inside = np.abs(Zfun(p, c + mid * e)) < kappa
        lo = np.where(inside, mid, lo)
        hi2 = np.where(inside, hi2, mid)
    return c, c + 0.5 * (lo + hi2) * e, int(np.count_nonzero(multi != 1))


def dual(p, s, zb):
    """kernel direction k and dual vector m along boundary points zb"""
    Z = Zfun(p, zb)
    k = np.exp(1j * (np.angle(Z) / 2 + np.pi / 2))
    d = zb[:, None] - p
    r = np.abs(d)
    u = d / r
    m = r ** (-s - 1) * (u * k[:, None].conjugate()).real
    return k, m


def charges(p, z, s):
    lam = bary(p, z)
    r = np.abs(z[..., None] - p)
    q = lam * r ** (s + 2)
    return q / q.sum(axis=-1, keepdims=True)


def invariants(p, s, n=4096):
    c, zb, nonstar = trace(p, s, n)
    k, m = dual(p, s, zb)
    qc = charges(p, np.array([c]), s)[0]
    m = m * qc
    m /= np.linalg.norm(m, axis=1, keepdims=True)
    flip = np.ones(n)
    for i in range(1, n):
        if np.dot(m[i], m[i - 1]) < 0:
            m[i] = -m[i]
            k[i] = -k[i]
    one_sided = np.dot(m[-1], m[0]) < 0
    sg = -1.0 if one_sided else 1.0
    ext = np.vstack([sg * m[-2:], m, sg * m[:2]])
    d1 = (ext[3:-1] - ext[1:-3]) / 2
    d2 = ext[3:-1] - 2 * ext[2:-2] + ext[1:-3]
    tv = d1 - np.sum(d1 * m, axis=1, keepdims=True) * m
    dots = np.einsum('ij,ij->i', tv[1:], tv[:-1])
    infl = int(np.count_nonzero(dots < 0)) + int(np.dot(tv[-1], sg * tv[0]) < 0)
    det = np.einsum('ij,ij->i', m, np.cross(d1, d2))
    cusps = int(np.count_nonzero(np.sign(det[1:]) * np.sign(det[:-1]) < 0)) + int(np.sign(det[-1]) * np.sign(sg * det[0]) < 0)
    speed = np.linalg.norm(tv, axis=1)
    return dict(one_sided=bool(one_sided), inflections=infl, cusps=cusps, nonstar_rays=nonstar,
                min_rel_speed=float(speed.min() / np.median(speed)), c=c, zb=zb, m=m, k=k)
