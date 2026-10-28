"""Collinear charges q1 at -1, q3 < 0 at c (|c| < 1), q2 at +1 (Coulomb, s = 1).
Count on-axis equilibria (zeros of E_x on the axis) and rings (zeros of (E_x, E_rho) with rho > 0).
After a small non-collinear perturbation each nondegenerate axis point persists and each ring leaves
exactly 2 in-plane equilibria, so  N_pert = n_axis + 2 n_rings."""
import sys
import numpy as np
from scipy.optimize import brentq, fsolve

s = 1.0


def Ex_axis(x, q, p):
    d = x - p
    return np.sum(q * np.sign(d) / np.abs(d) ** (s + 1))


def axis_zeros(q, p, n=20001):
    zs = []
    pts = np.sort(p)
    edges = [-60.0] + list(pts) + [60.0]
    for a, b in zip(edges[:-1], edges[1:]):
        xs = np.linspace(a, b, n)[1:-1]
        f = np.array([Ex_axis(x, q, p) for x in xs])
        for i in np.nonzero(np.sign(f[1:]) != np.sign(f[:-1]))[0]:
            zs.append(brentq(Ex_axis, xs[i], xs[i + 1], args=(q, p)))
    return zs


def ring_zeros(q, p, R=12.0, n=300):
    """solve S0 = sum q_i/r_i^(s+2) = 0 and S1 = sum q_i p_i / r_i^(s+2) = 0 for rho > 0"""
    def F(v):
        x, rho = v
        r = np.sqrt((x - p) ** 2 + rho ** 2)
        w = q / r ** (s + 2)
        return [np.sum(w) * 1.0, np.sum(w * p)]
    xs = np.linspace(-R, R, n); rs = np.linspace(1e-3, R, n)
    X, Rr = np.meshgrid(xs, rs)
    r = np.sqrt((X[..., None] - p) ** 2 + Rr[..., None] ** 2)
    w = q / r ** (s + 2)
    S0 = w.sum(-1); S1 = (w * p).sum(-1)
    found = []
    # seeds where both change sign nearby
    cand = np.nonzero((np.sign(S0[1:, 1:]) != np.sign(S0[:-1, :-1])) & (np.sign(S1[1:, 1:]) != np.sign(S1[:-1, :-1])))
    for i, j in zip(*cand):
        sol, info, ier, msg = fsolve(F, [X[i, j], Rr[i, j]], full_output=True)
        if ier == 1 and sol[1] > 1e-6 and np.max(np.abs(F(sol))) < 1e-10:
            if all(np.hypot(sol[0] - u[0], sol[1] - u[1]) > 1e-6 for u in found):
                found.append(sol)
    return found


if __name__ == '__main__':
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    best = 0
    for k in range(int(sys.argv[2]) if len(sys.argv) > 2 else 300):
        c = rng.uniform(-0.95, 0.95)
        p = np.array([-1.0, 1.0, c])
        q = np.array([1.0, 10 ** rng.uniform(-2, 2), -10 ** rng.uniform(-2, 2)])
        na = len(axis_zeros(q, p)); nr = len(ring_zeros(q, p))
        tot = na + 2 * nr
        if tot > best or tot > 4:
            best = max(best, tot)
            print(f'q={q}, c={c:.3f}: axis={na}, rings={nr}  ->  perturbed count {tot}', flush=True)
    print('max perturbed count', best)
