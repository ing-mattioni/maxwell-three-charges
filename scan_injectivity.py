"""Scan triangle shapes: is the charge map x -> q(x) injective on the minimum region M?

A failure of injectivity = two distinct local minima for the same charges = 6 equilibria,
i.e. a counterexample to Maxwell's bound for three positive charges.
"""
import sys
import numpy as np
from scipy.spatial import cKDTree
from core import triangle, barycentric, charges_at, detT


def grid_in_triangle(P, n):
    us = (np.arange(n) + 0.5) / n
    A, B = np.meshgrid(us, us, indexing='ij')
    m = A + B < 1
    lam = np.stack([A[m], B[m], 1 - A[m] - B[m]], axis=-1)
    return lam @ P, lam


def analyse(a, b, s, n=500):
    P = triangle(a, b)
    X, lam = grid_in_triangle(P, n)
    D = detT(P, X)
    thr = (s + 1) / (s + 2) ** 2
    inM = D > thr
    if inM.sum() < 2:
        return dict(area=inM.mean(), pairs=0, maxsep=0.0)
    Xm = X[inM]
    q = charges_at(P, Xm, s)
    Y = np.log(q[:, :2] / q[:, 2:3])
    h = 1.0 / n
    tree = cKDTree(Y)
    # image spacing scale: typical |dY| for neighbouring domain points
    t2 = cKDTree(Xm)
    dd, ii = t2.query(Xm, k=2)
    dY = np.linalg.norm(Y - Y[ii[:, 1]], axis=1)
    eps = 2.0 * np.median(dY)
    pairs = tree.query_pairs(eps, output_type='ndarray')
    if len(pairs) == 0:
        return dict(area=inM.mean(), pairs=0, maxsep=0.0)
    sep = np.linalg.norm(Xm[pairs[:, 0]] - Xm[pairs[:, 1]], axis=1)
    return dict(area=inM.mean(), pairs=int((sep > 20 * h).sum()), maxsep=float(sep.max()))


if __name__ == '__main__':
    s = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    rng = np.random.default_rng(1)
    worst = []
    for k in range(400):
        while True:
            a, b = rng.uniform(0, 0.5), rng.uniform(0, 1)
            if (a - 1) ** 2 + b ** 2 <= 1:
                break
        r = analyse(a, b, s, n=300)
        worst.append((r['maxsep'], a, b, r['area'], r['pairs']))
    worst.sort(reverse=True)
    for w in worst[:10]:
        print('maxsep=%.4f a=%.4f b=%.4f areaM=%.4f farpairs=%d' % w)
