"""Test property P3: for every x in M, x minimises V_{q(x)} over M  (S|M is a convex cap).
Also P2: central projection of S|M,  x -> ((r3/r1)^s, (r3/r2)^s), has convex image.
"""
import sys
import numpy as np
from scipy.spatial import ConvexHull
from core import triangle, charges_at, detT
from scan_injectivity import grid_in_triangle


def test(a, b, s, n=220, nsamp=300, rng=None):
    P = triangle(a, b)
    X, lam = grid_in_triangle(P, n)
    inM = detT(P, X) > (s + 1) / (s + 2) ** 2
    Xm = X[inM]
    if len(Xm) < 10:
        return None
    r = np.linalg.norm(Xm[:, None, :] - P, axis=-1)
    Phi = r ** (-s)
    idx = rng.choice(len(Xm), size=min(nsamp, len(Xm)), replace=False)
    Q = charges_at(P, Xm[idx], s)
    Vall = Q @ Phi.T  # (nsamp, |M|)
    V0 = Vall[np.arange(len(idx)), idx]
    viol = (Vall.min(axis=1) - V0) / V0
    # P2 convexity of central projection image
    Y = np.stack([(r[:, 2] / r[:, 0]) ** s, (r[:, 2] / r[:, 1]) ** s], -1)
    hull = ConvexHull(Y)
    return viol.min(), len(Xm), Y, hull


if __name__ == '__main__':
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    s_list = [float(sys.argv[2])] if len(sys.argv) > 2 else None
    worst = 0
    for k in range(40):
        while True:
            a, b = rng.uniform(0, 0.5), 10 ** rng.uniform(-1.5, 0)
            if (a - 1) ** 2 + b ** 2 <= 1:
                break
        s = s_list[0] if s_list else 10 ** rng.uniform(-1, 1.3)
        res = test(a, b, s, rng=rng)
        if res is None:
            continue
        v, m, Y, hull = res
        worst = min(worst, v)
        print(f'a={a:.3f} b={b:.3f} s={s:.3f} |M|={m:6d}  min_y [V(y)-V(x)]/V(x) = {v:.3e}')
    print('worst', worst)
