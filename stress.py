"""Stress test over triangle shapes and exponents s:
  (1) number of connected components of M and its convexity,
  (2) self-intersections of the caustic F(dM) (=> possible 6 equilibria).
"""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from core import triangle, barycentric, charges_at, detT


def contours(a, b, s, n=900):
    P = triangle(a, b)
    xs = np.linspace(min(0, a), max(1, a), n)
    ys = np.linspace(0.0, b, n)
    Xg, Yg = np.meshgrid(xs, ys)
    X = np.stack([Xg, Yg], -1)
    lam = barycentric(P, X)
    D = detT(P, X)
    D[lam.min(-1) <= 0] = 0.0
    thr = (s + 1) / (s + 2) ** 2
    fig = plt.figure()
    cs = plt.contour(Xg, Yg, D, levels=[thr])
    segs = [sg for sg in cs.allsegs[0] if len(sg) > 5]
    plt.close(fig)
    return P, segs, D.max(), thr


def poly_area(p):
    x, y = p[:, 0], p[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def self_intersections(Y):
    """count proper crossings of a closed polyline (O(n^2) vectorised in chunks)"""
    A = Y
    B = np.roll(Y, -1, axis=0)
    n = len(A)
    cnt = 0
    for i in range(n):
        j = np.arange(i + 2, n)
        if i == 0:
            j = j[j != n - 1]
        if len(j) == 0:
            continue
        p, r = A[i], B[i] - A[i]
        qv, sv = A[j], B[j] - A[j]
        den = r[0] * sv[:, 1] - r[1] * sv[:, 0]
        ok = np.abs(den) > 1e-300
        t = ((qv[:, 0] - p[0]) * sv[:, 1] - (qv[:, 1] - p[1]) * sv[:, 0]) / np.where(ok, den, 1)
        u = ((qv[:, 0] - p[0]) * r[1] - (qv[:, 1] - p[1]) * r[0]) / np.where(ok, den, 1)
        cnt += int((ok & (t > 0) & (t < 1) & (u > 0) & (u < 1)).sum())
    return cnt


if __name__ == '__main__':
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    bad = 0
    for k in range(N):
        while True:
            a, b = rng.uniform(0, 0.5), 10 ** rng.uniform(-2.5, 0)
            if (a - 1) ** 2 + b ** 2 <= 1:
                break
        s = 10 ** rng.uniform(-1.3, 1.7)
        P, segs, Dmax, thr = contours(a, b, s)
        if not segs:
            print(f'a={a:.4f} b={b:.4f} s={s:.3f}: M empty (maxD={Dmax:.5f}, thr={thr:.5f})')
            continue
        info = []
        for sg in segs:
            closed = np.allclose(sg[0], sg[-1])
            hull = ConvexHull(sg)
            conv = poly_area(sg) / hull.volume
            q = charges_at(P, sg, s)
            Y = np.log(q[:, :2] / q[:, 2:3])
            si = self_intersections(Y[:-1] if closed else Y)
            info.append((closed, round(conv, 5), si))
        flag = len(segs) != 1 or any(i[1] < 0.999 or i[2] > 0 or not i[0] for i in info)
        bad += flag
        print(f'{"!!" if flag else "  "} a={a:.4f} b={b:.4f} s={s:.3f} comps={len(segs)} (closed,convexity,caustic-crossings)={info}')
    print('flagged:', bad)
