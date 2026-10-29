"""Theory-independent check: count planar equilibria of three positive Coulomb charges by
multistart damped Newton on grad V over the triangle, for random triangles (incl. very flat)
and random charges (log-uniform over 6 decades).  Prints the histogram of counts."""
import sys
import numpy as np
from multiprocessing import Pool


def count(args):
    seed, n, spread = args
    rng = np.random.default_rng(seed)
    hist = {}
    worst = None
    for _ in range(n):
        while True:
            a, b = rng.uniform(0, 0.5), 10 ** rng.uniform(-2.5, 0)
            if (a - 1) ** 2 + b ** 2 <= 1:
                break
        p = np.array([0, 1, a + 1j * b])
        q = 10 ** rng.uniform(-spread, spread, 3)
        # starts: barycentric grid
        m = 70
        us = (np.arange(m) + 0.5) / m
        A, B = np.meshgrid(us, us)
        k = A + B < 1
        z = A[k] * p[0] + B[k] * p[1] + (1 - A[k] - B[k]) * p[2]
        # multi-scale starts near each vertex (equilibria can hug a small charge)
        extra = []
        for i in range(3):
            j, l = (i + 1) % 3, (i + 2) % 3
            for eps in 10 ** np.linspace(-7, -0.5, 40):
                for t in np.linspace(0.02, 0.98, 25):
                    dirn = (1 - t) * (p[j] - p[i]) + t * (p[l] - p[i])
                    extra.append(p[i] + eps * dirn)
        z = np.concatenate([z, np.array(extra)])
        for _ in range(150):
            d = z[:, None] - p
            r = np.abs(d)
            w = q / r ** 3
            g = np.sum(w * d, 1)                       # = -grad V (complex)
            # Jacobian of (Re g, Im g): d/dz and d/dzbar of g
            u2 = d / d.conjugate()
            gz = np.sum(w * (1 - 1.5), 1)              # d g / dz      = sum w (1 - 3/2)
            gzb = np.sum(-1.5 * w * u2, 1)              # d g / d zbar  = -3/2 sum w u^2
            # solve gz*dz + gzb*conj(dz) = -g
            det = np.abs(gz) ** 2 - np.abs(gzb) ** 2
            det = np.where(np.abs(det) < 1e-300, 1e-300, det)
            dz = (np.conj(gz) * (-g) - gzb * np.conj(-g)) / det
            nr = np.abs(dz)
            lim = 0.25 * r.min(1)
            dz = np.where(nr > lim, dz * lim / np.maximum(nr, 1e-300), dz)
            z = z + dz
        d = z[:, None] - p
        r = np.abs(d)
        g = np.sum(q / r ** 3 * d, 1)
        ok = np.isfinite(z) & (np.abs(g) * r.min(1) ** 2 / q.max() < 1e-9) & (r.min(1) > 1e-9 * b)
        sols = []
        for w_ in z[ok]:
            if all(abs(w_ - u) > 1e-7 * b for u in sols):
                sols.append(w_)
        c = len(sols)
        hist[c] = hist.get(c, 0) + 1
        if c > 4 or c % 2:
            worst = (a, b, q, sols)
    return hist, worst


if __name__ == '__main__':
    nw = int(sys.argv[1]) if len(sys.argv) > 1 else 16
    per = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
    spread = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0
    with Pool(nw) as pool:
        res = pool.map(count, [(s, per, spread) for s in range(nw)])
    tot = {}
    for h, w in res:
        for k, v in h.items():
            tot[k] = tot.get(k, 0) + v
        if w is not None:
            print('ANOMALY', w)
    print('histogram of equilibrium counts:', dict(sorted(tot.items())))
