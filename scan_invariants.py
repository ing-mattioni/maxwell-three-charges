"""Scan (shape, s) for the projective invariants of the caustic.
Shape space: p1=0, p2=1, p3=a+ib with 0<=a<=1/2, (a-1)^2+b^2<=1 (p1p2 longest side).
Counts: cusps (sign changes of det(m,m',m'') away from dual cusps), inflections (dual cusps).
"""
import sys
import numpy as np
from fast import tri, trace, dual, charges


def robust_invariants(p, s, n=8192):
    c, zb, nonstar = trace(p, s, n)
    k, m = dual(p, s, zb)
    qc = charges(p, np.array([c]), s)[0]
    m = m * qc
    m /= np.linalg.norm(m, axis=1, keepdims=True)
    # continuous orientation
    sgn = np.ones(n)
    dots = np.einsum('ij,ij->i', m[1:], m[:-1])
    flips = np.concatenate([[1.0], np.where(dots < 0, -1.0, 1.0)])
    m = m * np.cumprod(flips)[:, None]
    one_sided = np.dot(m[-1], m[0]) < 0
    sg = -1.0 if one_sided else 1.0
    ext = np.vstack([sg * m[-2:], m, sg * m[:2]])
    d1 = (ext[3:-1] - ext[1:-3]) / 2
    d2 = ext[3:-1] - 2 * ext[2:-2] + ext[1:-3]
    tv = d1 - np.sum(d1 * m, axis=1, keepdims=True) * m
    sp = np.linalg.norm(tv, axis=1)
    tvn = tv / sp[:, None]
    dd = np.einsum('ij,ij->i', tvn[1:], tvn[:-1])
    infl_idx = np.nonzero(dd < 0)[0]
    infl = len(infl_idx) + int(np.dot(tvn[-1], sg * tvn[0]) < 0)
    # geodesic curvature sign of dual curve: det(m, t, t')
    det = np.einsum('ij,ij->i', m, np.cross(d1, d2))
    sc = np.nonzero(np.sign(det[1:]) * np.sign(det[:-1]) < 0)[0]
    # discard det sign changes within 3 samples of a dual cusp
    bad = set()
    for i in infl_idx:
        bad.update(range(i - 3, i + 4))
    cusp_idx = [i for i in sc if i not in bad]
    cusps = len(cusp_idx) + int(np.sign(det[-1]) * np.sign(sg * det[0]) < 0)
    return dict(one_sided=bool(one_sided), cusps=cusps, inflections=infl, nonstar=nonstar,
                min_rel_speed=float(sp.min() / np.median(sp)))


def sample_shape(rng):
    while True:
        a, b = rng.uniform(0, 0.5), 10 ** rng.uniform(-2.5, 0)
        if (a - 1) ** 2 + b ** 2 <= 1:
            return a, b


if __name__ == '__main__':
    seed = int(sys.argv[1]); N = int(sys.argv[2])
    s_fixed = float(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3] != 'rand' else None
    rng = np.random.default_rng(seed)
    out = []
    for j in range(N):
        a, b = sample_shape(rng)
        s = s_fixed if s_fixed else 10 ** rng.uniform(-1, 1)
        try:
            r = robust_invariants(tri(a, b), s)
        except Exception as e:
            print('ERR', a, b, s, e); continue
        out.append((a, b, s, r['one_sided'], r['cusps'], r['inflections'], r['nonstar'], r['min_rel_speed']))
        if not (r['one_sided'] and r['cusps'] == 3 and r['inflections'] == 0):
            print(f'!! a={a:.5f} b={b:.5f} s={s:.4f} {r}')
    out = np.array(out, dtype=float)
    np.save(f'data_inv_{seed}_{sys.argv[3] if len(sys.argv) > 3 else "rand"}.npy', out)
    good = (out[:, 3] == 1) & (out[:, 4] == 3) & (out[:, 5] == 0)
    print(f'{N} samples, hedgehog-type (one-sided, 3 cusps, 0 inflections): {good.sum()}')
