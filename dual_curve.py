"""Chart-free (projective) analysis of the caustic via its dual curve.

C   = caustic = F(dM) in the projective plane of charges.
C*  = dual curve = { [m(x0)] : x0 in dM },  m_i = r_i^{-s-1} (u_i . k(x0)).
Projective dictionary:
  inflection of C  <=>  cusp of C*        (velocity of [m] vanishes and reverses)
  cusp of C        <=>  inflection of C*  (det(m, m', m'') changes sign)
  C one-sided      <=>  continuous lift of m along dM closes with m_end = -m_start
  n(q) = # tangent lines of C through q = # zeros of <q, m> along dM.
We rescale charges diagonally by c = q(x*) (a projective map) for conditioning.
"""
import sys
import numpy as np
from core import triangle, barycentric, charges_at, Tmat
from inspect_case import trace_boundary, argmax_detT


def dual_curve(P, s, npts=6000):
    c, B = trace_boundary(P, s, m=npts)
    T = Tmat(P, B)
    w, V = np.linalg.eigh(T)
    k = V[..., :, 0]
    d = B[:, None, :] - P
    r = np.linalg.norm(d, axis=-1)
    u = d / r[..., None]
    m = r ** (-s - 1) * np.einsum('nij,nj->ni', u, k)
    qc = charges_at(P, c[None], s)[0]
    m = m * qc  # projective rescaling
    m /= np.linalg.norm(m, axis=1, keepdims=True)
    sign = 1.0
    for i in range(1, len(m)):
        if np.dot(m[i], m[i - 1]) < 0:
            m[i] = -m[i]
            k[i] = -k[i]
    return B, k, m, qc


def invariants(B, k, m):
    n = len(m)
    one_sided = np.dot(m[-1], m[0]) < 0
    mm = m.copy()
    # closed periodic version for derivatives
    ext = np.vstack([m[-2:] * (-1 if one_sided else 1), m, m[:2] * (-1 if one_sided else 1)])
    d1 = (ext[3:-1] - ext[1:-3]) / 2
    d2 = ext[3:-1] - 2 * ext[2:-2] + ext[1:-3]
    # tangent of the spherical curve, projected
    tvec = d1 - np.sum(d1 * m, axis=1, keepdims=True) * m
    tn = np.linalg.norm(tvec, axis=1)
    # cusp of C* (inflection of C): tangent reverses direction
    flips = np.einsum('ij,ij->i', tvec[1:], tvec[:-1]) < 0
    infl_C = int(flips.sum()) + int(np.dot(tvec[-1], tvec[0] * (-1 if one_sided else 1)) < 0)
    det = np.einsum('ij,ij->i', m, np.cross(d1, d2))
    cusps_C = int(np.count_nonzero(np.sign(det[1:]) * np.sign(det[:-1]) < 0))
    # closing the loop: the lift continues with -m, and det(-m,-m',-m'') = -det(m,m',m'')
    if np.sign(det[-1]) * np.sign(det[0]) * (-1 if one_sided else 1) < 0:
        cusps_C += 1
    return dict(one_sided=bool(one_sided), inflections=infl_C, cusps=cusps_C,
                min_speed=float(tn.min() / tn.mean()))


def tangent_count(m, q):
    v = m @ q
    return int(np.count_nonzero(np.sign(v[1:]) * np.sign(v[:-1]) < 0) + (np.sign(v[-1]) * np.sign(-v[0] if np.dot(m[-1], m[0]) < 0 else v[0]) < 0))


if __name__ == '__main__':
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    fixed_s = float(sys.argv[3]) if len(sys.argv) > 3 else None
    nbad = 0
    for j in range(N):
        while True:
            a, b = rng.uniform(0, 0.5), 10 ** rng.uniform(-2, 0)
            if (a - 1) ** 2 + b ** 2 <= 1:
                break
        s = fixed_s if fixed_s else 10 ** rng.uniform(-1, 0.7)
        P = triangle(a, b)
        B, k, m, qc = dual_curve(P, s, npts=4000)
        inv = invariants(B, k, m)
        ok = inv['one_sided'] and inv['inflections'] == 0 and inv['cusps'] == 3
        nbad += not ok
        print(f"{'  ' if ok else '!!'} a={a:.4f} b={b:.4f} s={s:.3f} {inv}")
    print('bad', nbad)
