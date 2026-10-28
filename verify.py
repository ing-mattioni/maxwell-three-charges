"""Computer-assisted verification of the three boundary conditions for s = 1 (Coulomb).

Domain (WLOG by cyclic relabelling, theta3 is the largest angle):
    theta1, theta2 in [0, pi],  theta3 = 2pi - theta1 - theta2 in [max(theta1, theta2), pi],
    lambda1, lambda2 >= 0, lambda3 = 1 - lambda1 - lambda2 >= 0.
Conditions to certify:
  (J)  jac  < 0          on  { |zeta| <= kappa }         (Z is an orientation-reversing local diffeo on M-bar)
  (I)  Ih   > 0          on  { |zeta|  = kappa }         (caustic has no inflections)
  (T)  turn < 1/2        on  { |zeta|  = kappa }         (exactly three simple cusps)
where zeta = sum lambda_i v_i^2 and kappa = s/(s+2).
"""
import sys
import time
import numpy as np
from interval import I, C, sin_i, cos_i, PI_LO, PI_HI

TWO_PI = I(2 * PI_LO, 2 * PI_HI)


def thetas_from_y(y1, y2):
    """rectangle coordinates: y1 = theta3 in [2pi/3, pi], y2 = tau in [-1, 1];
    theta1,2 = (2pi - theta3)/2 +- tau (3 theta3 - 2pi)/2   (so theta3 = max, sum = 2pi exactly)"""
    base = (TWO_PI - y1) * 0.5
    half = (3.0 * y1 - TWO_PI) * 0.5
    phi = y2 * half
    return base + phi, base - phi, y1


def enclose_y(y1, y2, l1, l2, s, need_bdry=True):
    th1, th2, th3 = thetas_from_y(y1, y2)
    return enclose(th1, th2, l1, l2, s, need_bdry, th3=th3)


def enclose(th1, th2, l1, l2, s, need_bdry=True, th3=None):
    if th3 is None:
        th3 = TWO_PI - th1 - th2
    l3 = 1.0 - l1 - l2
    s1, s2, s3 = sin_i(th1), sin_i(th2), sin_i(th3)
    c1, c2, c3 = cos_i(th1), cos_i(th2), cos_i(th3)
    n = th1.lo.shape[-1]
    zero = I(np.zeros(n))
    v = [C(I(np.ones(n)), zero), C(c3, s3), C(c2, -s2)]
    lam = [l1, l2, l3]
    rho = [l1 / s1, l2 / s2, None]
    rhoh = [s3 * rho[0], s3 * rho[1], l3]
    v2 = [vi.sq() for vi in v]
    zeta = v2[0] * l1 + v2[1] * l2 + v2[2] * l3
    zeta2 = zeta.abs2()

    def s3rr(a, b):
        if a == 2:
            return l3 * rho[b]
        if b == 2:
            return l3 * rho[a]
        return s3 * rho[a] * rho[b]

    muh = []
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        w = v[j].conj() * s3rr(i, k) - v[k].conj() * s3rr(i, j)
        muh.append(C(-0.5 * w.im, 0.5 * w.re))  # (i/2) * w
    Zz = None
    Zzb = None
    for i in range(3):
        a = muh[i] * v2[i] - v[i] * (lam[i] * rhoh[i])
        b = muh[i].conj() * v2[i] + (v2[i] * v[i]) * (lam[i] * rhoh[i])
        Zz = a if Zz is None else Zz + a
        Zzb = b if Zzb is None else Zzb + b
    jac = Zz.abs2() - Zzb.abs2()
    out = dict(zeta2=zeta2, jac=jac)
    if not need_bdry:
        return out
    Z = zeta
    iZ = C(-Z.im, Z.re)
    num = Zz.conj() * iZ - Zzb * iZ.conj()
    t = num / jac
    absz = zeta2.sqrt()
    Kb = C(-Z.re, Z.im) / absz          # -conj(zeta)/|zeta|
    dot = lambda a, b: a.re * b.re + a.im * b.im

    def Q(a, b):
        ab = a * b * Kb
        return 0.5 * dot(a, b) + 0.5 * ab.re

    v1, v3 = v[0], v[2]
    Q13 = Q(v1, v3)
    Ih = 0.5 * (s1 * s2) + (s + 2) * Q13 * (l1 * s3 * dot(v1, t) - s1 * l3 * dot(v3, t)) \
        - l1 * s3 * Q(v3, t) + s1 * l3 * Q(v1, t)
    Zzz = None; Zzzb = None; Zzbzb = None
    for i in range(3):
        a = (muh[i] * v[i]) * (-2.0 * rhoh[i])
        b = (muh[i] * (v2[i] * v[i])) * rhoh[i] - (muh[i].conj() * v[i]) * rhoh[i] - v2[i] * (lam[i] * rhoh[i].sq())
        c = (muh[i].conj() * (v2[i] * v[i])) * (2.0 * rhoh[i]) + v2[i].sq() * (2.0 * lam[i] * rhoh[i].sq())
        Zzz = a if Zzz is None else Zzz + a
        Zzzb = b if Zzzb is None else Zzzb + b
        Zzbzb = c if Zzbzb is None else Zzbzb + c
    rhs = -Z - (Zzz * t.sq() + Zzzb * (2.0 * t.abs2()) + Zzbzb * t.conj().sq())
    tp = (Zz.conj() * rhs - Zzb * rhs.conj()) / jac
    tpt = tp * t.conj()
    turn = tpt.im / t.abs2()
    # cusp locus: t parallel to k  <=>  conj(t)^2 * zeta is real and negative
    cz = t.conj().sq() * Z
    # (W): |d arg u_1 / d omega| = sigma3 * |Im(th * conj(v1))| * rho_1, with v1 = 1
    Wq = s3 * t.im * rho[0]
    out.update(Ih=Ih, turn=turn, czre=cz.re, czim=cz.im, Wq=Wq)
    return out


KEYS = ('zeta2', 'jac', 'Ih', 'turn', 'czre', 'czim')


def enclose_mvf(blo, bhi, s):
    """mean-value-form enclosures (intersected with naive ones) of all quantities on boxes"""
    from ad import DI, mvf
    c = 0.5 * (blo + bhi)
    ec = enclose(*(I(c[:, j], c[:, j]) for j in range(4)), s)
    eb = enclose(*(DI.var(I(blo[:, j], bhi[:, j]), j) for j in range(4)), s)
    return {k: mvf(eb[k], ec[k], blo, bhi) for k in KEYS}


def run(s=1.0, maxit=60, batch=200000, verbose=True, tag=''):
    kappa = s / (s + 2)
    k2lo, k2hi = kappa * kappa * (1 - 1e-15), kappa * kappa * (1 + 1e-15)
    # initial boxes: theta in [1.0, pi] (smaller theta are excluded analytically below), lambda in [0,1]^2
    lo = np.array([[1.0, 1.0, 0.0, 0.0]])
    hi = np.array([[PI_HI, PI_HI, 1.0, 1.0]])
    # subdivide initially
    n0 = 16
    g = np.linspace(0, 1, n0 + 1)
    L, H = [], []
    for a in range(n0):
        for b in range(n0):
            for c in range(n0):
                for d in range(n0):
                    L.append([1.0 + (PI_HI - 1.0) * g[a], 1.0 + (PI_HI - 1.0) * g[b], g[c], g[d]])
                    H.append([1.0 + (PI_HI - 1.0) * g[a + 1], 1.0 + (PI_HI - 1.0) * g[b + 1], g[c + 1], g[d + 1]])
    lo, hi = np.array(L), np.array(H)
    stats = dict(verified_bdry=0, verified_int=0, discarded=0, boxes=0)
    vol_ver = 0.0
    it = 0
    failed = []
    t0 = time.time()
    while len(lo) and it < maxit:
        it += 1
        newlo, newhi = [], []
        for k in range(0, len(lo), batch):
            blo, bhi = lo[k:k + batch], hi[k:k + batch]
            stats['boxes'] += len(blo)
            th1, th2, l1, l2 = (I(blo[:, j], bhi[:, j]) for j in range(4))
            th3 = TWO_PI - th1 - th2
            l3 = 1.0 - l1 - l2
            # structural discards
            disc = (th3.hi < th1.lo) | (th3.hi < th2.lo) | (th3.lo > PI_HI) | (l3.hi < 0)
            keep = ~disc
            stats['discarded'] += int(disc.sum())
            blo, bhi = blo[keep], bhi[keep]
            if len(blo) == 0:
                continue
            th1, th2, l1, l2 = (I(blo[:, j], bhi[:, j]) for j in range(4))
            e = enclose_mvf(blo, bhi, s)
            z2 = e['zeta2']
            outside = z2.lo > k2hi            # |zeta| > kappa: not in M-bar
            inside = z2.hi < k2lo             # strictly inside M: only (J) needed
            ok_int = inside & (e['jac'].hi < 0)
            ok_bd = (~outside) & (~inside) & (e['jac'].hi < 0) & (e['Ih'].lo > 0) & ((e['turn'].hi < 0.5) | ~((e['czim'].lo <= 0) & (e['czim'].hi >= 0) & (e['czre'].lo <= 0)))
            done = outside | ok_int | ok_bd
            stats['discarded'] += int(outside.sum())
            stats['verified_int'] += int(ok_int.sum())
            stats['verified_bdry'] += int(ok_bd.sum())
            blo, bhi = blo[~done], bhi[~done]
            if len(blo) == 0:
                continue
            # bisect along the dimension with largest width (scaled)
            w = (bhi - blo) * np.array([1.0, 1.0, 2.0, 2.0])
            dim = np.argmax(w, axis=1)
            mid = 0.5 * (blo[np.arange(len(blo)), dim] + bhi[np.arange(len(blo)), dim])
            b1hi = bhi.copy(); b1hi[np.arange(len(blo)), dim] = mid
            b2lo = blo.copy(); b2lo[np.arange(len(blo)), dim] = mid
            newlo += [blo, b2lo]
            newhi += [b1hi, bhi]
        lo = np.concatenate(newlo) if newlo else np.zeros((0, 4))
        hi = np.concatenate(newhi) if newhi else np.zeros((0, 4))
        if verbose:
            print(f'{tag} it={it:3d} pending={len(lo):9d} stats={stats} min_width={(hi-lo).min(axis=0) if len(lo) else None} t={time.time()-t0:.0f}s', flush=True)
    return lo, hi, stats


if __name__ == '__main__':
    s = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
    maxit = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    lo, hi, st = run(s, maxit)
    np.save('pending_lo.npy', lo); np.save('pending_hi.npy', hi)
    print('FINAL', st, 'pending', len(lo))
