"""Parallel branch-and-bound certification of (J), (I), (T), (W) for exponent s.

Usage:  python certify.py s nworkers [maxit] [n0]
Writes certificate summary to cert_s{s}.json and any unresolved boxes to unresolved_s{s}.npy.
"""
import sys
import json
import time
import numpy as np
from multiprocessing import Pool
from interval import I, PI_HI
from verify import enclose_y as enclose, thetas_from_y
from ad import DI
from domains import simplex_boxes
from provenance import metadata

KEYS = ('zeta2', 'jac', 'Ih', 'turn', 'czre', 'czim', 'Wq')


def zeta2_naive(blo, bhi):
    from interval import sin_i, cos_i, C, TRIG_ARG_LO, TRIG_ARG_HI
    y1, y2, l1, l2 = (I(blo[:, j], bhi[:, j]) for j in range(4))
    th1, th2, th3 = thetas_from_y(y1, y2)
    l3 = 1.0 - l1 - l2
    # v_i^2 = e^{2 i beta_i}: v1^2 = 1, v2^2 = e^{2 i th3}, v3^2 = e^{-2 i th2}
    a2 = 2.0 * th3; a3 = 2.0 * th2
    def cs(x):
        # cos/sin for intervals possibly outside [0, 2pi]: reduce by multiples of 2pi (exact enough via interval)
        k = np.floor(x.lo / (2 * 3.141592653589793))
        xr = x - I(3.141592653589793, PI_HI) * (2.0 * k)
        # The reduction puts the lower endpoint in [0, 2pi), but a box that is wide in the
        # angles can still reach beyond the range where sin_i/cos_i are valid.  Use the
        # trivial enclosure [-1, 1] there.  The comparison is False for NaN, so NaN stays NaN.
        wide = (xr.lo < TRIG_ARG_LO) | (xr.hi > TRIG_ARG_HI)
        xs = I(np.where(wide, 0.0, xr.lo), np.where(wide, 0.0, xr.hi))
        c, s = cos_i(xs), sin_i(xs)
        return (I(np.where(wide, -1.0, c.lo), np.where(wide, 1.0, c.hi)),
                I(np.where(wide, -1.0, s.lo), np.where(wide, 1.0, s.hi)))
    c2, s2 = cs(a2); c3, s3 = cs(a3)
    re = l1 + l2 * c2 + l3 * c3
    im = l2 * s2 - l3 * s3
    return re.sq() + im.sq()


def enclose_mvf(blo, bhi, s, keys=KEYS, need_bdry=True):
    c = 0.5 * (blo + bhi)
    ec = enclose(*(I(c[:, j], c[:, j]) for j in range(4)), s, need_bdry=need_bdry)
    eb = enclose(*(DI.var(I(blo[:, j], bhi[:, j]), j) for j in range(4)), s, need_bdry=need_bdry)
    out, grad = {}, {}
    half = 0.5 * (bhi - blo)
    # rigorous enclosure of B - c (c is the float midpoint actually used above)
    dev = [I(blo[:, j], bhi[:, j]) - I(c[:, j], c[:, j]) for j in range(4)]
    for k in keys:
        acc = ec[k]
        contrib = np.zeros_like(blo)
        for j in range(4):
            dj = eb[k].d(j)
            acc = acc + dj * dev[j]
            contrib[:, j] = np.maximum(np.abs(dj.lo), np.abs(dj.hi)) * half[:, j]
        nv = eb[k].val()
        lower, upper = np.maximum(acc.lo, nv.lo), np.minimum(acc.hi, nv.hi)
        if np.any(lower > upper):
            raise ArithmeticError('Empty interval intersection for ' + k)
        out[k] = I(lower, upper)
        grad[k] = np.nan_to_num(contrib, nan=1e300, posinf=1e300)
    return out, grad


def process(args):
    blo, bhi, s, maxit, min_width = args
    kappa = s / (s + 2)
    k2lo, k2hi = kappa * kappa * (1 - 1e-12), kappa * kappa * (1 + 1e-12)
    st = dict(boxes=0, disc_struct=0, disc_outside=0, ver_int=0, ver_bdry=0, ver_bdry_nocusp=0)
    worst = dict(jac=-np.inf, Ih=np.inf, turn=-np.inf, W=-np.inf)
    unresolved = []
    for it in range(maxit):
        if len(blo) == 0:
            break
        st['boxes'] += len(blo)
        l3 = 1.0 - I(blo[:, 2], bhi[:, 2]) - I(blo[:, 3], bhi[:, 3])
        disc = (l3.hi < 0)
        st['disc_struct'] += int(disc.sum())
        blo, bhi = blo[~disc], bhi[~disc]
        if len(blo) == 0:
            break
        with np.errstate(all='ignore'):
            z2n = zeta2_naive(blo, bhi)
        if np.any(z2n.lo > z2n.hi):
            raise ArithmeticError('Invalid naive zeta2 interval')
        far = z2n.lo > k2hi
        st['disc_outside'] += int(far.sum())
        blo, bhi = blo[~far], bhi[~far]
        if len(blo) == 0:
            break
        with np.errstate(all='ignore'):
            e, g = enclose_mvf(blo, bhi, s)
        z2 = e['zeta2']
        outside = z2.lo > k2hi
        inside = z2.hi < k2lo
        jac_ok = e['jac'].hi < 0
        ih_ok = e['Ih'].lo > 0
        # NaN-safe: every 'ok' must be a positively established strict inequality
        no_cusp = (e['czim'].lo > 0) | (e['czim'].hi < 0) | (e['czre'].lo > 0)
        cusp_poss = ~no_cusp
        turn_ok = (e['turn'].hi < 0.5) | no_cusp
        ok_int = inside & jac_ok
        bd = ~outside & ~inside
        w_ok = (e['Wq'].lo > -0.5) & (e['Wq'].hi < 0.5)
        ok_bd = bd & jac_ok & ih_ok & turn_ok & w_ok
        st['disc_outside'] += int(outside.sum())
        st['ver_int'] += int(ok_int.sum())
        st['ver_bdry'] += int(ok_bd.sum())
        st['ver_bdry_nocusp'] += int((ok_bd & ~cusp_poss).sum())
        if ok_bd.any():
            worst['Ih'] = min(worst['Ih'], float(e['Ih'].lo[ok_bd].min()))
            tsel = ok_bd & cusp_poss
            if tsel.any():
                worst['turn'] = max(worst['turn'], float(e['turn'].hi[tsel].max()))
        if ok_bd.any():
            worst['W'] = max(worst['W'], float(np.maximum(-e['Wq'].lo, e['Wq'].hi)[ok_bd].max()))
        okj = ok_int | ok_bd
        if okj.any():
            worst['jac'] = max(worst['jac'], float(e['jac'].hi[okj].max()))
        done = outside | ok_int | ok_bd
        keep = ~done
        blo, bhi = blo[keep], bhi[keep]
        if len(blo) == 0:
            break
        # choose split dimension: largest contribution to the first failing quantity
        sc = np.zeros_like(blo)
        fj = ~jac_ok[keep]; fi = (~ih_ok[keep]) & bd[keep]; ft = (~turn_ok[keep]) & bd[keep]
        fw = (~w_ok[keep]) & bd[keep]
        fz = ~(fj | fi | ft | fw)
        sc[fw] = g['Wq'][keep][fw]
        sc[fz] = g['zeta2'][keep][fz]
        sc[ft] = g['turn'][keep][ft]
        sc[fi] = g['Ih'][keep][fi]
        sc[fj] = g['jac'][keep][fj]
        w = bhi - blo
        nonfinite = ~(sc < 1e299).all(axis=1)
        sc[nonfinite] = (w * np.array([1.0, 0.5, 2.0, 2.0]))[nonfinite]
        sc = np.where(w > min_width, sc, -1.0)
        dim = np.argmax(sc, axis=1)
        tiny = sc.max(axis=1) < 0
        if tiny.any():
            unresolved.append(np.concatenate([blo[tiny], bhi[tiny]], 1))
            blo, bhi, dim = blo[~tiny], bhi[~tiny], dim[~tiny]
        ar = np.arange(len(blo))
        mid = 0.5 * (blo[ar, dim] + bhi[ar, dim])
        b1hi = bhi.copy(); b1hi[ar, dim] = mid
        b2lo = blo.copy(); b2lo[ar, dim] = mid
        blo = np.concatenate([blo, b2lo]); bhi = np.concatenate([b1hi, bhi])
    if len(blo):
        unresolved.append(np.concatenate([blo, bhi], 1))
    unres = np.concatenate(unresolved) if unresolved else np.zeros((0, 8))
    return st, worst, unres


def initial_boxes(n0=12):
    simplex_lo, simplex_hi = simplex_boxes(n0)
    a = np.linspace(2 * 3.141592653589793 / 3, 3.141592653589793, n0 + 1)
    a[-1] = PI_HI  # make sure theta3 = pi is covered
    a[0] = np.nextafter(2 * 3.141592653589793 / 3, 0)
    b = np.linspace(-1.0, 1.0, n0 + 1)
    L, H = [], []
    for i in range(n0):
        for j in range(n0):
            for low, high in zip(simplex_lo, simplex_hi):
                L.append([a[i], b[j], *low])
                H.append([a[i + 1], b[j + 1], *high])
    return np.array(L), np.array(H)


if __name__ == '__main__':
    s = float(sys.argv[1]); nw = int(sys.argv[2]); maxit = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    n0 = int(sys.argv[4]) if len(sys.argv) > 4 else 12
    t0 = time.time()
    L, H = initial_boxes(n0)
    chunks = np.array_split(np.random.default_rng(n0).permutation(len(L)), nw * 8)
    jobs = [(L[c], H[c], s, maxit, 1e-9) for c in chunks]
    with Pool(nw) as pool:
        res = pool.map(process, jobs, chunksize=1)
    st = {k: sum(r[0][k] for r in res) for k in res[0][0]}
    worst = dict(jac=max(r[1]['jac'] for r in res), Ih=min(r[1]['Ih'] for r in res), turn=max(r[1]['turn'] for r in res), W=max(r[1]['W'] for r in res))
    unres = np.concatenate([r[2] for r in res])
    tag = f's{s}' + ('' if n0 == 12 else f'_n{n0}')
    np.save(f'unresolved_{tag}.npy', unres)
    summary = dict(s=s, n0=n0, initial_boxes=int(len(L)), stats=st, worst_certified=worst,
                   unresolved=int(len(unres)), seconds=time.time() - t0)
    summary.update(metadata(['certify.py', 'verify.py', 'interval.py', 'ad.py', 'domains.py'],
                            s=s, n0=n0, workers=nw, maxit=maxit, min_width=1e-9))
    summary['worst_certified'] = {k: v if np.isfinite(v) else None for k, v in worst.items()}
    json.dump(summary, open(f'cert_{tag}.json', 'w'), indent=1)
    print(json.dumps(summary, indent=1))
    if len(unres):
        raise SystemExit(1)
