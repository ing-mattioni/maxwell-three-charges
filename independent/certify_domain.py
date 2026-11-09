"""Independent certification of items 1-2 of Theorem 5.1 (s = 1), from scratch.

Usage:  python independent/certify_domain.py WORKERS [GRID]
Writes independent_domain.json in the current directory.  Exit status 1 if a box is left
unresolved.  Only NumPy is required.

Own code throughout: evaluator.py (plain interval evaluation, no automatic differentiation,
no mean-value forms), own initial grid, own bisection rule.  The domain is
    theta3 in [a0, b0] with a0 < 2pi/3 and b0 > pi,   tau in [-1, 1],   (lambda1, lambda2) in [0, 1]^2,
which contains the domain D of the manuscript.  A box is closed when one of these holds on
the whole box, as a strict inequality between enclosure endpoints (so a NaN closes nothing):
    (o)   lambda3 < 0                       (the box lies outside D);
    (i)   |Z|^2 > 1/9                       (no point of the closed minimum region);
    (ii)  |Z|^2 < 1/9 and sigma3^2 J < 0    (interior points; (J) holds);
    (iii) sigma3^2 J < 0, Icheck > 0, |rho_1 Im t| < 1/2, and either nu < 1/2 or no cusp is
          possible (Im(conj(t)^2 Z) of constant strict sign, or Re(conj(t)^2 Z) > 0).
Every other box is bisected at the floating-point midpoint of one coordinate, so the two
halves share a face.  Boxes narrower than 1e-13 in every coordinate count as unresolved.
"""
import sys
import os
import time
import json
import hashlib
import platform
import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import evaluator as ev
from evaluator import Iv

WEIGHTS = np.array([1.0, 0.8, 1.6, 1.6])
MIN_W = 1e-13
BATCH = 60000
MAX_DEPTH = 160


def step(lo, hi, st):
    """one round on a batch: returns the boxes that are still open"""
    y = [Iv(lo[:, j], hi[:, j]) for j in range(4)]
    l3 = 1.0 - y[2] - y[3]
    dom = l3.b < 0
    Z, _ = ev.zeta(*y)
    z2 = ev.cabs2(Z)
    out = z2.a > ev.NINTH.b
    closed = dom | out
    st['o'] += int(dom.sum()); st['i'] += int((out & ~dom).sum())
    lo, hi = lo[~closed], hi[~closed]
    if len(lo) == 0:
        return lo, hi
    q = ev.quantities(*(Iv(lo[:, j], hi[:, j]) for j in range(4)))
    jneg = q['jac'].b < 0
    ok2 = (q['z2'].b < ev.NINTH.a) & jneg
    nocusp = (q['czim'].a > 0) | (q['czim'].b < 0) | (q['czre'].a > 0)
    ok3 = jneg & (q['Ich'].a > 0) & (q['W'].a > -0.5) & (q['W'].b < 0.5) & ((q['nu'].b < 0.5) | nocusp)
    st['ii'] += int(ok2.sum()); st['iii'] += int((ok3 & ~ok2).sum())
    done = ok2 | ok3
    return lo[~done], hi[~done]


def check(lo, hi):
    """branch and bound on the given boxes until every sub-box is closed or unresolved"""
    st = dict(evals=0, o=0, i=0, ii=0, iii=0, unresolved=0, maxdepth=0)
    depth = 0
    pend_lo, pend_hi = lo, hi
    while len(pend_lo):
        st['evals'] += len(pend_lo)
        nl, nh = [], []
        for k in range(0, len(pend_lo), BATCH):
            with np.errstate(all='ignore'):
                a, b = step(pend_lo[k:k + BATCH], pend_hi[k:k + BATCH], st)
            if len(a):
                nl.append(a); nh.append(b)
        if not nl:
            break
        lo, hi = np.concatenate(nl), np.concatenate(nh)
        w = hi - lo
        sc = np.where(w > MIN_W, w * WEIGHTS, -1.0)
        dim = np.argmax(sc, axis=1)
        tiny = sc.max(axis=1) < 0
        st['unresolved'] += int(tiny.sum())
        lo, hi, dim = lo[~tiny], hi[~tiny], dim[~tiny]
        ar = np.arange(len(lo))
        mid = 0.5 * (lo[ar, dim] + hi[ar, dim])
        h1 = hi.copy(); h1[ar, dim] = mid
        l2 = lo.copy(); l2[ar, dim] = mid
        pend_lo, pend_hi = np.concatenate([lo, l2]), np.concatenate([h1, hi])
        depth += 1
        st['maxdepth'] = depth
        if depth > MAX_DEPTH:
            st['unresolved'] += len(pend_lo)
            break
    return st


def initial(grid):
    a0 = float(np.nextafter(np.nextafter(2 * 3.141592653589793 / 3, 0.0), 0.0))
    b0 = float(ev.PI.b)
    # a0 < 2pi/3 and b0 > pi, established with the interval arithmetic itself
    if not ((Iv(a0) * 3.0).b < ev.TWO_PI.a and ev.PI.a < b0):
        raise ArithmeticError('the theta3 range does not provably contain [2pi/3, pi]')
    e1 = a0 + (b0 - a0) * (np.arange(grid + 1) / grid); e1[0] = a0; e1[-1] = b0
    e2 = -1.0 + 2.0 * (np.arange(grid + 1) / grid); e2[0] = -1.0; e2[-1] = 1.0
    e3 = np.arange(grid + 1) / grid                      # exact when grid is a power of two
    if not (np.all(np.diff(e1) > 0) and np.all(np.diff(e2) > 0)):
        raise ArithmeticError('degenerate initial grid')
    L, H = [], []
    for i in range(grid):
        for j in range(grid):
            for k in range(grid):
                for m in range(grid):
                    # consecutive cells use the same floating-point value as a shared endpoint
                    L.append([e1[i], e2[j], e3[k], e3[m]]); H.append([e1[i + 1], e2[j + 1], e3[k + 1], e3[m + 1]])
    return np.array(L), np.array(H), a0, b0


def work(job):
    lo, hi = job
    t0 = time.time()
    st = check(lo, hi)
    st['initial'] = int(len(lo)); st['seconds'] = time.time() - t0
    return st


def provenance(names):
    """runtime versions and SHA-256 of the source files, named relative to the repository root"""
    hashes = {}
    for n in sorted(names):
        with open(os.path.join(HERE, n), 'rb') as f:
            hashes['independent/' + n] = hashlib.sha256(f.read()).hexdigest()
    return dict(runtime=dict(python=sys.version, numpy=np.__version__, platform=platform.platform()),
                source_sha256=hashes)


if __name__ == '__main__':
    nproc = int(sys.argv[1]); grid = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    if grid < 1 or grid & (grid - 1):
        raise SystemExit('GRID must be a power of two (exact lambda endpoints)')
    L, H, a0, b0 = initial(grid)
    perm = np.random.default_rng(1).permutation(len(L))
    chunks = np.array_split(perm, min(len(L), nproc * 16))
    t0 = time.time()
    res = []
    with Pool(nproc) as pool:
        for r in pool.imap_unordered(work, [(L[c], H[c]) for c in chunks], chunksize=1):
            res.append(r)
            print(json.dumps(r), flush=True)
    tot = {k: sum(r[k] for r in res) for k in ('initial', 'evals', 'o', 'i', 'ii', 'iii', 'unresolved')}
    closed = tot['o'] + tot['i'] + tot['ii'] + tot['iii']
    summary = dict(
        schema_version=1,
        description='independent certification of (J), (I), (T), (W) on the chart domain, s = 1',
        s=1.0, initial_boxes=tot['initial'], boxes=tot['evals'], unresolved=tot['unresolved'],
        closed=dict(outside_domain=tot['o'], outside_minimum_region=tot['i'], interior=tot['ii'], boundary=tot['iii']),
        max_depth=max(r['maxdepth'] for r in res),
        binary_tree_identity=bool(tot['evals'] == 2 * (closed + tot['unresolved']) - tot['initial']),
        theta3_range=[a0, b0], seconds=time.time() - t0,
        parameters=dict(workers=nproc, grid=grid, split_weights=WEIGHTS.tolist(), min_width=MIN_W, max_depth=MAX_DEPTH))
    summary.update(provenance(['evaluator.py', 'certify_domain.py']))
    json.dump(summary, open('independent_domain.json', 'w'), indent=1)
    print(json.dumps(summary, indent=1), flush=True)
    if summary['unresolved'] or not summary['binary_tree_identity']:
        raise SystemExit(1)
