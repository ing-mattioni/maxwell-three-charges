"""Sanity test: interval/AD mean-value enclosures must contain point evaluations of the
independent float implementation (closedform.py) at random points of random boxes."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import numpy as np
from certify import enclose_mvf, zeta2_naive
from closedform import quantities
from verify import thetas_from_y
from interval import I


def test_containment(seed=0, nbox=3000, npts=10, s=1.0):
    rng = np.random.default_rng(seed)
    blo = np.column_stack([rng.uniform(2.1, 3.1415, nbox), rng.uniform(-1, 1, nbox),
                           rng.uniform(0, .6, nbox), rng.uniform(0, .6, nbox)])
    w = 10 ** rng.uniform(-5, -2, (nbox, 1)) * np.array([1, 2, 1, 1])
    bhi = blo + w
    with np.errstate(all='ignore'):
        e, _ = enclose_mvf(blo, bhi, s)
        zn = zeta2_naive(blo, bhi)
    bad = 0
    for _ in range(npts):
        p = blo + (bhi - blo) * rng.uniform(size=blo.shape)
        t1, t2, _t3 = thetas_from_y(I(p[:, 0]), I(p[:, 1]))
        with np.errstate(all='ignore'):
            q = quantities(t1.lo, t2.lo, p[:, 2], p[:, 3], s)
        vals = dict(jac=q['jac'], Ih=q['Ih'], turn=q['turn'], zeta2=np.abs(q['zeta']) ** 2)
        for k, v in vals.items():
            m = np.isfinite(v) & np.isfinite(e[k].lo) & np.isfinite(e[k].hi)
            bad += int((m & ~((e[k].lo <= v) & (v <= e[k].hi))).sum())
        z = vals['zeta2']
        bad += int((~((zn.lo <= z) & (z <= zn.hi))).sum())
    assert bad == 0, f'{bad} containment violations'
    return True


if __name__ == '__main__':
    for sd in range(3):
        for s in (1.0, 0.5, 1.4):
            test_containment(sd, s=s)
    print('all containment tests passed')
