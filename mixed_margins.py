"""Margins of the four hedgehog conditions along dM_A for mixed signs (+,+,-)."""
import sys
import numpy as np
from fast import tri
from mixed_caustic import centers_RA, trace_RA
from generic_invariants import invariants

s = 1.0
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
N = int(sys.argv[2]) if len(sys.argv) > 2 else 40
rows = []
for k in range(N):
    a, b = rng.uniform(-3.0, 4.0), 10 ** rng.uniform(-2.5, 0.7)
    p = tri(a, b)
    cs = centers_RA(p, s, ngrid=500, R=12.0)
    if not cs:
        continue
    for c, _ in cs:
        zb, ns = trace_RA(p, s, c, 4096, Rmax=200.0)
        g = invariants(p, zb, s)
        cuspsgn = np.sign(g['cusp'])
        ncusp = int(np.count_nonzero(cuspsgn[1:] != cuspsgn[:-1]) + (cuspsgn[-1] != -cuspsgn[0]))
        rows.append((a, b, len(cs), ns, g['jac'].max(), g['jac'].min(), g['In'].min(), g['In'].max(),
                     g['turn'].max(), g['turn'].min(), g['W'][:, 0].max(), g['W'][:, 1].max(), g['W'][:, 2].max(), ncusp,
                     (-2 * g['lam'][:, 2] + 1).max()))
        print(f"a={a:+.3f} b={b:.4f} comps={len(cs)} nonstar={ns} jac[{g['jac'].min():.2e},{g['jac'].max():.2e}] "
              f"In[{g['In'].min():+.3f},{g['In'].max():+.3f}] turn[{g['turn'].min():+.2f},{g['turn'].max():+.2f}] "
              f"W=({g['W'][:,0].max():.2f},{g['W'][:,1].max():.2f},{g['W'][:,2].max():.2f}) cusps~{ncusp} Lmax={(-2*g['lam'][:,2]+1).max():.1f}", flush=True)
np.save(f'data/mixed_margins_{sys.argv[1] if len(sys.argv) > 1 else 0}.npy', np.array(rows))
