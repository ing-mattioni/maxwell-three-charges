"""Plot the minimum region M (det T > (s+1)/(s+2)^2) and its image (caustic) in charge space."""
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from core import triangle, barycentric, charges_at, detT

shapes = [(0.5, np.sqrt(3) / 2), (0.5, 0.4), (0.5, 1.6), (0.2, 0.6), (0.05, 0.3), (0.3, 0.95)]
svals = [0.3, 1.0, 3.0, 10.0]
fig, axes = plt.subplots(len(shapes), 2, figsize=(9, 3.2 * len(shapes)))
for row, (a, b) in enumerate(shapes):
    P = triangle(a, b)
    xs = np.linspace(min(0, a) - 0.01, max(1, a) + 0.01, 700)
    ys = np.linspace(0.0, b, 700)
    Xg, Yg = np.meshgrid(xs, ys)
    X = np.stack([Xg, Yg], -1)
    lam = barycentric(P, X)
    D = detT(P, X)
    D[lam.min(-1) <= 0] = np.nan
    ax, bx = axes[row]
    ax.plot(*np.vstack([P, P[:1]]).T, 'k-')
    for s, col in zip(svals, ['C0', 'C1', 'C2', 'C3']):
        thr = (s + 1) / (s + 2) ** 2
        cs = ax.contour(Xg, Yg, D, levels=[thr], colors=col)
        for seg in cs.allsegs[0]:
            if len(seg) < 3:
                continue
            q = charges_at(P, seg, s)
            bx.plot(np.log(q[:, 0] / q[:, 2]), np.log(q[:, 1] / q[:, 2]), color=col, lw=1)
    ax.set_aspect('equal'); ax.set_title(f'M, a={a:.2f}, b={b:.2f}')
    bx.set_title('caustic F(dM): log q1/q3 vs log q2/q3')
plt.tight_layout()
plt.savefig(sys.argv[1] if len(sys.argv) > 1 else 'caustic.png', dpi=80)
