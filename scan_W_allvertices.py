"""Check whether the (W) bound |d arg u_j / d omega| < 1/2 holds on dM for ALL three vertices j,
not only for the vertex labelled 1 in the max-angle frame.

Samples the chart directly: (theta3, tau) uniform in the certificate domain, lambda on random
segments of the simplex, solving |zeta| = kappa along each segment (|zeta|^2 is quadratic).
d arg u_j / d omega = -Im(t conj(v_j)) / r_j, and with t = sigma3 * th (closedform frame):
   j=1,2:  sigma3 * Im(th conj v_j) * lambda_j / sigma_j,     j=3:  Im(th conj v_3) * lambda_3.
"""
import sys
import numpy as np
from closedform import quantities

rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
N = int(sys.argv[2]) if len(sys.argv) > 2 else 2_000_000
s = 1.0
kappa = s / (s + 2)

th3 = rng.uniform(2 * np.pi / 3, np.pi, N)
tau = rng.uniform(-1, 1, N)
th1 = (2 * np.pi - th3) / 2 + tau * (3 * th3 - 2 * np.pi) / 2
th2 = (2 * np.pi - th3) / 2 - tau * (3 * th3 - 2 * np.pi) / 2
la = rng.dirichlet([1, 1, 1], N)
lb = rng.dirichlet([1, 1, 1], N)
v = np.stack([np.ones(N) + 0j, np.exp(1j * th3), np.exp(-1j * th2)], 1)
za = np.sum(la * v ** 2, 1)
zb = np.sum(lb * v ** 2, 1)
d = zb - za
# |za + t d|^2 = kappa^2
A = np.abs(d) ** 2
B = 2 * (za * d.conjugate()).real
Cc = np.abs(za) ** 2 - kappa ** 2
disc = B ** 2 - 4 * A * Cc
rows = []
for sgn in (-1, 1):
    t = (-B + sgn * np.sqrt(np.maximum(disc, 0))) / (2 * A)
    ok = (disc > 0) & (t > 0) & (t < 1)
    lam = la[ok] + t[ok, None] * (lb[ok] - la[ok])
    rows.append((th1[ok], th2[ok], th3[ok], lam))
th1 = np.concatenate([r[0] for r in rows]); th2 = np.concatenate([r[1] for r in rows])
th3 = np.concatenate([r[2] for r in rows]); lam = np.concatenate([r[3] for r in rows])
q = quantities(th1, th2, lam[:, 0], lam[:, 1], s)
th = q['t']
s1, s2, s3 = np.sin(th1), np.sin(th2), np.sin(th3)
v = [np.ones_like(th) , np.exp(1j * th3), np.exp(-1j * th2)]
W1 = np.abs(s3 * (th * v[0].conjugate()).imag * lam[:, 0] / s1)
W2 = np.abs(s3 * (th * v[1].conjugate()).imag * lam[:, 1] / s2)
W3 = np.abs((th * v[2].conjugate()).imag * lam[:, 2])
print('boundary samples:', len(th))
for name, W in (('W1', W1), ('W2', W2), ('W3', W3)):
    i = np.nanargmax(W)
    print(f'{name}: max {W[i]:.4f}  frac>0.5 {np.mean(W > 0.5):.2e}  at th=({th1[i]:.4f},{th2[i]:.4f},{th3[i]:.4f}) lam={lam[i]}')
Wmax = np.maximum(np.maximum(W1, W2), W3)
print('max over all vertices:', np.nanmax(Wmax))
# per-sample: does at least one of the two far vertices (1,2) satisfy the bound?
print('max of min(W1,W2):', np.nanmax(np.minimum(W1, W2)))
