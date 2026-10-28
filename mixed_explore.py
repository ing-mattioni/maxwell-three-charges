"""Mixed signs q = (+,+,-): extrema live in R_A = {l1,l2>0>l3} (minima) and R_B = {l1,l2<0<l3} (maxima),
in both cases where |Z| < kappa.  Count extrema directly for charges sampled from F(M_A) and F(M_B)."""
import sys
import numpy as np
from fast import tri, bary, Zfun


def field(p, q, s, z):
    d = z[..., None] - p
    r = np.abs(d)
    return np.sum(q * r ** (-s - 2) * d, axis=-1)


def crit_points(p, q, s, R=6.0, n=160):
    """all planar critical points of sum q_i r_i^-s in a disc of radius R (multistart Newton)"""
    xs = np.linspace(-R, R, n)
    X, Y = np.meshgrid(xs, xs)
    z = (X + 1j * Y).ravel() + p.mean()
    h = 1e-7
    for _ in range(60):
        f = field(p, q, s, z)
        fx = (field(p, q, s, z + h) - field(p, q, s, z - h)) / (2 * h)
        fy = (field(p, q, s, z + 1j * h) - field(p, q, s, z - 1j * h)) / (2 * h)
        # Jacobian of (Re f, Im f) wrt (x, y)
        a, b, c, d = fx.real, fy.real, fx.imag, fy.imag
        det = a * d - b * c
        det = np.where(np.abs(det) < 1e-300, 1e-300, det)
        dx = (d * f.real - b * f.imag) / det
        dy = (-c * f.real + a * f.imag) / det
        step = dx + 1j * dy
        nr = np.abs(step)
        step = np.where(nr > 0.2, step * 0.2 / np.maximum(nr, 1e-300), step)
        z = z - step
        z = np.where(np.abs(z - p.mean()) > 4 * R, np.nan, z)
    f = field(p, q, s, z)
    r = np.abs(z[:, None] - p).min(axis=1)
    ok = np.isfinite(z) & (np.abs(f) * r ** (s + 1) < 1e-9) & (r > 1e-6)
    z = z[ok]
    out = []
    for w in z:
        if all(abs(w - u) > 1e-6 for u in out):
            out.append(w)
    return np.array(out)


if __name__ == '__main__':
    s = 1.0; kappa = s / (s + 2)
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    worst = 0
    for trial in range(int(sys.argv[2]) if len(sys.argv) > 2 else 30):
        a, b = rng.uniform(-0.5, 1.5), 10 ** rng.uniform(-1, 0.3)
        p = tri(a, b)
        xs = np.linspace(-4, 5, 500)
        X, Y = np.meshgrid(xs, np.linspace(-5, 5, 500))
        z = (X + 1j * Y).ravel()
        lam = bary(p, z)
        Z = Zfun(p, z)
        A = (lam[:, 0] > 0) & (lam[:, 1] > 0) & (lam[:, 2] < 0) & (np.abs(Z) < kappa)
        B = (lam[:, 0] < 0) & (lam[:, 1] < 0) & (lam[:, 2] > 0) & (np.abs(Z) < kappa)
        r = np.abs(z[:, None] - p)
        qq = lam * r ** (s + 2)
        qq = qq / qq[:, :2].sum(axis=1, keepdims=True)      # normalise so q1+q2 = 1 (sign class (+,+,-))
        cand = np.concatenate([np.nonzero(A)[0][::max(1, A.sum() // 15)], np.nonzero(B)[0][::max(1, B.sum() // 15)]])
        mx = 0
        for i in cand:
            q = qq[i]
            cps = crit_points(p, q, s)
            mx = max(mx, len(cps))
        worst = max(worst, mx)
        print(f'tri a={a:.3f} b={b:.3f}: |M_A| pts={A.sum()}, |M_B| pts={B.sum()}, max #equilibria over sampled charges = {mx}', flush=True)
    print('overall max', worst)
