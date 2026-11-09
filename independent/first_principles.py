"""Closed forms of Section 6 against first principles, in 50-digit arithmetic.

Usage:  python independent/first_principles.py [POINTS_PER_REGIME]      (needs mpmath)

For a chart point on dM (|Z| = kappa to working precision) the actual triangle
p_i = r_i v_i is built and the following are computed WITHOUT any formula of Section 6:
   Z, Z_z, Z_zbar by finite differences of the definition, and J;
   x0(omega) by Newton's method on Z(x) = kappa e^{i omega}; t = dx0/domega; nu = d arg t/domega;
   m(omega) from its definition, det(1, m, m'), d arg u_1/domega;
   the gradient and the Hessian of V_q for q = F(x0).
They are compared with a direct transcription of (6.1) and (6.4)-(6.10) of the manuscript.
The agreement is limited by the finite differences (about 30 digits).  This is a numerical
check of the formulas at sample points, not a proof.
"""
import sys
import random
import mpmath as mp

mp.mp.dps = 50
I_ = mp.mpc(0, 1)


def dot(a, b):
    return (a * mp.conj(b)).real


# ----------------------------------------------------------------------------- first principles
def bary(x, p):
    A = (mp.conj(p[1] - p[0]) * (p[2] - p[0])).imag
    return [(mp.conj(p[(i + 1) % 3] - x) * (p[(i + 2) % 3] - x)).imag / A for i in range(3)]


def Zdef(x, p):
    lam = bary(x, p)
    return sum(lam[i] * ((x - p[i]) / abs(x - p[i])) ** 2 for i in range(3))


STEN1 = [(-3, mp.mpf(-1) / 60), (-2, mp.mpf(3) / 20), (-1, mp.mpf(-3) / 4), (1, mp.mpf(3) / 4), (2, mp.mpf(-3) / 20), (3, mp.mpf(1) / 60)]
STEN2 = [(-3, mp.mpf(1) / 90), (-2, mp.mpf(-3) / 20), (-1, mp.mpf(3) / 2), (0, mp.mpf(-49) / 18), (1, mp.mpf(3) / 2), (2, mp.mpf(-3) / 20), (3, mp.mpf(1) / 90)]


def d1(f, h):
    return sum(c * f(j * h) for j, c in STEN1) / h


def grad_Z(x, p, h=mp.mpf(10) ** -7):
    Zx = d1(lambda e: Zdef(x + e, p), h)
    Zy = d1(lambda e: Zdef(x + I_ * e, p), h)
    return (Zx - I_ * Zy) / 2, (Zx + I_ * Zy) / 2      # Z_z, Z_zbar


def solve_x0(target, p, x):
    for _ in range(60):
        F = Zdef(x, p) - target
        if abs(F) < mp.mpf(10) ** -44:
            return x
        Zz, Zzb = grad_Z(x, p)
        J = abs(Zz) ** 2 - abs(Zzb) ** 2
        x = x - (mp.conj(Zz) * F - Zzb * mp.conj(F)) / J
    raise RuntimeError('Newton failed')


def first_principles(p, kappa, om, s):
    h = mp.mpf(10) ** -5
    xs = {}
    x = mp.mpc(0)
    for j in (0, 1, 2, 3):
        x = solve_x0(kappa * mp.e ** (I_ * (om + j * h)), p, x)
        xs[j] = x
    x = xs[0]
    for j in (-1, -2, -3):
        x = solve_x0(kappa * mp.e ** (I_ * (om + j * h)), p, x)
        xs[j] = x
    t = sum(c * xs[j] for j, c in STEN1) / h
    tp = sum(c * xs[j] for j, c in STEN2) / h ** 2
    nu = (tp / t).imag

    def mvec(j):
        k = I_ * mp.e ** (I_ * (om + j * h) / 2)
        out = []
        for i in range(3):
            d = xs[j] - p[i]
            r = abs(d)
            out.append(-s * r ** (-s - 1) * dot(d / r, k))
        return out
    m0 = mvec(0)
    ms = {j: mvec(j) for j in (-3, -2, -1, 1, 2, 3)}
    mprime = [sum(c * ms[j][i] for j, c in STEN1) / h for i in range(3)]
    det = mp.det(mp.matrix([[1, 1, 1], m0, mprime]))
    Zz, Zzb = grad_Z(xs[0], p)
    J = abs(Zz) ** 2 - abs(Zzb) ** 2
    darg_u1 = (t / (xs[0] - p[0])).imag
    return dict(Z=Zdef(xs[0], p), Zz=Zz, Zzb=Zzb, J=J, t=t, nu=nu, det=det, W=darg_u1, m=m0, mprime=mprime)


# ----------------------------------------------------------------------------- manuscript formulas
def paper_formulas(th, lam, s):
    """Direct transcription of (6.1), (6.4)-(6.10) of the manuscript (unrescaled)."""
    th1, th2, th3 = th
    sig = [mp.sin(a) for a in th]
    v = [mp.mpc(1), mp.e ** (I_ * th3), mp.e ** (-I_ * th2)]
    rho = [lam[i] / sig[i] for i in range(3)]
    mu = []
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        mu.append(I_ / 2 * rho[i] * (rho[k] * mp.conj(v[j]) - rho[j] * mp.conj(v[k])))
    Z = sum(lam[i] * v[i] ** 2 for i in range(3))
    Zz = sum(mu[i] * v[i] ** 2 - lam[i] * rho[i] * v[i] for i in range(3))
    Zzb = sum(mp.conj(mu[i]) * v[i] ** 2 + lam[i] * rho[i] * v[i] ** 3 for i in range(3))
    Zzz = -2 * sum(mu[i] * rho[i] * v[i] for i in range(3))
    Zzzb = sum(mu[i] * rho[i] * v[i] ** 3 - mp.conj(mu[i]) * rho[i] * v[i] - lam[i] * rho[i] ** 2 * v[i] ** 2 for i in range(3))
    Zzbzb = 2 * sum(mp.conj(mu[i]) * rho[i] * v[i] ** 3 + lam[i] * rho[i] ** 2 * v[i] ** 4 for i in range(3))
    J = abs(Zz) ** 2 - abs(Zzb) ** 2
    t = (mp.conj(Zz) * (I_ * Z) - Zzb * mp.conj(I_ * Z)) / J
    R = -Z - (Zzz * t ** 2 + 2 * Zzzb * abs(t) ** 2 + Zzbzb * mp.conj(t) ** 2)
    tp = (mp.conj(Zz) * R - Zzb * mp.conj(R)) / J
    nu = (tp / t).imag
    kb2 = -mp.conj(Z) / abs(Z)

    def Q(a, b):
        return (a * mp.conj(b)).real / 2 + (a * b * kb2).real / 2
    tc = t / sig[2]
    Ich = sig[0] * sig[1] / 2 + (s + 2) * Q(v[0], v[2]) * (lam[0] * sig[2] * dot(v[0], tc) - sig[0] * lam[2] * dot(v[2], tc)) \
        - lam[0] * sig[2] * Q(v[2], tc) + sig[0] * lam[2] * Q(v[0], tc)
    # (6.8) directly
    om = mp.arg(Z)
    k = I_ * mp.e ** (I_ * om / 2)
    kp = I_ / 2 * k
    u = [-vi for vi in v]
    a = [sig[i] * dot(u[i], k) for i in range(3)]
    b = [sig[i] * dot(u[i], kp) - lam[i] * ((s + 2) * dot(u[i], t) * dot(u[i], k) - dot(t, k)) for i in range(3)]
    W = -rho[0] * t.imag
    return dict(Z=Z, Zz=Zz, Zzb=Zzb, J=J, t=t, nu=nu, Ich=Ich, a=a, b=b, W=W, sig=sig, rho=rho, v=v, cz=mp.conj(t) ** 2 * Z)


def boundary_point(th3, tau, om, kappa):
    """chart point with Z = kappa e^{i omega}: the barycentric coordinates are those of
    kappa e^{i omega} in the triangle (v_1^2, v_2^2, v_3^2)"""
    th1 = (2 * mp.pi - th3) / 2 + tau * (3 * th3 - 2 * mp.pi) / 2
    th2 = (2 * mp.pi - th3) / 2 - tau * (3 * th3 - 2 * mp.pi) / 2
    w = [mp.mpc(1), mp.e ** (2 * I_ * th3), mp.e ** (-2 * I_ * th2)]
    zeta = kappa * mp.e ** (I_ * om)
    den = (mp.conj(w[1] - w[2]) * (w[0] - w[2])).imag
    l1 = (mp.conj(w[1] - w[2]) * (zeta - w[2])).imag / den
    l2 = (mp.conj(zeta - w[2]) * (w[0] - w[2])).imag / den
    return (th1, th2, th3), [l1, l2, 1 - l1 - l2]


def rel(a, b):
    return abs(a - b) / max(abs(a), abs(b), mp.mpf(10) ** -30)


def run(npts=40, s=mp.mpf(1), seed=1, flat=False):
    rng = random.Random(seed)
    kappa = s / (s + 2)
    worst = {}
    done = 0
    while done < npts:
        if flat:
            th3 = mp.pi - mp.mpf(10) ** (-rng.uniform(1, 6))
        else:
            th3 = mp.mpf(rng.uniform(2.0944, 3.1415))
        tau = mp.mpf(rng.uniform(-1, 1))
        om = mp.mpf(rng.uniform(0, 6.2832))
        th, lam = boundary_point(th3, tau, om, kappa)
        if min(lam) <= mp.mpf('1e-4') or min(th) <= 0:
            continue
        P = paper_formulas(th, lam, s)
        p = [P['sig'][i] / lam[i] * P['v'][i] for i in range(3)]
        try:
            fp = first_principles(p, kappa, mp.arg(P['Z']), s)
        except RuntimeError:
            print('newton failed, skipping'); continue
        done += 1
        errs = dict(
            Z=rel(fp['Z'], P['Z']), absZ=abs(abs(P['Z']) - kappa), Zz=rel(fp['Zz'], P['Zz']), Zzb=rel(fp['Zzb'], P['Zzb']),
            J=rel(fp['J'], P['J']), t=rel(fp['t'], P['t']), nu=abs(fp['nu'] - P['nu']), W=abs(fp['W'] - P['W']),
        )
        # (6.8): D m = -s a, D m' = -s b with D = diag(lam_i r_i^{s+2})
        q = [lam[i] * (P['sig'][i] / lam[i]) ** (s + 2) for i in range(3)]
        errs['a'] = max(abs(q[i] * fp['m'][i] + s * P['a'][i]) for i in range(3))
        errs['b'] = max(abs(q[i] * fp['mprime'][i] + s * P['b'][i]) for i in range(3))
        errs['sum_a'] = abs(sum(P['a'])); errs['sum_b'] = abs(sum(P['b']))
        # det(1,a,b) = 3 sigma3 Icheck ; det(1,m,m') = s^2 sigma3 Icheck sum(q)/prod(q)
        deta = mp.det(mp.matrix([[1, 1, 1], P['a'], P['b']]))
        errs['det_ab'] = rel(deta, 3 * P['sig'][2] * P['Ich'])
        errs['det_m'] = rel(fp['det'], s ** 2 * P['sig'][2] * P['Ich'] * sum(q) / (q[0] * q[1] * q[2]))
        # Lami + Hessian: gradient of V_q at x0 = 0 and kernel of the Hessian

        def V(x):
            return sum(q[i] * abs(x - p[i]) ** (-s) for i in range(3))
        h = mp.mpf(10) ** -7
        gx = d1(lambda e: V(mp.mpc(e, 0)), h); gy = d1(lambda e: V(mp.mpc(0, e)), h)
        k = I_ * mp.e ** (I_ * mp.arg(P['Z']) / 2)

        def dd(dirn):
            return sum(c * V(j * h * dirn) for j, c in STEN2) / h ** 2
        Hkk, Hnn = dd(k), dd(I_ * k)
        scale = sum(q[i] * (P['sig'][i] / lam[i]) ** (-s - 2) for i in range(3))
        errs['grad'] = abs(mp.mpc(gx, gy)) / scale
        errs['Hkk'] = abs(Hkk) / scale
        errs['Hnn_pos'] = mp.mpf(0) if Hnn > 0 else mp.mpf(1)
        for kk, e in errs.items():
            worst[kk] = max(worst.get(kk, 0), e)
    return worst


if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 25
    bad = False
    for label, kw in (('generic points, s = 1', dict()), ('near the flat face, s = 1', dict(flat=True, seed=2)),
                      ('generic points, s = 0.7', dict(s=mp.mpf('0.7'), seed=3))):
        w = run(n, **kw)
        print(label)
        for k, v in w.items():
            print(f'   max error {k:8s} {mp.nstr(v, 5)}')
        bad = bad or max(w.values()) > mp.mpf(10) ** -25
    if bad:
        raise SystemExit('disagreement above 1e-25')
    print('closed forms agree with first principles to better than 1e-25 at all sampled points')
