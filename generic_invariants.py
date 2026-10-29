"""Chart-free evaluation of the boundary invariants at a point x0 of dM for ANY sign pattern.

Input: triangle vertices p (complex, shape (3,)), boundary points z (complex array), exponent s.
All formulas use the actual geometry: barycentric gradients, Z and its first/second
derivatives, t = dx0/domega, the kernel direction k, and the rescaled dual vectors
m~_i = q0_i m_i (q0 = lambda_i r_i^{s+2}, possibly with a negative entry).
Returns jac, In (normalised inflection function), turn = d arg t/d omega, cusp indicator,
and W_j = |d arg u_j / d omega| for j = 0,1,2.
"""
import numpy as np


def bary_grad(p):
    """d lambda_i / dz (complex Wirtinger derivative), constant for a triangle"""
    p1, p2, p3 = p
    den = ((p2 - p3).conjugate() * (p1 - p3)).imag
    dl1 = (p2 - p3).conjugate() / (2j * den)
    dl2 = -(p1 - p3).conjugate() / (2j * den)
    return np.array([dl1, dl2, -dl1 - dl2])


def bary(p, z):
    p1, p2, p3 = p
    den = ((p2 - p3).conjugate() * (p1 - p3)).imag
    l1 = ((p2 - p3).conjugate() * (z - p3)).imag / den
    l2 = ((z - p3).conjugate() * (p1 - p3)).imag / den
    return np.stack([l1, l2, 1 - l1 - l2], -1)


def invariants(p, z, s):
    mu = bary_grad(p)
    lam = bary(p, z)
    d = z[:, None] - p
    db = d.conjugate()
    r = np.abs(d)
    u = d / r
    u2 = d / db
    Z = np.sum(lam * u2, 1)
    Zz = np.sum(mu * u2 + lam / db, 1)
    Zzb = np.sum(mu.conjugate() * u2 - lam * d / db ** 2, 1)
    jac = np.abs(Zz) ** 2 - np.abs(Zzb) ** 2
    iZ = 1j * Z
    t = (Zz.conjugate() * iZ - Zzb * iZ.conjugate()) / jac
    Zzz = np.sum(2 * mu / db, 1)
    Zzzb = np.sum(-mu * d / db ** 2 + mu.conjugate() / db - lam / db ** 2, 1)
    Zzbzb = np.sum(-2 * mu.conjugate() * d / db ** 2 + 2 * lam * d / db ** 3, 1)
    rhs = -Z - (Zzz * t ** 2 + 2 * Zzzb * np.abs(t) ** 2 + Zzbzb * t.conjugate() ** 2)
    tp = (Zz.conjugate() * rhs - Zzb * rhs.conjugate()) / jac
    turn = (tp / t).imag
    k = 1j * np.exp(0.5j * np.angle(Z))
    kp = 0.5j * k
    dot = lambda a, b: (a * b.conjugate()).real
    q0 = lam * r ** (s + 2)
    # m_i = k . grad phi_i  with grad phi_i = -s r^{-s-1} u_i ;  m~ = q0 * m / (-s)
    mt = lam * r * dot(u, k[:, None])
    mtp = lam * r * dot(u, kp[:, None]) - lam * ((s + 2) * dot(u, t[:, None]) * dot(u, k[:, None]) - dot(t, k)[:, None])
    one = np.ones_like(mt)
    I = np.linalg.det(np.stack([one, mt, mtp], -2))
    In = I / (np.linalg.norm(mt, axis=1) * np.linalg.norm(mtp, axis=1) * np.sqrt(3))
    cusp = (t * k.conjugate()).imag / np.abs(t)
    W = np.abs((t[:, None] * u.conjugate()).imag) / r
    return dict(Z=Z, jac=jac, In=In, turn=turn, cusp=cusp, W=W, t=t, lam=lam, q0=q0)
