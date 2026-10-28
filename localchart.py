"""Closed-form evaluation of all boundary invariants in the local chart (theta1, theta2, omega).

A point x0 on dM of some triangle is described by
  * the angles theta_i subtended at x0 by the sides (theta_i = angle p_j x0 p_k), sum = 2 pi,
  * omega = arg Z(x0)  (|Z(x0)| = kappa on dM).
Frame: x0 = 0, direction to p_i is v_i = exp(i beta_i), beta_1 = 0, beta_2 = theta_3, beta_3 = theta_3 + theta_1.
lambda = barycentric coords of zeta = kappa e^{i omega} in the triangle (v_1^2, v_2^2, v_3^2);
sigma_i = sin theta_i, r_i = sigma_i / lambda_i (Lami), p_i = r_i v_i.
Everything below is explicit (no root finding).  Works with numpy arrays.
"""
import numpy as np


def config(th1, th2, om, kappa):
    th3 = 2 * np.pi - th1 - th2
    b1 = np.zeros_like(th1)
    b2 = th3
    b3 = th3 + th1
    v = np.stack([np.exp(1j * b1), np.exp(1j * b2), np.exp(1j * b3)], -1)
    w = v ** 2
    zeta = kappa * np.exp(1j * om)
    # barycentric coordinates of zeta in triangle w
    w1, w2, w3 = w[..., 0], w[..., 1], w[..., 2]
    den = ((w2 - w3).conjugate() * (w1 - w3)).imag
    l1 = ((w2 - w3).conjugate() * (zeta - w3)).imag / den
    l2 = ((zeta - w3).conjugate() * (w1 - w3)).imag / den
    lam = np.stack([l1, l2, 1 - l1 - l2], -1)
    sig = np.stack([np.sin(th1), np.sin(th2), np.sin(th3)], -1)
    return v, lam, sig, zeta


def invariants(th1, th2, om, s):
    kappa = s / (s + 2)
    v, lam, sig, Z = config(th1, th2, om, kappa)
    r = sig / lam
    u = -v                  # unit vector from p_i to x0
    d = r * u               # x0 - p_i
    db = d.conjugate()
    # mu_i = d lambda_i / dz : sum mu = 0, sum mu d = -1, sum conj(mu) d = 0
    # real/imag parts:  sum a d = -1/2, sum b d = i/2  (mu = a + i b)
    A = np.stack([np.ones_like(d.real), d.real, d.imag], -2)  # (...,3,3) rows
    rhs_a = np.stack([np.zeros_like(Z.real), -0.5 * np.ones_like(Z.real), np.zeros_like(Z.real)], -1)
    rhs_b = np.stack([np.zeros_like(Z.real), np.zeros_like(Z.real), 0.5 * np.ones_like(Z.real)], -1)
    a = np.linalg.solve(A, rhs_a[..., None])[..., 0]
    b = np.linalg.solve(A, rhs_b[..., None])[..., 0]
    mu = a + 1j * b
    mub = mu.conjugate()
    u2 = u ** 2
    Zz = np.sum(mu * u2 + lam / db, -1)
    Zzb = np.sum(mub * u2 - lam * d / db ** 2, -1)
    jac = np.abs(Zz) ** 2 - np.abs(Zzb) ** 2
    iZ = 1j * Z
    t = (Zz.conjugate() * iZ - Zzb * iZ.conjugate()) / jac   # dx0/domega
    k = 1j * np.exp(1j * om / 2)
    kp = 0.5j * k
    dot = lambda p, q: (p * q.conjugate()).real
    mt = sig * dot(u, k[..., None])
    mtp = sig * dot(u, kp[..., None]) - lam * ((s + 2) * dot(u, t[..., None]) * dot(u, k[..., None]) - dot(t, k)[..., None])
    one = np.ones_like(mt)
    I = np.linalg.det(np.stack([one, mt, mtp], -2))
    # normalised version: sine of angle between mt and mtp
    In = I / (np.linalg.norm(mt, axis=-1) * np.linalg.norm(mtp, axis=-1) * np.sqrt(3))
    # second derivatives of Z
    Zzz = np.sum(2 * mu / db, -1)
    Zzzb = np.sum(-mu * d / db ** 2 + mub / db - lam / db ** 2, -1)
    Zzbzb = np.sum(-2 * mub * d / db ** 2 + 2 * lam * d / db ** 3, -1)
    rhs = -Z - (Zzz * t ** 2 + 2 * Zzzb * np.abs(t) ** 2 + Zzbzb * t.conjugate() ** 2)
    tp = (Zz.conjugate() * rhs - Zzb * rhs.conjugate()) / jac
    turn = (tp / t).imag          # d arg(t) / d omega
    cusp = (t * k.conjugate()).imag / np.abs(t)   # sin(angle between t and k)
    return dict(jac=jac, I=I, In=In, turn=turn, cusp=cusp, t=t, lam=lam, r=r, mt=mt, mtp=mtp)


def domain_mask(th1, th2, om, s):
    kappa = s / (s + 2)
    th3 = 2 * np.pi - th1 - th2
    v, lam, sig, Z = config(th1, th2, om, kappa)
    return (lam.min(-1) > 0) & (th3 > 0) & (th3 < np.pi) & (th1 < np.pi) & (th2 < np.pi)
