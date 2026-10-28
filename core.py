"""Core numerics for Maxwell's problem with three point charges.

Potential in the plane of the charges:  V(x) = sum_i q_i |x - p_i|^{-s}
(s = 1 is Coulomb in R^3; s = 2*alpha in the notation of Gabrielov-Novikov-Shapiro).

Key identities (derived in notes/reformulation.md):
  * x is a critical point for charges q  <=>  q_i ∝ lam_i(x) r_i(x)^{s+2}
    (lam = barycentric coordinates of x w.r.t. the triangle p_1 p_2 p_3).
  * at such a critical point the in-plane Hessian is a positive multiple of
    (s+2) T(x) - I,   T(x) = sum_i lam_i u_i u_i^T,  u_i = (x - p_i)/r_i.
"""
import numpy as np


def triangle(a, b):
    return np.array([[0.0, 0.0], [1.0, 0.0], [a, b]])


def barycentric(P, X):
    """X: (...,2) -> lam: (...,3)"""
    p1, p2, p3 = P
    M = np.array([[p1[0] - p3[0], p2[0] - p3[0]], [p1[1] - p3[1], p2[1] - p3[1]]])
    Minv = np.linalg.inv(M)
    d = X - p3
    l12 = d @ Minv.T
    l3 = 1.0 - l12[..., 0] - l12[..., 1]
    return np.stack([l12[..., 0], l12[..., 1], l3], axis=-1)


def charges_at(P, X, s):
    """Charges (normalised to sum 1) for which X is a critical point."""
    lam = barycentric(P, X)
    r = np.linalg.norm(X[..., None, :] - P, axis=-1)
    q = lam * r ** (s + 2)
    return q / q.sum(axis=-1, keepdims=True)


def Tmat(P, X):
    lam = barycentric(P, X)
    d = X[..., None, :] - P
    r = np.linalg.norm(d, axis=-1)
    u = d / r[..., None]
    T = np.einsum('...i,...ij,...ik->...jk', lam, u, u)
    return T


def detT(P, X):
    T = Tmat(P, X)
    return T[..., 0, 0] * T[..., 1, 1] - T[..., 0, 1] ** 2


def is_min_region(P, X, s):
    """lambda_min(T) > 1/(s+2)  <=>  det T > (s+1)/(s+2)^2 (trace T = 1)."""
    return detT(P, X) > (s + 1) / (s + 2) ** 2


def field(P, q, s, X):
    d = X[..., None, :] - P
    r = np.linalg.norm(d, axis=-1)
    w = q * r ** (-s - 2)
    return np.einsum('...i,...ij->...j', w, d)  # = -grad V / s


def hessV(P, q, s, X):
    d = X[..., None, :] - P
    r = np.linalg.norm(d, axis=-1)
    u = d / r[..., None]
    c = q * s * r ** (-s - 2)
    I = np.eye(2)
    return np.einsum('...i,...ijk->...jk', c, (s + 2) * u[..., :, None] * u[..., None, :] - I)


def critical_points(P, q, s, n=60, tol=1e-12):
    """All critical points inside the triangle by multistart Newton on the field."""
    q = np.asarray(q, float)
    pts = []
    us = np.linspace(0, 1, n + 2)[1:-1]
    starts = []
    for i, a in enumerate(us):
        for b in us[: n - i]:
            c = 1 - a - b
            if c <= 0:
                continue
            starts.append(a * P[0] + b * P[1] + c * P[2])
    X = np.array(starts)
    for _ in range(80):
        F = field(P, q, s, X)
        J = -hessV(P, q, s, X) / s  # Jacobian of field
        det = J[:, 0, 0] * J[:, 1, 1] - J[:, 0, 1] * J[:, 1, 0]
        ok = np.abs(det) > 1e-300
        dx = np.zeros_like(X)
        dx[ok, 0] = (J[ok, 1, 1] * F[ok, 0] - J[ok, 0, 1] * F[ok, 1]) / det[ok]
        dx[ok, 1] = (-J[ok, 1, 0] * F[ok, 0] + J[ok, 0, 0] * F[ok, 1]) / det[ok]
        # damp
        nrm = np.linalg.norm(dx, axis=1)
        dx *= np.minimum(1, 0.1 / np.maximum(nrm, 1e-300))[:, None]
        X = X - dx
        lam = barycentric(P, X)
        bad = (lam.min(axis=1) < -0.5) | ~np.isfinite(X).all(axis=1)
        X[bad] = np.nan
    lam = barycentric(P, X)
    F = field(P, q, s, X)
    r = np.linalg.norm(X[:, None, :] - P, axis=-1)
    good = np.isfinite(X).all(axis=1) & (lam.min(axis=1) > 0) & (r.min(axis=1) > 1e-9)
    X = X[good]
    Fn = np.linalg.norm(field(P, q, s, X), axis=1) * r[good].min(axis=1) ** (s + 1)
    X = X[Fn < 1e-9]
    out = []
    for x in X:
        if all(np.linalg.norm(x - y) > 1e-7 for y in out):
            out.append(x)
    out = np.array(out)
    if len(out) == 0:
        return out, np.array([])
    H = hessV(P, q, s, out)
    ev = np.linalg.eigvalsh(H)
    kind = np.where(ev[:, 0] > 0, 0, np.where(ev[:, 1] < 0, 2, 1))  # 0=min,1=saddle,2=max
    return out, kind
