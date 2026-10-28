"""Closed-form (division-free up to a well-conditioned Jacobian) boundary invariants.

Variables: theta1, theta2 (theta3 = 2pi - theta1 - theta2, assumed to be the largest angle),
           lambda1, lambda2 (lambda3 = 1 - lambda1 - lambda2).
Frame: v1 = 1, v2 = e^{i theta3}, v3 = e^{-i theta2}  (unit vectors from x0 to p_i).
rho_i = lambda_i / sin(theta_i) = 1/r_i.  All length-dependent quantities are rescaled by
sigma3 = sin(theta3) (so that the flat face theta3 -> pi is regular):
   rhoh_i = sigma3 * rho_i,  muh_i = sigma3 * dlambda_i/dz,  Zh_z = sigma3 Z_z, ...,  th = t/sigma3.
Returns: jac (sign of Jacobian of Z), Ih = I/(3 sigma3) (inflection function), turn = d arg t / d omega.
"""
import numpy as np


def quantities(th1, th2, l1, l2, s, xp=np):
    th3 = 2 * np.pi - th1 - th2
    l3 = 1 - l1 - l2
    s1, s2, s3 = xp.sin(th1), xp.sin(th2), xp.sin(th3)
    c1, c2, c3 = xp.cos(th1), xp.cos(th2), xp.cos(th3)
    v = [1.0 + 0j, c3 + 1j * s3, c2 - 1j * s2]
    lam = [l1, l2, l3]
    sig = [s1, s2, s3]
    rho = [l1 / s1, l2 / s2, None]        # rho3 = l3/s3 is singular at the flat face
    rhoh = [s3 * rho[0], s3 * rho[1], l3]  # sigma3 * rho_i
    zeta = l1 * v[0] ** 2 + l2 * v[1] ** 2 + l3 * v[2] ** 2

    def s3rr(a, b):
        """sigma3 * rho_a * rho_b, finite"""
        if a == 2:
            return l3 * rho[b]
        if b == 2:
            return l3 * rho[a]
        return s3 * rho[a] * rho[b]

    conj = np.conj
    muh = []
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        # mu_i = (i/2) rho_i (rho_{i+2} conj(v_{i+1}) - rho_{i+1} conj(v_{i+2}))
        muh.append(0.5j * (s3rr(i, k) * conj(v[j]) - s3rr(i, j) * conj(v[k])))
    Zz = sum(muh[i] * v[i] ** 2 - lam[i] * rhoh[i] * v[i] for i in range(3))
    Zzb = sum(conj(muh[i]) * v[i] ** 2 + lam[i] * rhoh[i] * v[i] ** 3 for i in range(3))
    jac = abs(Zz) ** 2 - abs(Zzb) ** 2
    Z = zeta
    iZ = 1j * Z
    t = (conj(Zz) * iZ - Zzb * conj(iZ)) / jac
    # K = conj(k)^2 = -conj(zeta)/|zeta|
    Kb = -conj(zeta) / abs(zeta)
    dot = lambda a, b: (a * conj(b)).real
    Q = lambda a, b: 0.5 * (a * conj(b)).real + 0.5 * (a * b * Kb).real
    v1, v3 = v[0], v[2]
    Q13 = Q(v1, v3)
    Ih = 0.5 * s1 * s2 + (s + 2) * Q13 * (l1 * s3 * dot(v1, t) - s1 * l3 * dot(v3, t)) \
        - l1 * s3 * Q(v3, t) + s1 * l3 * Q(v1, t)
    # second derivatives (scaled by sigma3^2)
    Zzz = sum(-2 * muh[i] * rhoh[i] * v[i] for i in range(3))
    Zzzb = sum(muh[i] * rhoh[i] * v[i] ** 3 - conj(muh[i]) * rhoh[i] * v[i] - lam[i] * rhoh[i] ** 2 * v[i] ** 2 for i in range(3))
    Zzbzb = sum(2 * conj(muh[i]) * rhoh[i] * v[i] ** 3 + 2 * lam[i] * rhoh[i] ** 2 * v[i] ** 4 for i in range(3))
    rhs = -Z - (Zzz * t ** 2 + 2 * Zzzb * abs(t) ** 2 + Zzbzb * conj(t) ** 2)
    tp = (conj(Zz) * rhs - Zzb * conj(rhs)) / jac
    turn = (tp / t).imag
    return dict(jac=jac, Ih=Ih, turn=turn, zeta=zeta, t=t)
