"""Projective geometry of the caustic C = F(dM) in the linear charge chart.

For x0 in dM, the kernel direction k of the Hessian is the eigenvector of T(x0) for its
smaller eigenvalue.  The tangent line of C at F(x0) is  { q : <q, m(x0)> = 0 },
     m_i(x0) = r_i^{-s-1} (u_i . k).
We track along dM:
  * the tangent-line direction theta (mod pi) of C  -> inflections = sign changes of theta',
    total turning Theta;
  * cusps = points where k is tangent to dM;
  * n(q) = number of tangent lines through q = #zeros of <q, m(x0)> on dM.
"""
import sys
import numpy as np
from core import triangle, barycentric, charges_at, Tmat
from inspect_case import trace_boundary


def dual_data(P, B, s):
    T = Tmat(P, B)
    w, V = np.linalg.eigh(T)
    k = V[..., :, 0]  # eigenvector of smaller eigenvalue
    d = B[:, None, :] - P
    r = np.linalg.norm(d, axis=-1)
    u = d / r[..., None]
    m = r ** (-s - 1) * np.einsum('nij,nj->ni', u, k)
    return k, m, w


def analyse(a, b, s, npts=6000, verbose=False):
    P = triangle(a, b)
    c, B = trace_boundary(P, s, m=npts)
    k, m, w = dual_data(P, B, s)
    # orient m continuously (k defined up to sign)
    for i in range(1, len(m)):
        if np.dot(m[i], m[i - 1]) < 0:
            m[i] = -m[i]
            k[i] = -k[i]
    # tangent line direction in chart q3 = 1 - q1 - q2:  (m2-m3, m3-m1)
    dvec = np.stack([m[:, 1] - m[:, 2], m[:, 2] - m[:, 0]], -1)
    th = np.unwrap(2 * np.arctan2(dvec[:, 1], dvec[:, 0])) / 2  # mod pi lift
    dth = np.diff(np.concatenate([th, th[:1] + (th[-1] - th[0] + (th[1] - th[0])) * 0]))
    dth = np.diff(th)
    # total turning of the tangent LINE over the closed curve
    Theta = (th[-1] - th[0]) + ((th[0] - th[-1]) % np.pi if False else 0)
    # close the loop properly: compare last to first mod pi
    gap = (th[0] - th[-1] + np.pi / 2) % np.pi - np.pi / 2
    Theta = th[-1] - th[0] + gap
    sgn = np.sign(dth)
    infl = np.count_nonzero(sgn[1:] * sgn[:-1] < 0)
    # cusps: k tangent to dM  <=>  cross(k, dB) changes sign
    dB = np.roll(B, -1, axis=0) - np.roll(B, 1, axis=0)
    cr = k[:, 0] * dB[:, 1] - k[:, 1] * dB[:, 0]
    # k was made continuous along the curve except possibly at the closing point
    cusps = np.count_nonzero(np.sign(cr[1:]) * np.sign(cr[:-1]) < 0)
    if np.sign(cr[0]) * np.sign(cr[-1]) < 0 and np.dot(m[0], m[-1]) > 0:
        cusps += 1
    return dict(Theta_over_pi=Theta / np.pi, inflections=infl, cusps=cusps,
                min_abs_dtheta=np.abs(dth).min(), B=B, m=m, P=P)


if __name__ == '__main__':
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    N = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    fixed_s = float(sys.argv[3]) if len(sys.argv) > 3 else None
    for j in range(N):
        while True:
            a, b = rng.uniform(0, 0.5), 10 ** rng.uniform(-2, 0)
            if (a - 1) ** 2 + b ** 2 <= 1:
                break
        s = fixed_s if fixed_s else 10 ** rng.uniform(-1, 1.5)
        try:
            r = analyse(a, b, s, npts=3000)
        except Exception as e:
            print('fail', a, b, s, e)
            continue
        flag = '  ' if (r['inflections'] == 0 and r['cusps'] == 3 and abs(abs(r['Theta_over_pi']) - 1) < 1e-6) else '!!'
        print(f"{flag} a={a:.4f} b={b:.4f} s={s:.3f}  Theta/pi={r['Theta_over_pi']:+.6f}  cusps={r['cusps']}  inflections={r['inflections']}  min|dtheta|={r['min_abs_dtheta']:.2e}")
