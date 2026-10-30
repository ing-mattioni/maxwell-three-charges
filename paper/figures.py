"""Figures for the paper (SVG).  Run from the repository root:  python paper/figures.py"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from fast import trace, dual, charges, Zfun, bary  # noqa: E402
from core import critical_points, hessV  # noqa: E402

plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['cmr10'], 'mathtext.fontset': 'cm',
    'axes.formatter.use_mathtext': True, 'font.size': 9, 'svg.fonttype': 'path',
    'axes.linewidth': 0.6,
})
S = 1.0
KAPPA = S / (S + 2)
TEAL, TEAL_SOFT, RED, INK, GREY = '#1f6f6a', '#dcecea', '#a8322b', '#1a1a1a', '#9a9a9a'
OUT = os.path.dirname(os.path.abspath(__file__))
A3 = np.array([[0.0, 0.0], [1.0, 0.0], [0.5, np.sqrt(3) / 2]])  # equilateral embedding of the simplex


def chart(q, qstar):
    """projective chart of charge space centred at qstar: q -> q/qstar, normalised, equilateral embedding"""
    y = q / qstar
    y = y / y.sum(axis=-1, keepdims=True)
    return y @ A3


def oriented_dual(p, zb):
    k, m = dual(p, S, zb)
    for i in range(1, len(m)):
        if np.dot(m[i], m[i - 1]) < 0:
            m[i] = -m[i]
    return m


def tangent_params(q, m):
    """indices i where <q, m> changes sign along the loop (m is one-sided: flip at wrap)"""
    g = m @ q
    idx = list(np.nonzero(np.sign(g[1:]) * np.sign(g[:-1]) < 0)[0])
    if np.sign(g[-1]) * np.sign(-g[0]) < 0:
        idx.append(len(g) - 1)
    return idx


def tri_mask(p, X, Y):
    lam = bary(p, X + 1j * Y)
    return lam.min(-1) > 0


def classify(p, q, pts):
    P2 = np.stack([p.real, p.imag], 1)
    kinds = []
    for x in pts:
        H = hessV(P2, q, S, np.array(x))
        kinds.append('min' if np.linalg.det(H) > 0 else 'saddle')
    return kinds


def panel_plane(ax, p, q=None, show_eq=True, levels=None):
    P2 = np.stack([p.real, p.imag], 1)
    xmin, xmax = p.real.min(), p.real.max()
    ymin, ymax = p.imag.min(), p.imag.max()
    pad = 0.02 * max(xmax - xmin, ymax - ymin)
    xs = np.linspace(xmin - pad, xmax + pad, 700)
    ys = np.linspace(ymin - pad, ymax + pad, int(700 * (ymax - ymin + 2 * pad) / (xmax - xmin + 2 * pad)) + 2)
    Xg, Yg = np.meshgrid(xs, ys)
    inside = tri_mask(p, Xg, Yg)
    Zabs = np.abs(Zfun(p, (Xg + 1j * Yg).ravel())).reshape(Xg.shape)
    Zabs[~inside] = np.nan
    ax.contourf(Xg, Yg, Zabs, levels=[0, KAPPA], colors=[TEAL_SOFT])
    ax.contour(Xg, Yg, Zabs, levels=[KAPPA], colors=[TEAL], linewidths=1.0)
    if q is not None:
        V = sum(q[i] / np.abs(Xg + 1j * Yg - p[i]) ** S for i in range(3))
        V[~inside] = np.nan
        if levels is not None:
            ax.contour(Xg, Yg, V, levels=levels, colors=[GREY], linewidths=0.45)
    ax.plot(np.append(p.real, p.real[0]), np.append(p.imag, p.imag[0]), color=INK, lw=0.8)
    if q is not None:
        qq = q / q.max()
        for i in range(3):
            ax.plot(p.real[i], p.imag[i], 'o', color=RED, ms=3 + 5 * qq[i], zorder=5)
    ax.set_aspect('equal')
    ax.axis('off')


def fig1():
    p = np.array([0.0, 1.0, 0.36 + 0.78j])
    c, zb, _ = trace(p, S, 4096)
    m = oriented_dual(p, zb)
    qstar = charges(p, np.array([c]), S)[0]
    # a charge vector with 4 equilibria: image of a point of M near its centre
    x_in = c + 0.35 * (zb[700] - c)
    q_in = charges(p, np.array([x_in]), S)[0]
    P2 = np.stack([p.real, p.imag], 1)
    pts, kind = critical_points(P2, q_in, S)
    kinds = ['min' if k == 0 else 'saddle' for k in kind]
    Vc = [sum(q_in[i] / np.abs(complex(*x) - p[i]) ** S for i in range(3)) for x in pts]
    Vmin = min(Vc)
    levels = sorted(set([*(v * (1 + 1e-9) for v in Vc), *np.geomspace(Vmin * 1.0005, Vmin * 2.2, 9)]))

    fig = plt.figure(figsize=(6.6, 3.05))
    ax = fig.add_axes([0.0, 0.02, 0.47, 0.93])
    panel_plane(ax, p, q_in, levels=levels)
    for x, kd in zip(pts, kinds):
        if kd == 'min':
            ax.plot(x[0], x[1], 'o', color=INK, ms=4.2, zorder=6)
        else:
            ax.plot(x[0], x[1], 'D', mfc='white', mec=INK, ms=3.6, mew=0.9, zorder=6)
    for i, name in enumerate(['$p_1$', '$p_2$', '$p_3$']):
        off = [(-0.01, 0.035), (0.0, 0.035), (0.03, -0.01)][i]
        ax.text(p.real[i] + off[0], p.imag[i] + off[1], name, fontsize=9, clip_on=False,
                ha=['left', 'right', 'left'][i])
    ax.text(0.585, 0.205, '$M$', color=TEAL, fontsize=10, ha='left')
    ax.set_title('(a) plane of the charges', fontsize=9)

    bx = fig.add_axes([0.52, 0.06, 0.47, 0.86])
    Cpts = chart(charges(p, zb, S), qstar)
    bx.plot(np.append(Cpts[:, 0], Cpts[0, 0]), np.append(Cpts[:, 1], Cpts[0, 1]), color=RED, lw=1.1)
    lo, hi = Cpts.min(0), Cpts.max(0)
    span = (hi - lo).max()
    ctr = 0.5 * (lo + hi)
    for qv, mk, lab in ((q_in, 'o', '$q$'), (None, 's', "$q'$")):
        if qv is None:
            # a point outside the caustic
            qpt = ctr + np.array([0.58, 0.22]) * span
            # invert chart to get q' (barycentric w.r.t. A3, then multiply by qstar)
            Tm = np.vstack([A3.T, np.ones(3)])
            y = np.linalg.solve(Tm, np.append(qpt, 1.0))
            qv = y * qstar
        qpt = chart(qv, qstar)
        for i in tangent_params(qv, m):
            a = Cpts[i]
            d = a - qpt
            d /= np.linalg.norm(d)
            L = 3 * span
            bx.plot([qpt[0] - L * d[0], qpt[0] + L * d[0]], [qpt[1] - L * d[1], qpt[1] + L * d[1]], color=GREY, lw=0.5, zorder=1)
            bx.plot(*a, '.', color=INK, ms=3, zorder=4)
        bx.plot(*qpt, mk, color=INK, ms=4, zorder=5)
        bx.text(qpt[0] + 0.03 * span, qpt[1] + 0.03 * span, lab, fontsize=9)
    bx.set_xlim(ctr[0] - 0.62 * span, ctr[0] + 0.9 * span)
    bx.set_ylim(ctr[1] - 0.62 * span, ctr[1] + 0.62 * span)
    bx.set_aspect('equal')
    bx.set_xticks([]); bx.set_yticks([])
    bx.text(ctr[0] + 0.06 * span, ctr[1] - 0.40 * span, '$C=F(\\partial M)$', color=RED, fontsize=9)
    bx.set_title('(b) charge space (projective chart)', fontsize=9)
    fig.savefig(os.path.join(OUT, 'fig1.svg'), bbox_inches='tight', pad_inches=0.03)
    if os.environ.get('PREVIEW'):
        fig.savefig(os.path.join(os.environ['PREVIEW'], 'fig1.png'), dpi=160)
    plt.close(fig)
    print('fig1: equilibria', kinds, 'q_in', q_in)


def fig2():
    shapes = [('equilateral', np.array([0, 1, 0.5 + np.sqrt(3) / 2 * 1j])),
              ('right', np.array([0, 1, 0.0 + 0.8j])),
              ('obtuse', np.array([0, 1, 0.3 + 0.28j])),
              ('nearly flat', np.array([0, 1, 0.62 + 0.07j]))]
    fig, axes = plt.subplots(2, 4, figsize=(6.6, 3.4), gridspec_kw=dict(height_ratios=[1, 1.15]))
    for col, (name, p) in enumerate(shapes):
        c, zb, _ = trace(p, S, 4096)
        qstar = charges(p, np.array([c]), S)[0]
        panel_plane(axes[0, col], p)
        axes[0, col].set_title(name, fontsize=9)
        Cpts = chart(charges(p, zb, S), qstar)
        ax = axes[1, col]
        ax.plot(np.append(Cpts[:, 0], Cpts[0, 0]), np.append(Cpts[:, 1], Cpts[0, 1]), color=RED, lw=1.0)
        ax.plot(*chart(qstar, qstar), '+', color=INK, ms=5, mew=0.6)
        ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
        lo, hi = Cpts.min(0), Cpts.max(0)
        ctr, span = 0.5 * (lo + hi), (hi - lo).max()
        ax.set_xlim(ctr[0] - 0.62 * span, ctr[0] + 0.62 * span)
        ax.set_ylim(ctr[1] - 0.62 * span, ctr[1] + 0.62 * span)
    fig.tight_layout(h_pad=0.4, w_pad=0.6)
    fig.savefig(os.path.join(OUT, 'fig2.svg'), bbox_inches='tight', pad_inches=0.03)
    if os.environ.get('PREVIEW'):
        fig.savefig(os.path.join(os.environ['PREVIEW'], 'fig2.png'), dpi=160)
    plt.close(fig)


if __name__ == '__main__':
    fig1()
    fig2()
