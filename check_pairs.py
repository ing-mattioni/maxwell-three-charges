"""For suspicious shapes, directly count minima for charges q(x), x on a grid of M."""
import sys
import numpy as np
from core import triangle, charges_at, detT, critical_points
from scan_injectivity import grid_in_triangle


def max_minima(a, b, s, n=60, nstart=40):
    P = triangle(a, b)
    X, lam = grid_in_triangle(P, n)
    inM = detT(P, X) > (s + 1) / (s + 2) ** 2
    Q = charges_at(P, X[inM], s)
    best = (0, None)
    counts = {}
    for q in Q:
        C, kind = critical_points(P, q, s, n=nstart)
        nm = int((kind == 0).sum())
        counts[len(C)] = counts.get(len(C), 0) + 1
        if nm > best[0]:
            best = (nm, q)
        if nm >= 2:
            print('TWO MINIMA', a, b, s, q, C, kind)
    return best, counts


if __name__ == '__main__':
    for (a, b, s) in [(0.4860, 0.8377, 1.0), (0.4628, 0.8088, 1.0), (0.4860, 0.8377, 3.0), (0.3874, 0.7875, 3.0)]:
        best, counts = max_minima(a, b, s)
        print(a, b, s, 'max #minima =', best[0], 'critical-point count histogram', counts)
