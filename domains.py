"""Outward-rounded covers of the exact barycentric simplex.

For the exact rational grid k/n, the cells i+j<n cover the closed simplex.
Rounded endpoints from linspace do not necessarily retain that property.
Enclose each rational endpoint before assembling the same combinatorial cover.
"""
import numpy as np


def simplex_boxes(n):
    """Return lower/upper corners covering {(x,y)>=0: x+y<=1}."""
    if not isinstance(n, (int, np.integer)) or n < 1 or n >= 2**53:
        raise ValueError('n must be a positive integer smaller than 2**53')
    # Integers here are exactly represented. One outward step encloses the exact
    # quotient under IEEE binary64 round-to-nearest division.
    grid = np.arange(n + 1, dtype=np.float64) / float(n)
    lower = np.nextafter(grid, -np.inf)
    upper = np.nextafter(grid, np.inf)
    lower[0] = upper[0] = 0.0
    lower[-1] = upper[-1] = 1.0
    lo, hi = [], []
    for i in range(n):
        for j in range(n - i):
            lo.append([lower[i], lower[j]])
            hi.append([upper[i + 1], upper[j + 1]])
    return np.array(lo), np.array(hi)
