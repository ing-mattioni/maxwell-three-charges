"""Exact grid-coverage checks and high-precision arithmetic regressions.

Run directly with Python; requires NumPy and mpmath, not pytest/unittest.
The grid reference uses exact rational arithmetic, independent of interval.py.
"""
import sys
from fractions import Fraction
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import mpmath as mp
from domains import simplex_boxes
from interval import I, PI_LO, PI_HI, sin_i, cos_i, sincos_point
from certify import initial_boxes, zeta2_naive


def test_simplex_cover():
    for n in (1, 3, 7, 10, 12, 13, 48, 64):
        lo, hi = simplex_boxes(n)
        assert len(lo) == n*(n+1)//2
        row = 0
        for i in range(n):
            for j in range(n-i):
                # Every exact rational cell in the combinatorial cover is
                # contained in its floating-point box, including the diagonal.
                for axis, k in enumerate((i, j)):
                    assert Fraction(float(lo[row, axis])) <= Fraction(k, n)
                    assert Fraction(float(hi[row, axis])) >= Fraction(k+1, n)
                row += 1
        assert (lo >= 0).all() and (hi <= 1).all() and (lo <= hi).all()
    # Regression: the old rounded grid omitted a nonempty simplex sliver here.
    old = np.linspace(0, 1, 13)
    a, b = Fraction(float(old[5])), Fraction(float(old[7]))
    assert a+b < 1
    delta = (1-a-b)/4
    witness = (a+delta, b+delta)
    lo, hi = simplex_boxes(12)
    assert any(all(Fraction(float(low[k])) <= witness[k] <= Fraction(float(high[k]))
                   for k in (0, 1)) for low, high in zip(lo, hi))


def test_angular_cover():
    mp.mp.dps = 80
    for n in (10, 12):
        lo, hi = initial_boxes(n)
        assert mp.mpf(float(lo[:, 0].min())) <= 2*mp.pi/3
        assert mp.mpf(float(hi[:, 0].max())) >= mp.pi
        assert lo[:, 1].min() == -1 and hi[:, 1].max() == 1


def test_trigonometry():
    mp.mp.dps = 80
    rng = np.random.default_rng(20261030)
    special = np.array([0., PI_LO/2, PI_LO, 1.5*PI_LO, 2*PI_LO])
    arguments = np.concatenate((rng.uniform(-.1, 2*PI_HI+.1, 600), special,
                                np.nextafter(special, -np.inf),
                                np.nextafter(special, np.inf)))
    sin_bounds, cos_bounds = sincos_point(arguments)
    for index, argument in enumerate(arguments):
        for fn, bounds in ((mp.sin, sin_bounds), (mp.cos, cos_bounds)):
            value = fn(mp.mpf(float(argument)))
            assert mp.mpf(float(bounds.lo[index])) <= value <= mp.mpf(float(bounds.hi[index]))
    for center in special:
        low, high = np.nextafter(center, -np.inf), np.nextafter(center, np.inf)
        for fn, extension in ((mp.sin, sin_i), (mp.cos, cos_i)):
            bounds = extension(I(low, high))
            for argument in (mp.mpf(float(low)), mp.mpf(float(high))):
                assert mp.mpf(float(bounds.lo)) <= fn(argument) <= mp.mpf(float(bounds.hi))


def test_negative_range_reduction():
    # Outward angular covers can extend slightly beyond the physical endpoint.
    # Then theta2 is negative and the 2*pi reduction uses a negative integer.
    mp.mp.dps = 80
    point = np.array([[PI_HI, 1., .25, .25]])
    z2 = zeta2_naive(point, point)
    theta3, tau, l1, l2 = [mp.mpf(float(x)) for x in point[0]]
    theta2 = (2*mp.pi-theta3)/2-tau*(3*theta3-2*mp.pi)/2
    assert theta2 < 0
    z = l1+l2*mp.exp(2j*theta3)+(1-l1-l2)*mp.exp(-2j*theta2)
    assert z2.lo[0] <= z2.hi[0]
    assert mp.mpf(float(z2.lo[0])) <= abs(z)**2 <= mp.mpf(float(z2.hi[0]))


if __name__ == '__main__':
    test_simplex_cover()
    test_angular_cover()
    test_trigonometry()
    test_negative_range_reduction()
    print('all exact-coverage and high-precision arithmetic tests passed')
