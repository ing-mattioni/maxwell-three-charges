"""Exact grid-coverage checks and high-precision arithmetic regressions.

Run directly with Python; requires NumPy and mpmath, not pytest/unittest.
The grid reference uses exact rational arithmetic, independent of interval.py.
"""
import struct
import subprocess
import sys
from fractions import Fraction
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import mpmath as mp
from domains import simplex_boxes
from interval import I, PI_LO, PI_HI, TRIG_ARG_LO, TRIG_ARG_HI, sin_i, cos_i, sincos_point
from interval import arithmetic_cases, check_arithmetic
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


def test_trig_domain_guard():
    # sin_i/cos_i add only the extrema at 0, pi/2, pi, 3pi/2 and 2pi. Outside
    # (-pi/2, 5pi/2) the hull of the endpoint values is not an enclosure, so
    # such arguments must be refused instead of being evaluated.
    for low, high in ((-2., -1.), (7., 8.), (6., 9.5), (PI_LO, 13.)):
        for extension in (sin_i, cos_i):
            try:
                extension(I(np.array([low]), np.array([high])))
            except ArithmeticError:
                continue
            raise AssertionError('out-of-range trigonometric argument accepted')
    # The whole admissible range contains all handled extrema of both functions.
    whole = I(np.array([TRIG_ARG_LO]), np.array([TRIG_ARG_HI]))
    for extension in (sin_i, cos_i):
        bounds = extension(whole)
        assert bounds.lo[0] == -1 and bounds.hi[0] == 1
    assert mp.mpf(TRIG_ARG_LO) > -mp.pi/2 and mp.mpf(TRIG_ARG_HI) < 5*mp.pi/2


def test_wide_angle_boxes():
    # Regression: this box touches theta3 = pi, tau = 1 and is wide in the angles.
    # 2*theta2 then starts just below 0, is reduced with k = -1, and the reduced
    # interval [2*pi, 2*pi + 6.8] spans more than a period. The prefilter must
    # not treat sin and cos as monotone there.
    mp.mp.dps = 50
    low = np.array([[2.6179938779914944, -1., .25, 0.]])
    high = np.array([[PI_HI, 1., .5, .25]])
    z2 = zeta2_naive(low, high)
    assert z2.lo[0] <= z2.hi[0]

    def exact(point):
        theta3, tau, l1, l2 = [mp.mpf(float(x)) for x in point]
        theta2 = (2*mp.pi-theta3)/2-tau*(3*theta3-2*mp.pi)/2
        return abs(l1+l2*mp.exp(2j*theta3)+(1-l1-l2)*mp.exp(-2j*theta2))**2

    # A point of the closed minimum region inside the box (theta2 is close to pi/2).
    witness = exact((3.1, 0.013786943456877778, .25, .25))
    assert witness < mp.mpf(1)/9
    assert mp.mpf(float(z2.lo[0])) <= witness <= mp.mpf(float(z2.hi[0]))
    rng = np.random.default_rng(20261109)
    for point in low + (high-low)*rng.uniform(size=(200, 4)):
        assert mp.mpf(float(z2.lo[0])) <= exact(point) <= mp.mpf(float(z2.hi[0]))
    # The same holds on the single box of the coarsest admissible grid.
    low, high = initial_boxes(1)
    z2 = zeta2_naive(low, high)
    for point in low[0] + (high[0]-low[0])*rng.uniform(size=(200, 4)):
        assert mp.mpf(float(z2.lo[0])) <= exact(point) <= mp.mpf(float(z2.hi[0]))


def test_environment_check():
    # The check passes here, and every expected result in it is the correctly rounded value,
    # computed from exact rationals (CPython's int/int division rounds correctly, ties to even).
    assert check_arithmetic() is True
    word = lambda x: struct.unpack('<Q', struct.pack('<d', x))[0]
    exact = {np.add: lambda a, b: Fraction(a)+Fraction(b), np.subtract: lambda a, b: Fraction(a)-Fraction(b),
             np.multiply: lambda a, b: Fraction(a)*Fraction(b), np.divide: lambda a, b: Fraction(a)/Fraction(b)}
    for name, op, args, expected in arithmetic_cases():
        if op in exact:
            q = exact[op](*args)
            assert word(float(q) if q >= 0 else -float(-q)) == word(expected), name
        elif op is np.sqrt:
            x = Fraction(args[0])
            assert all(abs(Fraction(expected)**2-x) < abs(Fraction(c)**2-x)
                       for c in (np.nextafter(expected, 0), np.nextafter(expected, 10))), name
        else:
            a, b = args
            step = 1 if (b > a) == (a >= 0) else -1
            assert word(expected) == (word(a)+step if a != 0 else (1 if b > 0 else 2**63+1)), name
    # It must refuse to load in a wrong environment. The rounding mode and the FTZ/DAZ bits
    # can be set portably only through the C runtime, so this part runs on Windows only.
    if sys.platform != 'win32':
        return
    for value, mask in ((0x100, 0x300), (0x200, 0x300), (0x300, 0x300),
                        (0x01000000, 0x03000000), (0x02000000, 0x03000000), (0x03000000, 0x03000000)):
        code = ('import ctypes, sys\n'
                f'ctypes.cdll.ucrtbase._controlfp(ctypes.c_uint({value}), ctypes.c_uint({mask}))\n'
                f'sys.path.insert(0, {str(ROOT)!r})\n'
                'try:\n    import interval\nexcept ArithmeticError:\n    sys.exit(3)\nsys.exit(0)\n')
        result = subprocess.run([sys.executable, '-c', code], capture_output=True)
        assert result.returncode == 3, f'wrong environment {value:#x}/{mask:#x} accepted'


if __name__ == '__main__':
    test_simplex_cover()
    test_angular_cover()
    test_trigonometry()
    test_negative_range_reduction()
    test_trig_domain_guard()
    test_wide_angle_boxes()
    test_environment_check()
    print('all exact-coverage and high-precision arithmetic tests passed')
