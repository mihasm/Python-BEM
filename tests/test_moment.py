import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.append(str(ROOT))

from bending_inertia import PointsInCircum, pi, calculate_bending_inertia_2


def test_circle_inertia():
    r = 7777777.7777777777777
    x, y = PointsInCircum(r, 10000)
    inertia, _, _, _area = calculate_bending_inertia_2(x, y)
    expected = pi / 4 * r ** 4
    error = abs(inertia - expected) / expected * 100
    assert error <= 1e-4


def test_hollow_circle_inertia():
    r1 = 0.5
    r2 = 1
    x1, y1 = PointsInCircum(r1, 10000)
    x2, y2 = PointsInCircum(r2, 10000)
    inertia_1, _, _, _area_1 = calculate_bending_inertia_2(x1, y1)
    inertia_2, _, _, _area_2 = calculate_bending_inertia_2(x2, y2)
    expected = pi / 4 * (r2 ** 4 - r1 ** 4)
    inertia = inertia_2 - inertia_1
    error = abs(inertia - expected) / expected * 100
    assert error <= 1e-4
