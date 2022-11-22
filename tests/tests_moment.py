# Python BEM - Blade Element Momentum Theory Software.

# Copyright (C) 2022 M. Smrekar

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import sys 
sys.path.append('..')
from bending_inertia import PointsInCircum, pi, calculate_bending_inertia_2


def test_circle():
    r = 7777777.7777777777777
    x, y = PointsInCircum(r, 10000)
    I, _, _, A = calculate_bending_inertia_2(x, y)
    I_theor = pi / 4 * r ** 4
    error = abs(I - I_theor) / I_theor * 100
    print("I:", I, "m4", "A:", A, "m2")
    print("I_theor", I_theor)
    print("Napaka:", abs(I - I_theor) / I_theor * 100, " %")
    if error <= 1e-4:
        return True
    return False


def test_hollow_circle():
    r1 = 0.5
    r2 = 1
    x1, y1 = PointsInCircum(r1, 10000)
    x2, y2 = PointsInCircum(r2, 10000)
    I1, _, _, A1 = calculate_bending_inertia_2(x1, y1)
    I2, _, _, A2 = calculate_bending_inertia_2(x2, y2)
    I_theor = pi / 4 * (r2 ** 4 - r1 ** 4)
    Ix = I2 - I1
    error = abs(Ix - I_theor) / I_theor * 100
    print("Ix", Ix)
    print("Ix_theor", I_theor)
    print("Napaka:", abs(Ix - I_theor) / I_theor * 100, " %")
    if error <= 1e-4:
        return True
    return False

test_hollow_circle()
test_circle()