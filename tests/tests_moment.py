from bending_inertia import PointsInCircum, calculate_bending_inertia, pi, calculate_bending_inertia_2


def test_circle():
    r = 1
    x, y = PointsInCircum(r, 10000)
    I, A = calculate_bending_inertia(x, y, 100000)
    I_theor = pi / 4 * r ** 4
    error = abs(I - I_theor) / I_theor * 100
    print("Ix:", I, "m4", "A:", A, "m2")
    print("Ix_theor", I_theor)
    print("Napaka:", abs(I - I_theor) / I_theor * 100, " %")

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
