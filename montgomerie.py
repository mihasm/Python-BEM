import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import Slider, Button
from numpy import *
from numpy import pi as PI


class POLAR_CLASS:
    def __init__(self, x, y, alpha, Cl, Cd):
        self.x = x
        self.y = y
        self.m_Alpha = alpha
        self.m_Cl = Cl
        self.m_Cd = Cd


class Montgomerie:
    def __init__(self, x, y, alpha, Cl, Cd, Re=100000, A=-5, Am=8, B=5, Bm=5, m_CD90=2.0, slope=0.106):
        # coefficient of lift at AoA  == 0
        self.CLzero = np.interp(0, alpha, Cl)
        self.CL180 = 0  # lift coefficient at AoA == 180
        # angle of attack where lift  == 0
        self.alphazero = np.interp(0, Cl, alpha)

        self.deltaCD = 0
        self.deltaalpha = 1

        self.slope = slope  # ocitno od 0 do 1
        self.m_CD90 = m_CD90  # Drag coefficientat AoA == 90

        self.m_pctrlA = A  # -10, 30
        self.m_pctrlB = B  # 1-100
        self.m_pctrlAm = Am  # 1-80
        self.m_pctrlBm = Bm  # 1-70

        if len(y) > 0:
            self.m_fThickness = np.max(y)-np.min(y)
        else:
            self.m_fThickness = 0.1

        self.m_fCamber = 0.043

        self.m_pCurPolar = POLAR_CLASS(x, y, alpha, Cl, Cd)
        self.posalphamax = np.argmax(self.m_pCurPolar.m_Cl)

        self.reynolds = Re

    def CD90(self, alpha):
        res = self.m_CD90 - 1.46 * self.m_fThickness / 2 + \
            1.46 * self.m_fCamber * sin(alpha / 360 * 2 * PI)
        return res

    def PlateFlow(self, alphazero, CLzero, alpha):
        # tukaj bi CD90  rala biti funkcija
        # alpha in degrees
        a = (1 + CLzero / sin(PI / 4) * sin(alpha / 360 * 2 * PI)) * self.CD90(alpha) * sin((alpha - 57.6 * 0.08 * sin(
            alpha / 360 * 2 * PI) - alphazero * cos(alpha / 360 * 2 * PI)) / 360 * 2 * PI) * cos(
            (alpha - 57.6 * 0.08 * sin(alpha / 360 * 2 * PI) - alphazero * cos(alpha / 360 * 2 * PI)) / 360 * 2 * PI)
        return a

    def PotFlow(self, CLzero, slope, alpha):
        return CLzero + slope * alpha

    def CDPlate(self, alpha):
        res = self.CD90(alpha) * sin(alpha / 360 * 2 * PI) ** 2
        return res

    def calculate_extrapolation(self):
        m_Alpha = []
        m_Cl = []
        m_Cd = []
        # print("A",self.m_pctrlA,"B",self.m_pctrlB)
        # positive extrapolation

        if len(self.m_pCurPolar.m_Alpha) > self.m_pctrlA + self.posalphamax >= 0:
            a1plus = self.m_pCurPolar.m_Alpha[int(
                self.posalphamax + self.m_pctrlA)]
            CL1plus = self.m_pCurPolar.m_Cl[int(
                self.posalphamax + self.m_pctrlA)]
        else:
            a1plus = (self.posalphamax + self.m_pctrlA) * self.deltaalpha
            CL1plus = self.PlateFlow(
                self.alphazero, self.CLzero, a1plus) + 0.03

        if (self.posalphamax + self.m_pctrlB + self.m_pctrlA) < len(
                self.m_pCurPolar.m_Alpha) and self.posalphamax + self.m_pctrlB + self.m_pctrlA >= 0:
            a2plus = self.m_pCurPolar.m_Alpha[int(
                self.posalphamax + self.m_pctrlB + self.m_pctrlA)]
            CL2plus = self.m_pCurPolar.m_Cl[int(
                self.posalphamax + self.m_pctrlB + self.m_pctrlA)]
        else:
            a2plus = (self.posalphamax + self.m_pctrlB +
                      self.m_pctrlA) * self.deltaalpha
            CL2plus = self.PlateFlow(
                self.alphazero, self.CLzero, a2plus) + 0.03

        A = (self.PotFlow(self.CLzero, self.slope, a1plus) -
             self.PlateFlow(self.alphazero, self.CLzero, a1plus))
        if A == 0.:
            A = 1e-5
        f1plus = (
            (CL1plus - self.PlateFlow(self.alphazero, self.CLzero, a1plus)) / A)
        B = (self.PotFlow(self.CLzero, self.slope, a2plus) -
             self.PlateFlow(self.alphazero, self.CLzero, a2plus))
        if B == 0.:
            B = 1e-5
        f2plus = (
            (CL2plus - self.PlateFlow(self.alphazero, self.CLzero, a2plus)) / B)

        if (f1plus == 1):
            f1plus += 10e-6
            print("yes")
        if (f2plus == 1):
            print("yes2")
            f2plus += 10e-6

        G = pow((fabs((1 / f1plus - 1) / (1 / f2plus - 1))), 0.25)

        am = (a1plus - G * a2plus) / (1 - G)

        k_div = pow((a2plus - am), 4)
        if k_div == 0:
            k_div = 1e-5
        k = (1 / f2plus - 1) * 1 / k_div

        # rear end flying
        self.CL180 = self.PlateFlow(self.alphazero, self.CLzero, 180)

        self.slope2 = 0.8 * self.slope
        Re = self.reynolds  # Reynolds
        deltaCL = 1.324 * pow((1 - exp(Re / 1000000 * (-0.2))), 0.7262)

        CL1plus = self.CL180 - deltaCL
        a1plus = 170 + self.CL180 / self.slope2
        a2plus = a1plus - 15
        CL2plus = self.PlateFlow(self.alphazero, self.CLzero, a2plus) - 0.01

        f1plus = (CL1plus - self.PlateFlow(self.alphazero, self.CLzero, a1plus)) / (
            self.PotFlow(self.CL180, self.slope2, a1plus - 180) - self.PlateFlow(self.alphazero, self.CLzero,
                                                                                 a1plus))
        f2plus = (CL2plus - self.PlateFlow(self.alphazero, self.CLzero, a2plus)) / (
            self.PotFlow(self.CL180, self.slope2, a2plus - 180) - self.PlateFlow(self.alphazero, self.CLzero,
                                                                                 a2plus))

        G2 = pow(fabs(((1 / f1plus - 1) / (1 / f2plus - 1))), 0.25)

        am2 = (a1plus - G2 * a2plus) / (1 - G2)

        k2 = (1 / f2plus - 1) * 1 / pow((a2plus - am2), 4)

        alpha = int(1)

        while alpha <= 180:
            if alpha < (am2 - 70):
                if alpha < am:
                    delta = 0
                else:
                    delta = am - alpha
                f = 1 / (1 + k * pow(delta, 4))
                m_Alpha.append(alpha)
                m_Cl.append(f * self.PotFlow(self.CLzero, self.slope, alpha) + (1 - f) * self.PlateFlow(self.alphazero,
                                                                                                        self.CLzero,
                                                                                                        alpha))

            elif alpha < am2:
                delta = am2 - alpha
                f = 1 / (1 + k2 * pow(delta, 4))
                m_Alpha.append(alpha)
                m_Cl.append(
                    f * self.PotFlow(self.CL180, self.slope2, alpha - 180) + (1 - f) * self.PlateFlow(self.alphazero,
                                                                                                      self.CLzero,
                                                                                                      alpha))
            else:
                m_Alpha.append(alpha)
                m_Cl.append(self.PotFlow(self.CL180, self.slope2, alpha - 180))

            if alpha < am:
                delta = 0
            else:
                delta = am - alpha

            f = 1 / (1 + k * delta ** 4)
            self.deltaCD = 0.13 * ((f - 1) * self.PotFlow(self.CLzero, self.slope, alpha) - (1 - f) * self.PlateFlow(
                self.alphazero, self.CLzero, alpha))
            if (self.deltaCD <= 0):
                self.deltaCD = 0
            # tukaj nisem preprican kaj pomeni self.m_fThickness, je to v procentih tetive ali kaj?
            m_Cd.append(
                f * (self.deltaCD + 0.006 + 1.25 * self.m_fThickness ** 2 / 180 * abs(alpha)) + (1 - f) * self.CDPlate(
                    alpha) + 0.006)

            alpha += self.deltaalpha

        # negative extrapolation
        a1minus = (-float(self.m_pctrlAm) / 20 - self.CLzero) / self.slope - 4
        CL1minus = -float(self.m_pctrlAm) / 20

        a2minus = a1minus - self.m_pctrlBm * 2
        CL2minus = self.PlateFlow(self.alphazero, self.CLzero, a2minus) - 0.03

        f1minus = (CL1minus - self.PlateFlow(self.alphazero, self.CLzero, a1minus)) / (
            self.PotFlow(self.CLzero, self.slope, a1minus) - self.PlateFlow(self.alphazero, self.CLzero, a1minus))
        f2minus = (CL2minus - self.PlateFlow(self.alphazero, self.CLzero, a2minus)) / (
            self.PotFlow(self.CLzero, self.slope, a2minus) - self.PlateFlow(self.alphazero, self.CLzero, a2minus))

        G = abs((1 / f1minus - 1) / (1 / f2minus - 1)) ** 0.25

        am = (a1minus - G * a2minus) / (1 - G)

        k = (1 / f2minus - 1) * 1 / (a2minus - am) ** 4

        # rear end flying first

        CL1minus = self.CL180 + deltaCL
        a1minus = -170 + self.CL180 / self.slope2
        a2minus = a1minus + 15
        CL2minus = self.PlateFlow(self.alphazero, self.CLzero, a2minus) - 0.01

        f1minus = (CL1minus - self.PlateFlow(self.alphazero, self.CLzero, a1minus)) / (
            self.PotFlow(self.CL180, self.slope2, a1minus + 180) - self.PlateFlow(self.alphazero, self.CLzero,
                                                                                  a1minus))
        f2minus = (CL2minus - self.PlateFlow(self.alphazero, self.CLzero, a2minus)) / (
            self.PotFlow(self.CL180, self.slope2, a2minus + 180) - self.PlateFlow(self.alphazero, self.CLzero,
                                                                                  a2minus))

        G2 = abs(((1 / f1minus - 1) / (1 / f2minus - 1))) ** 0.25

        am2 = (a1minus - G2 * a2minus) / (1 - G2)

        k2 = (1 / f2minus - 1) * 1 / (a2minus - am2) ** 4

        alpha = 0

        while alpha >= -180:
            if (alpha > am2 + 70):
                if alpha > am:
                    delta = 0
                else:
                    delta = am - alpha
                f = 1 / (1 + abs(k * delta ** 4))
                m_Alpha.append(alpha)
                m_Cl.append(f * self.PotFlow(self.CLzero, self.slope, alpha) + (1 - f) * self.PlateFlow(self.alphazero,
                                                                                                        self.CLzero,
                                                                                                        alpha))
            elif alpha > am2:
                delta = am2 - alpha
                f = 1 / (1 + abs(k2 * delta ** 4))
                m_Alpha.append(alpha)
                m_Cl.append(
                    f * self.PotFlow(self.CL180, self.slope2, alpha + 180) + (1 - f) * self.PlateFlow(self.alphazero,
                                                                                                      self.CLzero,
                                                                                                      alpha))
            else:
                m_Alpha.append(alpha)
                m_Cl.append(self.PotFlow(self.CL180, self.slope2, alpha + 180))

            if alpha > am:
                delta = 0
            else:
                delta = am - alpha
            f = 1 / (1 + k * delta ** 4)
            self.deltaCD = 0.13 * (
                self.PotFlow(self.CLzero, self.slope, alpha) - f * self.PotFlow(self.CLzero, self.slope, alpha) - (
                    1 - f) * self.PlateFlow(self.alphazero, self.CLzero, alpha))

            if (self.deltaCD <= 0):
                self.deltaCD = 0
            m_Cd.append(
                f * (self.deltaCD + 0.006 + 1.25 * self.m_fThickness ** 2 / 180 * abs(alpha)) + (1 - f) * self.CDPlate(
                    alpha) + 0.006)
            alpha = alpha - self.deltaalpha
        return m_Alpha, m_Cl, m_Cd