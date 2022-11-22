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

from scipy.interpolate import interp1d

from montgomerie import Montgomerie


class Curve:
    """

    """
    def __init__(self):
        self.x = None
        self.y = None
        self.Re = None
        self.ncrit = None
        self.alpha = None
        self.cl = None
        self.cd = None
        self.A = None
        self.B = None
        self.Am = None
        self.Bm = None
        self.m_CD90 = None
        self.slope = None
        self.min_stable_aoa = None
        self.max_stable_aoa = None

    def create(self, x, y, Re, ncrit, alpha, cl, cd):
        """

        :param x:
        :param y:
        :param Re:
        :param ncrit:
        :param alpha:
        :param cl:
        :param cd:
        """
        self.x = x
        self.y = y
        self.Re = Re
        self.ncrit = ncrit
        self.alpha = alpha
        self.cl = cl
        self.cd = cd
        self.A = -5
        self.B = 5
        self.Am = 8
        self.Bm = 5
        self.m_CD90 = 1.5
        self.slope = 0.106
        self.max_stable_aoa = 10
        self.min_stable_aoa = -5

    def get_curves(self):
        """

        :return:
        """
        return self.alpha, self.cl, self.cd

    def get_extrapolated_curve(self):
        """

        :return:
        """
        M = Montgomerie(x=self.x, y=self.y, alpha=self.alpha, Cl=self.cl, Cd=self.cd, Re=self.Re, A=self.A, Am=self.Am,
                        B=self.B, Bm=self.Bm, m_CD90=self.m_CD90, slope=self.slope)
        alpha, cl, cd = M.calculate_extrapolation()
        return alpha, cl, cd

    def get_combined_curve(self):
        """

        :return:
        """
        M = Montgomerie(x=self.x, y=self.y, alpha=self.alpha, Cl=self.cl, Cd=self.cd, Re=self.Re, A=self.A, Am=self.Am,
                        B=self.B, Bm=self.Bm, m_CD90=self.m_CD90, slope=self.slope)
        _alpha, _cl, _cd = M.calculate_extrapolation()
        cl_out, cd_out = [], []
        f_cl = interp1d(self.alpha, self.cl, bounds_error=True)
        f_cd = interp1d(self.alpha, self.cd, bounds_error=True)
        for i in range(len(_alpha)):
            a = _alpha[i]
            try:
                cl = f_cl(a)
            except ValueError:
                cl = _cl[i]
            try:
                cd = f_cd(a)
            except ValueError:
                cd = _cd[i]
            cl_out.append(cl)
            cd_out.append(cd)
        return _alpha, cl_out, cd_out

    def save_curve(self):
        """

        :return:
        """
        out = {"x": list(self.x), "y": list(self.y), "Re": self.Re, "ncrit": self.ncrit, "alpha": list(self.alpha),
               "cl": list(self.cl), "cd": list(self.cd), "A": self.A, "B": self.B, "Am": self.Am, "Bm": self.Bm,
               "m_CD90": self.m_CD90, "slope": self.slope, "min_stable_aoa": self.min_stable_aoa,
               "max_stable_aoa": self.max_stable_aoa}
        return out

    def load_curve(self, out):
        """

        :param out:
        """
        self.x = out["x"]
        self.y = out["y"]
        self.Re = out["Re"]
        self.ncrit = out["ncrit"]
        self.alpha = out["alpha"]
        self.cl = out["cl"]
        self.cd = out["cd"]
        self.A = out["A"]
        self.B = out["B"]
        self.Am = out["Am"]
        self.Bm = out["Bm"]
        self.m_CD90 = out["m_CD90"]
        self.slope = out["slope"]
        try:
            self.min_stable_aoa = out["min_stable_aoa"]
            self.max_stable_aoa = out["max_stable_aoa"]
        except:
            self.max_stable_aoa = 10
            self.min_stable_aoa = -5
