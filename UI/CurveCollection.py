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

from numpy.core._multiarray_umath import array

from UI.Curve import Curve


class CurveCollection:
    """

    """
    def __init__(self):
        self.curve_list = []

    def add(self, curve):
        """

        :param curve:
        """
        self.curve_list.append(curve)

    def get_stall_angles(self):
        """

        :return:
        """
        out_re = []
        out_min_aoa_list = []
        out_max_aoa_list = []

        for c in self.curve_list:
            Re = c.Re
            min_stable_aoa = c.min_stable_aoa
            max_stable_aoa = c.max_stable_aoa

            out_re.append(Re)
            out_min_aoa_list.append(min_stable_aoa)
            out_max_aoa_list.append(max_stable_aoa)

        return out_re, out_min_aoa_list, out_max_aoa_list

    def get_curves_sorted(self):
        """

        :return:
        """
        return sorted(self.curve_list, key=lambda c: (c.ncrit, c.Re))

    def save_curves(self):
        """

        :return:
        """
        out_list = []
        for c in self.curve_list:
            data_curve = c.save_curve()
            out_list.append(data_curve)
        return out_list

    def load_curves(self, out):
        """

        :param out:
        """
        self.curve_list = []
        for data_curve in out:
            c = Curve()
            c.load_curve(data_curve)
            self.curve_list.append(c)

    def gather_curves(self, extrapolation=True):
        """

        :param extrapolation:
        :return:
        """
        out = []
        for curve in self.get_curves_sorted():
            if extrapolation:
                alpha, cl, cd = curve.get_combined_curve()
            else:
                alpha, cl, cd = curve.get_curves()
            for i in range(len(alpha)):
                Re = curve.Re
                ncrit = curve.ncrit
                _alpha = alpha[i]
                _cl = cl[i]
                _cd = cd[i]
                out.append([Re, ncrit, _alpha, _cl, _cd])
        out = array(out)
        return out

    def get_curve(self, re_in, ncrit_in):
        """

        :param re_in:
        :param ncrit_in:
        :return:
        """
        out = []
        for curve in self.curve_list:
            re = curve.Re
            ncrit = curve.ncrit
            if ncrit == ncrit_in and re == re_in:
                out.append(curve)

        if len(out) == 0:
            return None
        if len(out) == 1:
            return out[0]
        if len(out) > 1:
            for c in out:
                print(c.Re, c.ncrit)
            raise Exception("DataError: Multiple curves have same Reynolds and Ncrit...")

    def remove_curve(self, re_in, ncrit_in):
        """

        :param re_in:
        :param ncrit_in:
        :return:
        """
        out = []
        i = 0
        for curve in self.curve_list:
            re = curve.Re
            ncrit = curve.ncrit
            if ncrit == ncrit_in and re == re_in:
                j = i
                out.append(curve)
            i += 1

        if len(out) == 0:
            return None
        if len(out) > 1:
            for c in out:
                print(c.Re, c.ncrit)
            raise Exception("DataError: Multiple curves have same Reynolds and Ncrit...")

        del self.curve_list[j]
