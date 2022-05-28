from numpy.core._multiarray_umath import array

from UI.Curve import Curve


class Curves:
    def __init__(self):
        self.curve_list = []

    def add(self, curve):
        self.curve_list.append(curve)

    def get_stall_angles(self):
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
        return sorted(self.curve_list, key=lambda c: (c.ncrit, c.Re))

    def save_curves(self):
        out_list = []
        for c in self.curve_list:
            data_curve = c.save_curve()
            out_list.append(data_curve)
        return out_list

    def load_curves(self, out):
        self.curve_list = []
        for data_curve in out:
            c = Curve()
            c.load_curve(data_curve)
            self.curve_list.append(c)

    def gather_curves(self, extrapolation=True):
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
