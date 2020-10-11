from scipy.interpolate import interp1d

from montgomerie import Montgomerie


class Curve:
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

    def create(self, x, y, Re, ncrit, alpha, cl, cd):
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

    def get_curves(self):
        return self.alpha, self.cl, self.cd

    def get_extrapolated_curve(self):
        M = Montgomerie(x=self.x, y=self.y, alpha=self.alpha, Cl=self.cl, Cd=self.cd, Re=self.Re, A=self.A, Am=self.Am,
                        B=self.B, Bm=self.Bm, m_CD90=self.m_CD90, slope=self.slope)
        alpha, cl, cd = M.calculate_extrapolation()
        return alpha, cl, cd

    def get_combined_curve(self):
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
        out = {"x": list(self.x), "y": list(self.y), "Re": self.Re, "ncrit": self.ncrit, "alpha": list(self.alpha),
               "cl": list(self.cl), "cd": list(self.cd), "A": self.A, "B": self.B, "Am": self.Am, "Bm": self.Bm,
               "m_CD90": self.m_CD90, "slope": self.slope}
        return out

    def load_curve(self, out):
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