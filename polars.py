from scrape_polars import get_polars
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from pprint import pprint
from scipy.interpolate import Rbf
from montgomerie import Montgomerie
from xfoil import generate_polars_data
# from numpy import interp
from scipy.interpolate import interp1d
import scipy.interpolate
import matplotlib.cm as cm
from scipy.interpolate import griddata
import numpy


def scrape_data(link):
    out = []
    data = get_polars(link)
    for Re, value in data.items():
        for ncrit, value2 in value.items():
            for alpha, value3 in value2.items():
                cl = value3["cl"]
                cd = value3["cd"]
                out.append([Re, ncrit, alpha, cl, cd])
    out = np.array(out)
    return out


def get_extrapolated_data(data, airfoil_x=[], airfoil_y=[]):
    # imp_polar = np.loadtxt(open("foils/NACA_0015_polar.csv", "rb"), delimiter=",", skiprows=1)
    print("Getting inteprolation function")

    out = []

    x, y, z_cl, z_cd = [], [], [], []

    Re_list = np.unique(data[:, 0])
    ncrit_list = np.unique(data[:, 1])
    ncrit_selected = ncrit_list[0]

    for Re in Re_list:
        rows_with_Re = data[np.in1d(data[:, 0], Re)]
        rows_with_Re = rows_with_Re[np.in1d(
            rows_with_Re[:, 1], ncrit_selected)]

        _alpha = rows_with_Re[:, 2].flatten()
        _cl = rows_with_Re[:, 3].flatten()
        _cd = rows_with_Re[:, 4].flatten()

        M = Montgomerie(x=airfoil_x, y=airfoil_y,
                        alpha=_alpha, Cl=_cl, Cd=_cd, Re=Re)

        m_Alpha, m_Cl, m_Cd = M.calculate_extrapolation()
        f_cl = interp1d(_alpha, _cl, bounds_error=True)
        f_cd = interp1d(_alpha, _cd, bounds_error=True)
        for i in range(len(m_Alpha)):
            x.append(Re)
            y.append(m_Alpha[i])
            try:
                cl = f_cl(m_Alpha[i])
            except ValueError:
                cl = m_Cl[i]
            try:
                cd = f_cd(m_Alpha[i])
            except ValueError:
                cd = m_Cd[i]
            z_cl.append(cl)
            z_cd.append(cd)
            out.append([Re, ncrit_selected, m_Alpha[i], cl, cd])
    return np.array(out)
