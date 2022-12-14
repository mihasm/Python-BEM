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

import os
import re as regex
import sys
from math import isinf
from subprocess import Popen, PIPE
import threading
import time

import matplotlib.cm as cm
import numpy as np
import scipy.interpolate
import scipy.linalg
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

from bem import data_path

xfoil_path = os.path.join(data_path,"xfoil_executables", "xfoil.exe")
if sys.platform.startswith("darwin"):
    xfoil_path = os.path.join(data_path,"xfoil_executables", "xfoil")


def xfoil_runner(airfoil, reynolds, alpha, ncrit, printer=None, print_all=False):
    """

    :param airfoil:
    :param reynolds:
    :param alpha:
    :param printer:
    :param print_all:
    :return:
    """
    # alpha = round(degrees(alpha))
    # alpha in degrees
    if print_all:
        if printer != None:
            printer.print("        xfoil runner, alpha=",
                          alpha, "re=", reynolds)
    out = run_xfoil_analysis(airfoil, reynolds, alpha, ncrit)
    # if out == False:
    #    return xfoil_runner(airfoil, reynolds * 0.5, radians(alpha), printer, print_all=print_all)
    return out


def get_coefficients_from_output(output_str):
    """

    :param output_str:
    :return:
    """
    lines = output_str.splitlines()
    a, CL, Cm, CD, CDf, CDp = [None] * 6
    cL_last = None
    cD_last = None

    for l in lines:
        # print(l)
        cL_find = regex.match(r".+CL =([- \d][- \d][ \d]\.\d\d\d\d)", l)
        cD_find = regex.match(r".+CD =([- \d][- \d][ \d]\.\d\d\d\d)", l)
        if cL_find:
            cL = float(cL_find.groups()[0])
            if not isinf(cL):
                cL_last = cL
        if cD_find:
            cD = float(cD_find.groups()[0])
            if not isinf(cD):
                cD_last = cD
        if "VISCAL:  Convergence failed" in l:
            return False

    # Weird cases
    out = {"CL": cL_last,
           "CD": cD_last,
           "out": output_str}
    return out


def run_xfoil_analysis(airfoil, reynolds, alpha, ncrit, iterations=100, max_next_angle=2., print_output=False):
    """

    :param airfoil:
    :param reynolds:
    :param alpha:
    :param iterations:
    :param max_next_angle:
    :param print_output:
    :return:
    """
    # print("running xfoil for %s,Re=%s,alpha=%s" % (airfoil,reynolds,alpha))
    # alpha in degrees

    with Popen(os.path.abspath(xfoil_path), stdin=PIPE, stdout=PIPE, universal_newlines=True, shell=False) as process:

        def call(_str, proc=process):
            """

            :param _str:
            :param proc:
            """
            # print(_str,file=process.stdin)
            proc.stdin.write(_str + "\n")

        def kill():
            """

            """
            time.sleep(2)
            process.kill()

        thread = threading.Thread(target=kill)
        thread.start()

        # disables graphical preview
        call("plop")
        call("G F")

        # go back
        call("")

        if "naca" in airfoil.lower():
            # treat as standard naca airfoil
            airfoil = airfoil.replace(".dat", "").strip()
            call(airfoil)
        else:
            # open .dat file
            call("load %s" % os.path.join("foils", airfoil))

        # create back-design from dat file
        call("mdes")
        call("filt")
        call("exec")
        call("")

        # calculation mode
        call("pane")
        call("oper")

        # input reynolds
        call("re")
        call("%s" % reynolds)

        # viscid mode
        call("v")

        # set ncrit
        call("vpar")
        call("n")
        call("%s" % ncrit)
        call("")

        # num iterations
        call("iteration")
        call("%s" % iterations)

        # alfa
        call("alfa")
        call("%s" % alpha)

        # quit
        call("")
        call("")
        call("quit")

        output = process.communicate()[0]
        if print_output:
            for l in output.splitlines():
                print(l)  # print(output)
    return get_coefficients_from_output(output)


def draw_to_matplotlib(x, y, z, shrani=False, unit='CL'):
    """

    :param x:
    :param y:
    :param z:
    :param shrani:
    :param unit:
    """
    # Ne vem ali je to potrebno. Menda je.
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)

    # x/y tocke za risati ozadje
    xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(
        y.min(), y.max(), 100)
    xi, yi = np.meshgrid(xi, yi)

    # Linearna interpolacija ozadja, glede na x,y,z
    # interpolacija je lahko linear, cubic, itd.
    rbf = scipy.interpolate.Rbf(x, y, z, function='quintic')
    zi = rbf(xi, yi)

    # Barvno ozadje
    # plt.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower',
    #           extent=[x.min(), x.max(), y.min(), y.max()], cmap=cm.jet)

    plt.scatter(xi, yi, c=zi)
    plt.scatter(x, y, c=z, cmap=cm.jet)

    # ODKOMENTIRAJ ZA ABSOLUTNO SKALO
    # plt.clim(0,10)

    # nastavitev min/max vrednosti na osi
    # plt.xlim(-1000,0)
    # plt.ylim(-1000,0)

    # X Label
    # plt.xlabel('alfa' % (z.min(),z.max(),np.mean(z)))

    # Y Label
    plt.ylabel('re')

    # Title
    # plt.title('Hitrost')

    # Color bar
    cbar = plt.colorbar()
    cbar.ax.set_xlabel(unit)

    # Za shranjevanje
    if shrani:
        filename = "hitrosti.png"
        path = os.path.join(os.pwd(), filename)
        print(filename, "saved to", path)
        plt.savefig(path)

    plt.show()


def generate_polars(foil, alpha_from, alpha_to, alpha_num, reynolds_from, reynolds_to, reynolds_num, ncrit):
    """

    :param foil:
    :return:
    """
    all_ncrit = []
    all_re = []
    all_a = []
    all_cl = []
    all_cd = []

    for Re in np.linspace(reynolds_from, reynolds_to, reynolds_num):
        for a in np.linspace(alpha_from, alpha_to, alpha_num):
            Re = int(Re)
            cl = None
            cd = None
            print("Re", Re, "Alpha:", a)
            o = None
            o = xfoil_runner(foil, Re, a, ncrit)
            if o != False:
                cl = o["CL"]
                cd = o["CD"]
            if cl != None and cd != None:
                print("CL:", o["CL"], "CD", o["CD"])
                all_re.append(Re)
                all_a.append(a)
                all_cl.append(cl)
                all_cd.append(cd)
                all_ncrit.append(ncrit)
    return all_ncrit, all_a, all_re, all_cl, all_cd


def generate_polars_data(foil, alpha_from, alpha_to, alpha_num, reynolds_from, reynolds_to, reynolds_num, ncrit):
    """

    :param foil:
    :return:
    """
    all_ncrit, all_a, all_re, all_cl, all_cd = generate_polars(foil, alpha_from, alpha_to, alpha_num, reynolds_from,
                                                               reynolds_to, reynolds_num, ncrit)
    out = []
    for i in range(len(all_a)):
        out.append([all_re[i], all_ncrit[i], all_a[i], all_cl[i], all_cd[i]])
    out = np.array(out)
    return out
