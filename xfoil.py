#!/usr/bin/env python3
import os
import re as regex
import sys
from math import isinf
from subprocess import Popen, PIPE
from threading import Timer

import matplotlib.cm as cm
import numpy as np
import scipy.interpolate
import scipy.linalg
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

xfoil_path = os.path.join("xfoil_executables", "xfoil.exe")
if sys.platform.startswith("darwin"):
    xfoil_path = os.path.join("xfoil_executables", "xfoil")


def xfoil_runner(airfoil, reynolds, alpha, printer=None, print_all=False):
    # alpha = round(degrees(alpha))
    # alpha in degrees
    if print_all:
        if printer != None:
            printer.print("        xfoil runner, alpha=",
                          alpha, "re=", reynolds)
    out = run_xfoil_analysis(airfoil, reynolds, alpha)
    # if out == False:
    #    return xfoil_runner(airfoil, reynolds * 0.5, radians(alpha), printer, print_all=print_all)
    return out


def get_coefficients_from_output(output_str):
    lines = output_str.splitlines()
    last_lines = "\n".join(lines[-5:])
    if "VISCAL:  Convergence failed" in output_str:
        return False
    if "TRCHEK2: N2 convergence failed." in last_lines:
        return False
    _, __, a, CL, Cm, CD, CDf, CDp = [None] * 8
    for l in lines:
        a = regex.match('(.+)a =(.+)CL =(.+)', l)
        _, a, CL = a.groups() if a != None else [_, a, CL]
        b = regex.match('(.+)Cm =(.+)CD =(.+)=>.+CDf =(.+)CDp =(.+)', l)
        __, Cm, CD, CDf, CDp = b.groups() if b != None else [
            __, Cm, CD, CDf, CDp]
    # Weird cases
    if CL == None or CD == None:
        return False
    if "*" in CL or "*" in CD:
        return False
    if isinf(float(CL)) or isinf(float(CD)):
        return False

    out = {"CL": float(CL), "CD": float(CD), "Cm": float(
        Cm), "CDf": float(CDf), "CDp": float(CDp), "out": output_str}
    return out


def run_xfoil_analysis(airfoil, reynolds, alpha, iterations=70, max_next_angle=2., print_output=False):
    # print("running xfoil for %s,Re=%s,alpha=%s" % (airfoil,reynolds,alpha))
    # alpha in degrees
    with Popen(os.path.abspath(xfoil_path), stdin=PIPE, stdout=PIPE, universal_newlines=True, shell = True) as process:

        def call(_str, proc=process):
            # print(_str,file=process.stdin)
            proc.stdin.write(_str + "\n")

        call("plop")
        call("G F")
        call("")

        if "naca" in airfoil.lower():
            # treat as standard naca airfoil
            airfoil = airfoil.replace(".dat", "").strip()
            call(airfoil)
        else:
            # open .dat file
            call("load %s" % os.path.join("foils", airfoil))
        call("mdes")
        call("filt")
        call("exec")
        call("")
        call("pane")
        call("oper")
        call("re")
        # rint("settings reynolds")
        call("%s" % reynolds)
        call("v")
        call("iteration")
        call("%s" % iterations)

        call("alfa")
        call("%s" % alpha)
    
        call("")
        call("quit")
        # print(process.stdout.read())
        timer = Timer(1, process.kill)
        try:
            timer.start()
            output = process.communicate()[0]
            if print_output:
                for l in output.splitlines():
                    print(l)  # print(output)
        finally:
            timer.cancel()  # print("")  # return False

    return get_coefficients_from_output(output)


def draw_to_matplotlib(x, y, z, shrani=False, unit='CL'):
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


def generate_polars(foil):
    all_re = []
    all_a = []
    all_cl = []
    all_cd = []

    for Re in np.linspace(10e3, 1e6, 10):
        for a in np.linspace(-10, 20, 31):
            Re = int(Re)
            cl = None
            cd = None
            print("Re", Re, "Alpha:", a)
            o = None
            o = xfoil_runner(foil, Re, a)
            if o != False:
                cl = o["CL"]
                cd = o["CD"]
            if cl != None and cd != None:
                print("CL:", o["CL"], "CD", o["CD"])
                all_re.append(Re)
                all_a.append(a)
                all_cl.append(cl)
                all_cd.append(cd)
    return all_a, all_re, all_cl, all_cd


def generate_polars_data(foil):

    all_a, all_re, all_cl, all_cd = generate_polars(foil)
    out = []
    ncrit = 0
    for i in range(len(all_a)):
        out.append([all_re[i], ncrit, all_a[i], all_cl[i], all_cd[i]])
    out = np.array(out)
    return out


# all_a,all_re,all_cl,all_cd = generate_polars()
# data = generate_polars_data("s826.dat")
# print(data)
# from s826 import all_re,all_a,all_cl,all_cd


def draw_scatter(x, y, z):
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(all_a, all_re, all_cl)
    ax.set_xlabel('alpha')
    ax.set_ylabel('reynolds')
    ax.set_zlabel('lift')
    plt.show()

# X,Y,Z = create_approximation_plane(all_a,all_re,all_cl)
# surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2,label="interp")

# print("all_a",all_a)
# print("all_re",all_re)
# print("all_cl",all_cl)
# print("all_cd",all_cd)

# draw_to_matplotlib(all_a,all_re,all_cl)
# draw_scatter(all_a,all_re,all_cl)
#
