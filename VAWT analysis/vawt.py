from numpy import *
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
import numpy as np

alpha_CL,alpha_CD,CL,CD = [],[],[],[]

# Pridobitev podatkov o profilu
f = open("NACA_0012_T1_Re0.100_M0.00_N9.0_360_M.dat")
#f = open("NACA_6508_T1_Re0.100_M0.00_N9.0_360_M.dat")
for line in f.readlines()[14:]:
    _alpha,_cl,_cd = [l for l in line.strip().split(" ") if not l == ""]
    alpha_CL.append(float(_alpha))
    alpha_CD.append(float(_alpha))
    CL.append(float(_cl))
    CD.append(float(_cd))


# Interpolacijski funkciji za koeficienta vzgona in upora
get_cd = interp1d(alpha_CD, CD)
get_cl = interp1d(alpha_CL, CL)

R = 0.25  # radius
D = 2 * R  # diameter
c = 0.05  # chord
B = 3  # num blades
h = 0.3  # blade span [m]
rpm = 600
phi = linspace(0, 2 * pi, 360)
dphi = 2*pi/360
rho = 1.225  # kg/m^3
din_visc = 1.82e-5
theta = 0
delta = 0
U_min = 3
U_max = 20

max_iterations = 10

CP_all = []
TSR_all = []
U_all = []

print("TSR_min",2 * pi * rpm / 60*R/U_max)
print("TSR_max",2 * pi * rpm / 60*R/U_min)

for U in list(linspace(U_min, U_max, 10)):
    omega = 2 * pi * rpm / 60
    TSR = omega * R / U

    alpha_all = []
    Q_all = []
    CL_all = []
    CD_all = []
    U_rel_all = []
    phi_all = []
    Ct_all = []

    for p in phi:
        max_reached = False

        u=0.0
        F=1.5

        i = 0 #current interation
        while True:
            i += 1
            if i > max_iterations:
                max_reached = True
                break

            ui = (1 - u) * U
            U_rel = sqrt((ui * sin(p)) ** 2 + (ui * cos(p) + omega * R) ** 2)
            alpha = arctan2((1 - u) * sin(p), (1 - u) * cos(p) + TSR)

            CL = get_cl(degrees(alpha))
            CD = get_cd(degrees(alpha))

            
            Cn = CL * cos(alpha) + CD * sin(alpha)
            Ct = CL * sin(alpha) - CD * cos(alpha)
            
            Q = 0.5 * rho * U_rel ** 2 * h * c * Ct * R
            
            L = 0.5*CL*rho*U_rel**2*c*h
            D = 0.5*CD*rho*U_rel**2*c*h

            FT = L*sin(alpha)-D*cos(alpha)
            FN = L*cos(alpha)+D*sin(alpha)

            Fx = B * c / 8 / pi / R * (U_rel / U) ** 2 * (Cn * cos(theta) + Ct * sin(theta) / cos(delta))
            u_new = u ** 2 + Fx
            
            if u_new < 0.7:
                u_new=0.7

            if abs(u_new - u) <= 0.001:
                #print(u)
                break
            else:
                u = u_new

        if not max_reached:
            phi_all.append(p)
            U_all.append(U)
            CL_all.append(CL)
            CD_all.append(CD)
            U_rel_all.append(U_rel)
            Q_all.append(Q)
            alpha_all.append(alpha)
            Ct_all.append(Ct)

    plt.figure(1)
    plt.plot(degrees(phi_all),degrees(alpha_all),label = round(U,2))
    plt.figure(2)
    plt.plot(phi_all, U_rel_all, label=round(U,2))

    Qa = mean(Q_all)
    Qtot = Qa
    P =  Qtot * omega
    CQ = Qtot / (0.5 * rho * U ** 2 * D * h * R)
    CP = CQ * TSR

    TSR_all.append(TSR)
    CP_all.append(CP)

TSR_All = np.array(TSR_all)
U_all = np.array(U_all)
U_rel_all = np.array(U_rel_all)

plt.figure(1)
plt.legend(title="Hitrost vetra [m/s]")
plt.title("Napadni kot v odvisnosti od azimutnega kota")
plt.xlabel("Azimutni kot [°]")
plt.ylabel("Napadni kot [°]")

plt.figure(2)
plt.legend(title="Hitrost vetra [m/s]")
plt.title("Relativna hitrost vetra v odvisnosti od azimutnega kota")
plt.xlabel("Azimutni kot [m/s]")
plt.ylabel("Napadni kot [°]")

plt.figure(3)
plt.title("Izkoristek v odvisnosti od TSR")
plt.xlabel("lambda")
plt.ylabel("Cp")

plt.plot(TSR_all, CP_all)
plt.show()
