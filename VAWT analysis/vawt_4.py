from numpy import *
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
from numpy import arcsin as asin

alpha_CL,alpha_CD,CL,CD = [],[],[],[]
f = open("naca0015.0_360_M.dat")
for line in f.readlines():
    _alpha,_cl,_cd = line.split("\t")
    alpha_CL.append(float(_alpha))
    alpha_CD.append(float(_alpha))
    CL.append(float(_cl))
    CD.append(float(_cd))


rho = 1.225  # kg/m^3
visc = 1.64e-5
delta = 0


# Interpolacijski funkciji za koeficienta vzgona in upora
get_cd = interp1d(alpha_CD, CD)
get_cl = interp1d(alpha_CL, CL)

R = 0.25  # radius [m]
c = 0.05  # chord [m]
B = 3  # num B
h = 0.3  # blade height [m]
D = 2 * R  # diameter

elements = 1 # number of blade elements
m_radius_local = np.zeros(elements)+R
m_zeta = np.linspace(0.001,h,elements)/h #non dimensionles height
m_eta = m_radius_local / R #non dimensionless radius
m_theta = np.linspace(0,360-5,72) #azimuthal angle
m_theta = radians(m_theta)
m_twist = np.zeros(elements)
m_c_local = np.zeros(elements) + c
dh = h/elements #height difference



m_velocity_inf = 10 # does it matter??
lambda_global = 5

CP_all = []
TSR_all = []

for lambda_global in np.linspace(1,10,10):

    omega = lambda_global*m_velocity_inf/R
    epsilon     =   0.001
    iterations  =   100

    # Arrays of arrays
    m_u = []
    m_W = []
    m_alpha = []
    m_V = []
    m_Q = []

    # Arrays per element
    u_loc=np.zeros(72)
    W_loc=np.zeros(72)
    alpha_loc = np.zeros(72)
    V_loc=np.zeros(72)
    Q_loc=np.zeros(72)

    #loop over all elements for upwind part (0 - 180)
    for i in range(0,elements):
        zeta = m_zeta[i] #non dimensionless height
        eta = m_eta[i] #non dimensionless radius

        #loop over all upwind azi angles
        for l in  range(0,36):

            theta = m_theta[l]

            u = 0
            delta_u = 10000
            u_old = u

            count = 0

            while True:
                count += 1

                if delta_u <= epsilon:
                    break

                if count >= iterations:
                    max_iterations = True
                    break

                #local induced velocity and local tip speed ratio
                V = (1-u) * m_velocity_inf
                X = m_radius_local[i] * omega/V

                #upwind function
                #f = 0

                #tip loss correction
                F = 1.

                #local relative velocity
                W = sqrt(V*sin(theta)**2+(V*cos(theta)+omega*R))
                #W = V * sqrt(pow((X-sin(theta)),2) + pow(F,2) * pow(cos(theta),2)*pow(cos(delta),2))

                #alpha = arctan((1-u)*sin(theta)/((1-u)*cos(theta)+lambda_global))

                #alpha = asin(F * cos(theta) * cos(delta) * V/W);
                alpha = arctan((1-u)*sin(theta)/((1-u)*cos(theta)+lambda_global))

                RE = W*m_c_local[i]/visc
                CL = get_cl(degrees(alpha))
                CD = get_cd(degrees(alpha))

                Cn = CL * cos(alpha) + CD * sin(alpha)
                Ct = CL * sin(alpha) - CD * cos(alpha)

                Q = 0.5 * rho * W ** 2 * dh * m_c_local[i] * Ct * m_radius_local[i]

                FX = B*m_c_local[i]/8/pi/m_radius_local[i]*pow(W/m_velocity_inf,2)*((Cn*cos(theta)+Ct*sin(theta)/cos(delta)))
                u_old = u
                u = pow(u,2)+FX

                CTT = 4*u*(1-u)

                if (CTT>0.96*F):
                    u = (18*F-20-3*pow(fabs(CTT*(50-36*F)+12*F*(3*F-4)),0.5))/(36*F-50)

                #if (u<=0):
                #    u=0.01
                #if (u>=1):
                #    u=0.99

                delta_u = fabs(u_old-u)

            #convergence
            u = (1-u)

            #save results
            u_loc[l] = u
            W_loc[l] = W
            alpha_loc[l] = alpha
            V_loc[l] = V
            Q_loc[l] = Q


        #save upwind
        m_u.append(u_loc)
        m_W.append(W_loc)
        m_alpha.append(alpha_loc)
        m_V.append(V_loc)
        m_Q.append(Q_loc)



    #loop over all elements for downwind part (180 - 360)
    for i in range(0,elements):
        zeta = m_zeta[i]
        eta = m_eta[i]

        #loop over all downwind azi angles
        for l in  range(36,72):
            theta = m_theta[l]

            u2 = 0
            delta_u = 10000
            u2_old = u2

            count = 0

            while count <= iterations:
                count += 1

                if delta_u <= epsilon:
                    break

                #local induced velocity and local tip speed ratio
                #V = (1-u2) * m_velocity_inf*(2*m_u[i][36-(l-35)]-1)
                V = (1-u2) * m_velocity_inf

                #if V <= 0:
                #    V = 0.01
                X = m_radius_local[i] * omega/V

                #upwind function
                f = 0

                #tip loss correction
                F = 1

                #local relative velocity
                W = sqrt(V*sin(theta)**2+(V*cos(theta)+omega*R))
                #W = V * sqrt(pow((X-sin(theta)),2) + pow(F,2) * pow(cos(theta),2)*pow(cos(delta),2))
                #alpha = arctan((1-u2)*sin(theta)/((1-u2)*cos(theta)+lambda_global))
                #alpha = asin(F * cos(theta) * cos(delta) * V/W)
                alpha = arctan((1-u2)*sin(theta)/((1-u2)*cos(theta)+lambda_global))

                RE = W*m_c_local[i]/visc
                CL = get_cl(degrees(alpha))
                CD = get_cd(degrees(alpha))

                Cn = CL * cos(alpha) + CD * sin(alpha)
                Ct = CL * sin(alpha) - CD * cos(alpha)

                Q = 0.5 * rho * W ** 2 * dh * m_c_local[i] * Ct * m_radius_local[i]

                FX = B*m_c_local[i]/8/pi/m_radius_local[i]*pow(W/m_velocity_inf,2)*((Cn*cos(theta)+Ct*sin(theta)/cos(delta)))
                u2_old = u2
                u2 = pow(u2,2)+FX

                #if (u2>=1):
                #    u2=0.99
                #if (u2<=0):
                #    u2=0.01

                delta_u = fabs(u2_old-u2)

            #convergence
            u2 = (1-u2)

           #save results
            u_loc[l] = u
            W_loc[l] = W
            alpha_loc[l] = alpha
            V_loc[l] = V
            Q_loc[l] = Q


        #save upwind
        m_u.append(u_loc)
        m_W.append(W_loc)
        m_alpha.append(alpha_loc)
        m_V.append(V_loc)
        m_Q.append(Q_loc)



    plt.plot(degrees(m_theta),degrees(m_alpha[0]))

    Qa = mean(m_Q)
    Qtot = Qa * B
    Pkin =  0.5*rho*m_velocity_inf**3*h*D
    P = Qtot*omega
    CP = P/Pkin

    #https://www.academia.edu/30009924/Numerical_and_Analytical_Investigation_of_Vertical_Axis_Wind_Turbine

    #Qavg = B*np.sum((0.5*rho*np.array(U_rel_all)**2*(h*c)*np.array(Ct_all)*R)/360)
    #CQ = Qavg / (0.5 * rho * U ** 2 * D * h * R)
    #CP = CQ * TSR

    TSR_all.append(lambda_global)
    CP_all.append(CP)

plt.figure(2)
plt.plot(TSR_all,CP_all)

plt.show()