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

PI = np.pi

# Interpolacijski funkciji za koeficienta vzgona in upora
get_cd = interp1d(alpha_CD, CD)
get_cl = interp1d(alpha_CL, CL)

max_radius = 0.25  # radius
D = 2 * max_radius  # diameter
c = 0.05  # chord
B = 3  # num blades
blades = B
h = 0.3  # blade span [m]
rpm = 600
phi = linspace(0, 2 * pi, 360)
dphi = 2*pi/360
rho = 1.225  # kg/m^3
visc = 1.64e-5
#theta = 0
delta = 0
m_delta = 0 #blade inclination angle
m_velocity_inf = 10
CQ1=0
CQ2=0
m_bVariable = True
elements = 10
m_radius_local = np.zeros(elements)+max_radius
m_zeta = np.linspace(0.001,h,elements)/h #non dimensionles height
m_eta = m_radius_local / max_radius #non dimensionless radius
#m_theta = np.zeros(72)
m_theta = radians(np.linspace(0,360-5,72) )#azimuthal angle
epsilon     =   0.001
iterations  =   1000
relax       =   0.4
m_twist = np.zeros(elements)+0
m_bIsInverted = False
m_c_local = np.zeros(elements) + c

#final azi results
it_loc=np.zeros(72)
vrel_loc=np.zeros(72)
u_loc=np.zeros(72)
V_loc=np.zeros(72)
Re_loc=np.zeros(72)
DeltaRe_loc=np.zeros(72)
alpha_loc=np.zeros(72)
Ftip_loc=np.zeros(72)
CD_loc=np.zeros(72)
CL_loc=np.zeros(72)
LD_loc=np.zeros(72)
Cn_loc=np.zeros(72)
Ct_loc=np.zeros(72)


m_it_up=[]
m_Ftip_up=[]
m_u_up=[]
m_a_up=[]
m_velocity_up=[]
m_lambda_up=[]

m_it_dw=[]
m_Ftip_dw=[]
m_u_dw=[]
m_a_dw=[]
m_velocity_dw=[]
m_lambda_dw=[]

m_iterations=[]
m_vrel=[]
m_u=[]
m_V=[]
m_Re=[]
m_DeltaRe=[]
m_alpha=[]
m_Ftip=[]
m_CL=[]
m_CD=[]
m_LD=[]
m_Cn=[]
m_Ct=[]

m_iterations_down=[]
m_vrel_down=[]
m_u_down=[]
m_V_down=[]
m_Re_down=[]
m_DeltaRe_down=[]
m_alpha_down=[]
m_Ftip_down=[]
m_CL_down=[]
m_CD_down=[]
m_LD_down=[]
m_Cn_down=[]
m_Ct_down=[]

m_velocity_equil = []

lambda_global = 2
omega = lambda_global*m_velocity_inf/max_radius

if m_bVariable:

    #loop over all elements
    for i in range(0,elements):
        zeta = m_zeta[i]
        eta = m_eta[i]
        m_delta = 0

        #loop over all upwind azi angles
        for l in  range(0,36):
            theta = m_theta[l]
            u = 0
            delta_u = 10000
            u_old = u

            count = 0
            save = 0

            while delta_u > epsilon or u<=0 or save==0:
                if delta_u <= epsilon and u>0:
                    save=1

                count += 1

                if count == iterations:
                    save = 1
                elif count==iterations+1:
                    break

                #local induced velocity and local tip speed ratio
                V = (1-u) * m_velocity_inf
                X = m_radius_local[i] * omega/V

                #upwind function
                f = 0

                #tip loss correction
                F = 1

                #local relative velocity
                W = V * sqrt(pow((X-sin(theta)),2) + pow(F,2) * pow(cos(theta),2)*pow(cos(delta),2))

                #local angle of attack
                if (F*cos(theta)*cos(delta) * V/W >= 1):
                    alpha = asin(1)
                else:
                    alpha = asin(F*cos(theta)*cos(delta)*V/W)

                alpha_deg = asin(F * cos(theta) * cos(delta) * V/W)*180/PI+m_twist[i]*cos(delta)

                if ((X-sin(theta)<0)):
                    alpha = PI-alpha
                    alpha_deg = 180.0 - alpha_deg
                    alpha_deg = 180.0-asin(F * cos(theta) * cos(delta) * V/W)*180/PI+m_twist[i]*cos(delta)

                alpha_corrected = alpha_deg

                if (alpha_corrected > 180):
                    alpha_corrected-=360

                if (m_bIsInverted):
                    alpha_corrected *= -1

                RE = W*m_c_local[i]/visc
                CL = get_cl(alpha)
                CD = get_cd(alpha)

                if (m_bIsInverted):
                    CL *= -1

                #tip loss
                CL = CL
                CD = CD

                Cn = CL * cos(alpha) + CD * sin(alpha)
                Ct = CL * sin(alpha) - CD * cos(alpha)

                if save == 0:
                    FX = blades*m_c_local[i]/8/PI/m_radius_local[i]*pow(W/m_velocity_inf,2)*((Cn*cos(theta)+Ct*sin(theta)/cos(delta)))
                    u_old = u
                    u = pow(u,2)+FX

                    CTT = 4*u*(1-u)

                    if (CTT>0.96*F):
                        u = (18*F-20-3*pow(fabs(CTT*(50-36*F)+12*F*(3*F-4)),0.5))/(36*F-50)
                    

                    if (u<=0):
                        u=0.01
                    if (u>=1):
                        u=0.99

                    delta_u = fabs(u_old-u)

            #convergence
            u = (1-u)

            #save results
            it_loc[l] = count-1
            vrel_loc[l] = W
            u_loc[l] = u
            V_loc[l] = V
            Re_loc[l]= RE
            #DeltaRe_loc[l]= RE - REBlade
            alpha_loc[l] = alpha
            Ftip_loc[l] = F
            CD_loc[l]=CD
            CL_loc[l]=CL
            LD_loc[l]=CL/CD
            Cn_loc[l]=Cn
            Ct_loc[l]=Ct

        #average upwind interference and tiploss factor and local tipspeed ratio
        u = 0
        F = 0
        count = 0
        for l in range(0,36):
            u=u+u_loc[l]
            F = F+Ftip_loc[l]
            count=count+it_loc[l]
        u = u/36
        F = F/36
        count=count/36

        #//// save upwind dat
        m_it_up.append(count)
        #// average height dat
        m_Ftip_up.append(F)
        m_u_up.append(u)
        m_a_up.append(1-u)
        m_velocity_up.append(V)
        m_lambda_up.append(X)
        #//local azimuthal dat
        m_iterations.append(it_loc)
        m_vrel.append(vrel_loc)
        m_u.append(u_loc)
        m_V.append(V_loc)
        m_Re.append(Re_loc)
        m_DeltaRe.append(DeltaRe_loc)
        m_alpha.append(alpha_loc)
        m_Ftip.append(Ftip_loc)
        m_CL.append(CL_loc)
        m_CD.append(CD_loc)
        m_LD.append(LD_loc)
        m_Cn.append(Cn_loc)
        m_Ct.append(Ct_loc)

    #initialization for downwind part
    for i in range(0,elements):
        m_u_down.append(m_u_up[i])
        m_velocity_equil.append((2*m_u_up[i]-1)*m_velocity_inf)

    stop = 0

    #loop over all elements
    for i in range(0,elements):
        zeta = m_zeta[i]
        eta = m_eta[i]
        m_delta = 0

        #loop over all upwind azi angles
        for l in  range(36,72):
            theta = m_theta[l]
            u2 = 0
            delta_u = 10000
            u2_old = u2

            count = 0
            save = 0

            while delta_u > epsilon or u2<=0 or save==0:
                if delta_u <= epsilon and u2>0:
                    save=1

                count += 1

                if count == iterations:
                    save = 1
                elif count==iterations+1:
                    break

                #local induced velocity and local tip speed ratio
                V = (1-u2) * m_velocity_inf*(2*m_u[i][36-(l-35)]-1)
                if V <= 0:
                    V = 0.01
                X = m_radius_local[i] * omega/V

                #upwind function
                f = 0

                #tip loss correction
                F = 1

                #local relative velocity
                W = V * sqrt(pow((X-sin(theta)),2) + pow(F,2) * pow(cos(theta),2)*pow(cos(delta),2))

                #local angle of attack
                alpha = asin(F * cos(theta) * cos(delta) * V/W)
                alpha_deg = asin(F * cos(theta) * cos(delta) * V/W)*180/PI+m_twist[i]*cos(delta)
                
                if (X-sin(theta)) < 0:
                
                    alpha = -PI-alpha
                    alpha_deg = -180.0-asin(F * cos(theta) * cos(delta) * V/W)*180/PI+m_twist[i]*cos(delta)
                
                alpha_corrected = alpha_deg

                if (alpha_corrected > 180):
                    alpha_corrected -= 360

                if (m_bIsInverted):
                    alpha_corrected *= -1

                if (m_bIsInverted):
                    alpha_corrected *= -1

                RE = W*m_c_local[i]/visc
                CL = get_cl(alpha)
                CD = get_cd(alpha)

                if (m_bIsInverted):
                    CL *= -1

                #tip loss
                CL = CL
                CD = CD

                Cn = CL * cos(alpha) + CD * sin(alpha)
                Ct = CL * sin(alpha) - CD * cos(alpha)

                if save == 0:
                    FX = blades*m_c_local[i]/8/PI/m_radius_local[i]*pow(W/m_velocity_inf,2)*((Cn*cos(theta)+Ct*sin(theta)/cos(delta)))
                    u2_old = u2
                    u2 = pow(u2,2)+FX

                    """
                    CTT = 4*u*(1-u)

                    if (CTT>0.96*F):
                        u = (18*F-20-3*pow(fabs(CTT*(50-36*F)+12*F*(3*F-4)),0.5))/(36*F-50)
                    """

                    if (u2>=1):
                        u2=0.99
                    if (u2<=0):
                        u=0.01

                    delta_u = fabs(u2_old-u2)

            #convergence
            u2 = (1-u2)

            #save results
            it_loc[l] = count-1
            vrel_loc[l] = W
            u_loc[l] = u2
            V_loc[l] = V
            Re_loc[l]= RE
            alpha_loc[l] = alpha
            Ftip_loc[l] = F
            CD_loc[l]=CD
            CL_loc[l]=CL
            LD_loc[l]=CL/CD
            Cn_loc[l]=Cn
            Ct_loc[l]=Ct

        #average downwind interference and tiploss factor and local tipspeed ratio
        u2 = 0
        F = 0
        count = 0
        for l in range(36,72):
            u2=u2+u_loc[l]
            F = F+Ftip_loc[l]
            count=count+it_loc[l]
        u2 = u2/36
        F = F/36
        count=count/36

        #//// save downwind dat
        m_it_dw.append(count)
        #// average height dat
        m_Ftip_dw.append(F)
        m_u_dw.append(u)
        m_a_dw.append(1-u)
        m_velocity_dw.append(V)
        m_lambda_dw.append(X)

        #//local azimuthal dat
        m_iterations.append(it_loc)
        m_vrel.append(vrel_loc)
        m_u.append(u_loc)
        m_V.append(V_loc)
        m_Re.append(Re_loc)
        m_DeltaRe.append(DeltaRe_loc)
        m_alpha.append(alpha_loc)
        m_Ftip.append(Ftip_loc)
        m_CL.append(CL_loc)
        m_CD.append(CD_loc)
        m_LD.append(LD_loc)
        m_Cn.append(Cn_loc)
        m_Ct.append(Ct_loc)

from matplotlib import pyplot as plt
plt.plot(m_theta,m_alpha[0])
plt.show()