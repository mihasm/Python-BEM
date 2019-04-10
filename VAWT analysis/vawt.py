from numpy import *
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
import numpy as np

alpha_CL = [-8.774961369489471, -8.389547038327525, -7.980598899156693, -7.625662778366913, -7.330525677927586, -7.122193607029237, -6.861778518406301, -6.601363429783365, -6.3930313588850165, -6.149977276170277, -5.924284199363732, -5.750674140281776, -5.542342069383427, -5.316648992576882, -5.125677927586729, -4.934706862596576, -4.761096803514619, -4.604847750340857, -4.431237691258898, -4.257627632176941, -4.10137857900318, -3.945129525829419, -3.7888804726556575, -3.649992425390092, -3.493743372216329, -3.302772307226176, -3.077079230419633, -2.8340251477048923, -2.566169627978443, -2.313194970459021, -2.0354188759278884, -1.7923647932131477, -1.5195489860843594, -1.2715346159672762, -1.0284805332525355, -0.7854264505377966, -0.5597333737312518, -0.35140130283290283, -0.11082736381933245, 0.13470686259657683, 0.36039993940312165, 0.5860930162096665, 0.8117860930162113, 1.089562187547342, 1.3326162702620827, 1.5756703529768235, 1.8013634297833683, 2.041937368796937, 2.304832601121042, 2.5478866838357845, 2.7959010539528677, 3.047015603696412, 3.3291319497045926, 3.5895470383275274, 3.8326011210422664, 4.075655203757009, 4.318709286471748, 4.561763369186487, 4.822178457809425, 5.08259354643236, 5.345488778756465, 5.603423723678233, 5.881199818209364, 6.158975912740495, 6.463758016462155, 6.818694137251935, 7.200636267232241, 7.582578397212544, 7.9784093319194085, 8.346462657173156, 8.735349189516743, 9.110346917133768, 9.537427662475384, 9.915897591274053, 10.290895318891078, 10.676309650053025, 11.058637580164621, 11.436721708831996, 11.818663838812302, 12.200605968792608, 12.582548098772914, 12.976064232692021, 13.346432358733527, 13.745735494622028, 14.179760642326922, 14.551286168762314, 14.943644902287534, 15.32558703226784, 15.72681916881281, 16.089471292228453, 16.457524617482203, 16.86377215573398, 17.23529768216937, 17.606823208604762, 17.999181942129983, 18.38112407211029, 18.797788213906983, 19.0825087108014,90]
CL = [-0.317362802, -0.324851546, -0.317590619, -0.298646422, -0.278325572, -0.250585418, -0.224587401, -0.194100034, -0.165555819, -0.137011603, -0.113171153, -0.088566843, -0.058759101, -0.032148464, -0.005185571, 0.023465853, 0.049128845, 0.074931207, 0.105292604, 0.134801694, 0.170622675, 0.199858384, 0.230746441, 0.253310443, 0.279532928, 0.307500899, 0.334886882, 0.364407459, 0.395421295, 0.422185086, 0.451992828, 0.478699187, 0.507588001, 0.534811257, 0.561632481, 0.588396273, 0.612642583, 0.636831461, 0.66380584, 0.690100596, 0.715380701, 0.740589014, 0.766180215, 0.797107901, 0.824618322, 0.852530774, 0.876777085, 0.904794831, 0.931271458, 0.959298776, 0.986522032, 1.013199674, 1.043804299, 1.070718853, 1.095931952, 1.121719383, 1.147334514, 1.172203015, 1.199490883, 1.224352206, 1.251575462, 1.2758122, 1.303533209, 1.327674225, 1.353653099, 1.374315066, 1.386156713, 1.40043613, 1.415431894, 1.431111111, 1.447634594, 1.464443153, 1.481807246, 1.495717526, 1.507314299, 1.513004868, 1.521032091, 1.527050376, 1.528695049, 1.527488955, 1.529243273, 1.529608756, 1.5275986, 1.514572789, 1.50643714, 1.5031149, 1.496349811, 1.490904116, 1.482571105, 1.473872612, 1.464761123, 1.453544452, 1.441527374, 1.429020549, 1.413860318, 1.398656228, 1.379468375, 1.367911806,0]

alpha_CD = [-8.985305818050573, -8.899611358927679, -8.810489121439867, -8.70765577049239, -8.585398342143728, -8.488277955137779, -8.385444604190305, -8.255189026323503, -8.145500118646197, -8.025527875874143, -7.9055556331020895, -7.785583390330036, -7.665611147557982, -7.549066683150844, -7.374249986540137, -7.151444392820609, -6.894361015451922, -6.700120241440026, -6.637277638083236, -6.55158317896034, -6.423041490275997, -6.088833099696704, -5.725488593015627, -5.338149637780139, -4.957666239274483, -4.5806106191337435, -4.169277215343845, -3.802504930297852, -3.4151659750623633, -3.038110354921624, -2.6610547347808833, -2.301138006464722, -1.9754990617977182, -1.392776739762029, -1.015721119621288, -0.638665499480549, -0.26160987933980806, 0.10516240570618507, 0.49250136094167196, 0.9038347647315703, 1.2808903848723112, 1.6579460050130521, 1.9493071660308967, 2.628007282284228, 2.9947795673302213, 3.3718351874709622, 3.748890807611703, 4.125946427752444, 4.503002047893185, 4.880057668033922, 5.211409576642451, 5.634168908315404, 6.011224528456145, 6.388280148596886, 6.765335768737627, 7.142391388878364, 7.519447009019105, 7.896502629159846, 8.273558249300587, 8.650613869441328, 9.027669489582069, 9.404725109722806, 9.816058513512706, 10.193114133653447, 10.570169753794188, 10.93008648211035, 11.27286431860193, 11.598503263268933, 11.907003316111359, 12.181225585304622, 12.455447854497889, 12.746809015515735, 13.003892392884419, 13.226697986603947, 13.415225796674319, 13.569475823095532, 13.723725849516741, 13.895114767762534, 14.066503686008323, 14.237892604254116, 14.426420414324484, 14.597809332570277, 14.75205935899149, 14.966295506798728, 15.111976087307651, 15.249087221904285, 15.386198356500914, 15.540448382922127, 15.71183730116792, 15.86608732758913, 16.003198462185765, 16.140309596782394, 16.285990177291318, 16.405962420063375, 16.56878189239687, 16.740170810642667, 16.911559728888456, 17.048670863485086, 17.168643106257143, 17.305754240853773, 17.42572648362583, 17.5285598345733, 17.651959855710274, 17.775359876847247, 17.871337671064886, 17.991309913836943, 18.111282156608993, 18.21411550755647, 18.3512266421531, 18.45405999310058, 18.53975445222347, 18.676865586820107, 18.779698937767577, 18.899671180539634, 19.01621564494677, 19.088198990610003, 19.187604563192565,90]
CD = [0.0818635, 0.078751752, 0.074613041, 0.07123231, 0.067382034, 0.064386331, 0.061196773, 0.05807161, 0.055700627, 0.052734885, 0.050130343, 0.047192902, 0.044347454, 0.040770856, 0.038074992, 0.035152787, 0.032329762, 0.028848415, 0.025676968, 0.023347541, 0.020220509, 0.017500649, 0.015474406, 0.013606955, 0.012038062, 0.010636742, 0.008796366, 0.00752969, 0.007632868, 0.007940207, 0.008192664, 0.008330601, 0.008336926, 0.008697579, 0.008895154, 0.008993941, 0.008993941, 0.009135537, 0.009268351, 0.009531785, 0.009564714, 0.009612278, 0.009509832, 0.009996818, 0.010146463, 0.010296474, 0.010431849, 0.010596495, 0.010684307, 0.01085627, 0.010958716, 0.011134339, 0.011430702, 0.01177097, 0.012517365, 0.013863803, 0.01539684, 0.016864018, 0.018334856, 0.019684952, 0.020903333, 0.022304653, 0.024009653, 0.025535372, 0.027068409, 0.029375648, 0.031830702, 0.034485202, 0.037266702, 0.040327375, 0.043259002, 0.046240059, 0.049262162, 0.052144599, 0.054908212, 0.057081539, 0.060019554, 0.063126607, 0.065903635, 0.068450452, 0.071388315, 0.074025438, 0.076657578, 0.080063462, 0.083016571, 0.085485711, 0.087981012, 0.090972556, 0.09418908, 0.097276009, 0.099981264, 0.102753597, 0.105731859, 0.108755543, 0.111945964, 0.115914298, 0.119741767, 0.122842495, 0.125677163, 0.12898243, 0.132249128, 0.134430504, 0.137856848, 0.140300633, 0.143551769, 0.147015005, 0.150454093, 0.153056719, 0.156688321, 0.159892771, 0.162967827, 0.166985463, 0.170652952, 0.173857602, 0.17708137, 0.180276965, 0.183355845,0.9]

get_cd = interp1d(alpha_CD,CD)
get_cl = interp1d(alpha_CL,CL)

R = 4 #radius
D = 2*R #diameter
c=0.05 #chord
B=2 #num blades
h = 2 #blade span [m]
U=5 #wind speed
phi = linspace(0,2*pi,360)
rho = 1.225 #kg/m^3
din_visc = 1.82e-5
theta = 0
delta = 0

max_iterations = 10


CP_all = []
TSR_all = []

for rpm in list(linspace(40,200,10)):
    #rpm = 1500
    omega = 2*pi*rpm/60
    TSR = omega*R/U


    alpha_all = []
    Q_all = []
    CL_all = []
    CD_all = []
    U_rel_all = []
    for p in phi:
        max_reached = False

        u = 0.3

        i = 0
        while True:
            i+=1
            if i>max_iterations:
                #max_reached = True
                break

            if 0 <= degrees(p) < 180:
                upwind = True
            else:
                upwind = False
            
            if upwind:
                ui = (1-u)*U
                U_rel  = sqrt((ui*sin(p))**2+(ui*cos(p)+omega*R)**2)
                alpha = arctan2((1-u)*sin(p),(1-u)*cos(p)+TSR)
            else:
                uiprime = (1-u)*(1-2*u)*U
                U_rel = sqrt((uiprime*sin(p))**2+(uiprime*cos(p)+omega*R)**2)
                alpha = arctan2(-(1-u)*sin(p),(1-u)*cos(p)+TSR)

            CL = get_cl(degrees(alpha))
            CD = get_cd(degrees(alpha))

            if upwind:
                Cn = CL*cos(alpha)+CD*sin(alpha)
                Ct = CL*sin(alpha)-CD*cos(alpha)
                Q = 0.5*rho*U_rel**2*h*c*Ct*R
                #CT = -Ct*cos(p)/sin(p)+Cn
            else:
                Cn = CL*cos(alpha)+CD*sin(alpha)
                Ct = CL*sin(alpha)-CD*cos(alpha)
                Q = 0.5*rho*U_rel**2*h*c*Ct*R
                #CT = Ct*cos(p)/sin(p)+Cn


            # Original QBlade implementation is in
            # https://github.com/rmaugusto/qblade/blob/4540751660ffa790c9b8fb5c8c071233954c8b29/src/XDMS/DData.cpp
            
            F = 1 #tip loss factor

            Fx = B*c/8/pi/R*(U_rel/U)**2*((Cn*cos(theta)+Ct*sin(theta)/cos(delta)))

            u_new = u**2+Fx

            CTT = 4*u_new*(1-u_new)

            if CTT>0.96*F:
                u_new = (
                        18 * F - 20 - 3 * sqrt(CTT * (50 - 36 * F) +
                                               12 * F * (3 * F - 4))
                ) / (36 * F - 50)

            if u_new <=0:
                u_new=0.01

            if u_new >= 1:
                u_new=0.99

            if abs(u_new-u) <= 0.001:
                break
            else:
                u = u_new
            
        CL_all.append(CL)
        CD_all.append(CD)
        U_rel_all.append(U_rel)

        if not max_reached:
            Q_all.append(Q)
            alpha_all.append(alpha)
        else:
            Q_all.append(0)
            alpha_all.append(0)

    plt.plot(phi,Q_all)
    
    Qa = mean(Q_all)
    P = Qa * omega
    CQ = Qa/(0.5*rho*U**2*D*h*R)
    CP = CQ*TSR
    TSR_all.append(TSR)
    CP_all.append(CP)

TSR_All = np.array(TSR_all)
plt.legend(TSR_All.astype(int))
plt.figure(2)

plt.plot(TSR_all,CP_all)
plt.show()
