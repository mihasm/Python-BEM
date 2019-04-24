import numpy as np
#method from  Xiang, S., Liu, Y., Tong, G., Zhao, W., Tong, S., & Li, Y. (2018). An improved propeller design method for the electric aircraft. Aerospace Science and Technology, 78, 488–493.
"""
T = 245 # Thrust [N]
c1 = 0.1875 # max chord [m]
beta = 100 # coefficient
m = 0.001
n = 3 # blade number
#r = 0.5 #m current section radius
R = 0.8 #m tip radius
all_r = np.linspace(0,R,10)
epsilon_i = all_r/R
rho = 1.25 #kg/m^3 air density
Omega = 1600 #rpm
V = 100/3.6 #km/s
Lambda = V/(Omega*R)
k_p = 2/np.pi*np.arccos(np.e**(-n/2*(1-all_r/R)*np.sqrt(1+(Omega*R/V)**2)))

#determining Lagrange factor https://www.researchgate.net/publication/263656781_Aerodynamic_Performance_of_Propellers_with_Parametric_Considerations_on_the_Optimal_Design
F1 = np.sum((Omega*all_r)**2*k_p*all_r/(V**2+(Omega*all_r)**2))
K = T/(F1*rho*4*np.pi*V**2)

b = c1*beta**(-epsilon_i-m)**2
delta = np.arctan(Lambda/epsilon_i*(1+K))
"""

def chord_distribution(B,Cl,rpm,V,R):
	out = []
	for _r in np.linspace(0,1,11):
		lambda_r = 2*np.pi*R*rpm/60/V*(_r/R)
		phi = np.arctan(2/(3*lambda_r))
		_c = 8*np.pi*_r*np.sin(phi)/(3*B*Cl*lambda_r)
		out.append((_c*R,_r*R))
	[print(str(r)+'\t'+str(c)) for c,r in out]
	return

chord_distribution(3,1.2,1500,10,1.5)