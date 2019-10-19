from turbine_data import SET_INIT
from calculation import Calculator
from numpy import linspace
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

SET_INIT["convergence_limit"] = 1e-2
SET_INIT["relaxation_factor"] = 0.1
SET_INIT["cascade_correction"] = False
SET_INIT["max_iterations"] = 1000

calc = Calculator(SET_INIT["airfoils"])
v = 10
yaw_angles = []
max_cp_ar = []

plt.figure(1)

for yaw_angle in linspace(0,90,10):
	out_tsr = []
	out_cp = []
	for rpm in linspace(250,3000,10):
		print(rpm)
		out = calc.run_array(v=v,rpm=rpm,yaw_angle=yaw_angle,skewed_wake_correction=True,**SET_INIT)
		if out != None:
			out_tsr.append(out["TSR"])
			out_cp.append(out["cp_w"])
	if len(out_tsr) > 0:
		yaw_angles.append(yaw_angle)
		max_cp_ar.append(np.max(out_cp))
		plt.plot(out_tsr,out_cp,label="Yaw:"+str(yaw_angle)+"°")

plt.legend()
plt.xlabel("Tip-Speed Ratio lambda")
plt.ylabel("Efficiency Cp")
plt.title("Efficiency vs. TSR at different yaw angles")


plt.figure(2)
yaw_angles = np.array(yaw_angles)
theor_cp_ar = np.cos(np.radians(yaw_angles))**2*np.max(max_cp_ar)
plt.plot(yaw_angles,theor_cp_ar,label="cos^2")
plt.plot(yaw_angles,max_cp_ar,label="calculated")
plt.xlabel("Yaw angle [°]")
plt.ylabel("Efficiency Cp")
plt.title("Efficiency vs. Yaw angle")
plt.legend()
plt.show()