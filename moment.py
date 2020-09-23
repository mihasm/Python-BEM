from polars import get_x_y_from_link
from matplotlib import pyplot as plt
import numpy as np
from scipy.interpolate import interp1d


def airfoil_bending_inertia(x,y):
	cross = np.where(np.diff(np.signbit(np.gradient(x))))[0][0]

	x_up = x[:cross+2]
	y_up = y[:cross+2]
	x_down = x[cross+1:]
	y_down = y[cross+1:]

	interp_up = interp1d(x_up,y_up,'cubic')
	interp_down = interp1d(x_down,y_down,'cubic')

	_x=np.arange(0,1,0.0001)
	_y_up=interp_up(_x,)
	_y_down=interp_down(_x)
	return calculate_bending_inertia(_x,_y_up,_y_down)

def calculate_bending_inertia(x,y_up,y_down):
	"""
	Function that calculates the bending inertia (Second moment of Area)
	for a given airfoil, defined by its top and bottom line.

	Calculation result is given in p.u.^4 -> (per unit of length)^4

	Sources:https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/systems-labs-06/spl10b.pdf
			http://www.wingbike.nl/Wingbike_Hydrofoil/Background_files/Method%20to%20determine%20Section%20Modulus%20and%20Bending%20Inertia%20equations.pdf

	Example data:
	#x_up = [0,0.75,1.5,2.63,3.75,5.63,7.5,11.25,15,22.5,30,37.5,45,60,75,90,105,120,135,150]
	#y_up = [5.98,8.38,9.49,10.51,11.28,12.39,13.5,15.13,16.41,18.26,19.42,19.83,20,19.49,17.98,15.64,12.56,8.92,4.79,0.21]

	#x_down=[0,0.75,1.5,2.63,3.75,5.63,7.5,11.25,15,22.5,30,37.5,45,60,75,90,105,120,135,150]
	#y_down = [5.98,4.79,3.93,3.16,2.51,2.05,1.5,0.72,0.26,0.05,0,0,0,0,0,0,0,0,0,0]
	"""

	A=0
	zsum=0
	I=0

	for i in range(len(x)-1):
		dx = (x[i+1]-x[i])
		dy_up = (y_up[i+1]+y_up[i])
		dy_down = (y_down[i+1]+y_down[i])
		
		A_up=dy_up*dx*0.5
		A_down=dy_down*dx*0.5
		A+=A_up-A_down
		
		zsum+=0.5*(y_up[i+1]**2-y_down[i+1]**2)*dx

	z=zsum/A
	I=0
	for i in range(len(x)-1):
		dx = (x[i+1]-x[i])
		I+=1/3*( (y_up[i+1]-z)**3 - (y_down[i+1]-z)**3 )*dx

	return I,A

x,y = get_x_y_from_link('http://airfoiltools.com/airfoil/details?airfoil=n2414-il')

print(airfoil_bending_inertia(x,y))