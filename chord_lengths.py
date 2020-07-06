import numpy as np

#R=0.45
#radiuses = np.linspace(0,0.45,10)
#Cl_max = 1.4
#B=3
#TSR = 7

def generate_chord_lengths_betz(radiuses,R,Cl_max,B,TSR):
	"""
	Source: http://wflportal.amcplaza.com/Research/DYNAM/Resource%20Documents/WT_Theory_2009.pdf
	"""
	chords = 16*np.pi*R/(9*B*Cl_max)*(TSR*np.sqrt(TSR**2*(radiuses/R)**2+4/9))**-1
	return chords

def generate_chord_lengths_schmitz(radiuses,R,Cl_max,B,TSR):
	"""
	Source: http://wflportal.amcplaza.com/Research/DYNAM/Resource%20Documents/WT_Theory_2009.pdf
	"""
	chords = 16*np.pi*radiuses/(B*Cl_max)*np.sin(1/3*np.arctan(R/(TSR*radiuses)))**2
	return chords

def generate_twists_betz(radiuses,R,TSR,alpha_d):
	"""
	Source: http://wflportal.amcplaza.com/Research/DYNAM/Resource%20Documents/WT_Theory_2009.pdf
	"""
	thetas = np.rad2deg(np.arctan(2*R/(3*radiuses*TSR)))+alpha_d
	return thetas

def generate_twists_schmitz(radiuses,R,TSR,alpha_d):
	"""
	Source: http://wflportal.amcplaza.com/Research/DYNAM/Resource%20Documents/WT_Theory_2009.pdf
	"""
	thetas = 2/3*np.rad2deg(np.arctan(R/(radiuses*TSR)))-alpha_d
	return thetas

#c = generate_chord_lengths_schmitz(radiuses=radiuses,R=R,Cl_max=Cl_max,B=B,TSR=TSR)

#c = generate_chord_lengths_betz(radiuses=radiuses,R=R,Cl_max=Cl_max,B=B,TSR=TSR)

#thetas = generate_twists_betz(radiuses=radiuses,R=R,TSR=TSR,alpha_d=3)

#thetas = generate_twists_schmitz(radiuses=radiuses,R=R,TSR=TSR,alpha_d=7)

#print(thetas)