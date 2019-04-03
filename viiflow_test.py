import viiflow as vf
import viiflowtools.vf_tools as vft
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
#%matplotlib inline
#%config InlineBackend.figure_format = 'svg'
import logging
logging.getLogger().setLevel(logging.WARNING)
matplotlib.rcParams['figure.figsize'] = [11, 6]

# Read Airfoil Data
S805 = vft.repanel_spline(vft.read_selig("foils/s826.dat"),50)
plt.plot(S805[0,:],S805[1,:],color='black')
plt.axis('equal');
plt.title('The NREL S826 profile.');
#plt.show()

"""
results = {} # Dictionary of results
AOARange = np.arange(-5,18.5,.5)
# Go over RE range
for RE in [5e5, 7e5, 10e5, 15e5, 20e5]:
    
    # Settings
    ncrit = 10.2
    Mach = 0.04*RE/5e5 # c = 0.5m, assuming 20°C
    s = vf.setup(Re=RE,Ma=Mach,ncrit=ncrit,alpha=AOARange[0],tolerance = 1e-3)
    s.silent = True

    # (Maximum) internal iterations
    s.itermax = 100

    results[RE] = {} # Sub-Dictionary of results
    results[RE]["AOA"] = []
    results[RE]["CL"] = []
    results[RE]["CD"] = []
    results[RE]["TRUP"] = []
    results[RE]["TRLO"] = []
    
    # Go over AOA range
    faults = 0
    init = True
    for alpha in AOARange:
        
        # Set current alpha and set res/grad to None to tell viiflow that they are not valid
        s.alpha = alpha
        res = None
        grad = None
        
        # Set-up and initialize based on inviscid panel solution
        # This calculates panel operator
        if init:
            [p,bl,x] = vf.init([S805],s)
            init = False

        # Run viiflow
        [x,flag,res,grad,_] = vf.iter(x,bl,p,s,res,grad)
        # If converged add to cl/cd vectors (could check flag as well, but this allows custom tolerance 
        # to use the results anyways)
        if flag:
            results[RE]["AOA"].append(alpha)
            results[RE]["CL"].append(p.CL)
            results[RE]["CD"].append(bl[0].CD)
            # Calculate transition position based on BL variable
            results[RE]["TRUP"].append( \
                vft.interpLinear(p.foils[0].S,p.foils[0].X[0,:],bl[0].ST-bl[0].bl_fl.node_tr_up.xi[0]))
            results[RE]["TRLO"].append( \
                vft.interpLinear(p.foils[0].S,p.foils[0].X[0,:],bl[0].ST+bl[0].bl_fl.node_tr_lo.xi[0]))
            faults = 0
        else:
            faults+=1
            init = True
            
        # Skip current polar if 4 unconverged results in a row
        if faults>3:
            break
"""
results = {}
AoA = 0
RE = 100e5
ncrit = 10
#Mach = 0.1
Mach = 0.04*RE/5e5 #20°C
#Mach=0
#print(dir(vf))
s = vf.setup(Re=RE,Ma=Mach,ncrit=ncrit,alpha=AoA)
s.silent = True
s.itermax = 100

# Set current alpha and set res/grad to None to tell viiflow that they are not valid
s.alpha = AoA
res = None
grad = None

# Set-up and initialize based on inviscid panel solution
# This calculates panel operator
[p,bl,x] = vf.init([S805],s)


# Run viiflow
[x,flag,res,grad,_] = vf.iter(x,bl,p,s,res,grad)

if flag == True:
	print("AoA",s.alpha)
	print("CL",p.CL)
	print("CD",bl[0].CD)
else:
	print("no convergence")
	print("AoA",s.alpha)
	print("CL",p.CL)
	print("CD",bl[0].CD)

