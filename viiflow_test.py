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
plt.axis('equal')
plt.title('The NREL S826 profile.')
#plt.show()

results = {} # Dictionary of results
AOARange = np.arange(-90,90,5.)

# Go over RE range
for RE in [1e5, 5e5]:
#for RE in [5e5]:
    print("RE",RE)
    # Settings
    ncrit = 10.2
    Mach = 0.04*RE/5e5 # c = 0.5m, assuming 20°C
    s = vf.setup(Re=RE,Ma=Mach,ncrit=ncrit,alpha=AOARange[0],tolerance = 1e-2)
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
        print("alpha",alpha)
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
        try:
                [x,flag,res,grad,_] = vf.iter(x,bl,p,s,res,grad)
        except:
            #faults += 1
            #init = True
            flag  = False
            #del vf
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

#EXPRES=np.genfromtxt("S805Polars.csv",delimiter=",",names=True)
fix,ax = plt.subplots(1,1)
ax.plot(results[1e5]["CD"],results[5e5]["CL"], color="tab:blue")
#ax.plot(EXPRES['EXPPOLAR5_X'],EXPRES['EXPPOLAR5_Y'],marker=".",linestyle = 'None', color="tab:blue")

#ax.plot(results[5e5]["CD"],results[7e5]["CL"],'tab:orange')
#ax.plot(EXPRES['EXPPOLAR7_X'],EXPRES['EXPPOLAR7_Y'],marker=".", linestyle = 'None',color='tab:orange')

#ax.plot(results[10e5]["CD"],results[10e5]["CL"],'tab:green')
#ax.plot(EXPRES['EXPPOLAR10_X'],EXPRES['EXPPOLAR10_Y'],marker=".", linestyle = 'None',color='tab:green')

#ax.plot(results[15e5]["CD"],results[15e5]["CL"],'tab:red')
#ax.plot(EXPRES['EXPPOLAR15_X'],EXPRES['EXPPOLAR15_Y'],marker=".", linestyle = 'None',color='tab:red')

#ax.plot(results[20e5]["CD"],results[20e5]["CL"],'tab:purple')
#ax.plot(EXPRES['EXPPOLAR20_X'],EXPRES['EXPPOLAR20_Y'],marker=".", linestyle = 'None',color='tab:purple')

ax.set_xlabel('CD')
ax.set_ylabel('CL')
ax.set_xlim([0,0.015])
ax.set_ylim([-.1,1.3])
ax.grid(True)
ax.set_title("Polars at different Reynolds Numbers");
plt.show()