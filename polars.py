from scrape_polars import get_polars
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from pprint import pprint
from scipy.interpolate import Rbf
from montgomerie import POLAR_CLASS, Montgomerie, draw_matplotlib
#from numpy import interp
from scipy.interpolate import interp1d
import scipy.interpolate
import matplotlib.cm as cm
from scipy.interpolate import griddata
import numpy

def interp_at(x, y, v, xp, yp, algorithm='cubic', extrapolate=False):
    """
    Interpolate data onto the specified points.
 
    Parameters:
 
    * x, y : 1D arrays
        Arrays with the x and y coordinates of the data points.
    * v : 1D array
        Array with the scalar value assigned to the data points.
    * xp, yp : 1D arrays
        Points where the data values will be interpolated
    * algorithm : string
        Interpolation algorithm. Either ``'cubic'``, ``'nearest'``,
        ``'linear'`` (see scipy.interpolate.griddata)
    * extrapolate : True or False
        If True, will extrapolate values outside of the convex hull of the data
        points.
 
    Returns:
 
    * v : 1D array
        1D array with the interpolated v values.
 
    """
    if algorithm not in ['cubic', 'linear', 'nearest']:
        raise ValueError("Invalid interpolation algorithm: " + str(algorithm))
    grid = scipy.interpolate.griddata((x, y), v, (xp, yp),
                                      method=algorithm).ravel()
    if extrapolate and algorithm != 'nearest' and numpy.any(numpy.isnan(grid)):
        grid = extrapolate_nans(xp, yp, grid)
    return grid

def extrapolate_nans(x, y, v):
    """
    Extrapolate the NaNs or masked values in a grid INPLACE using nearest
    value.
 
    .. warning:: Replaces the NaN or masked values of the original array!
 
    Parameters:
 
    * x, y : 1D arrays
        Arrays with the x and y coordinates of the data points.
    * v : 1D array
        Array with the scalar value assigned to the data points.
 
    Returns:
 
    * v : 1D array
        The array with NaNs or masked values extrapolated.
 
    """
    if numpy.ma.is_masked(v):
        nans = v.mask
    else:
        nans = numpy.isnan(v)
    notnans = numpy.logical_not(nans)
    v[nans] = scipy.interpolate.griddata((x[notnans], y[notnans]), v[notnans],
                                         (x[nans], y[nans]),
                                         method='nearest').ravel()
    return v
    

def get_cl_cd_interpolation_function(link,airfoil_x=[],airfoil_y=[]):
    #imp_polar = np.loadtxt(open("foils/NACA_0015_polar.csv", "rb"), delimiter=",", skiprows=1)
    print("scraping")
    data = get_polars(link)
    #pprint(data)
    x,y,z_cl,z_cd = [],[],[],[]
    ncrit_set = 9.0
    for Re,value in data.items():
        _alpha = []
        _cl = []
        _cd = []
        for ncrit,value2 in value.items():
            if ncrit == ncrit_set:
                for alpha,value3 in value2.items():
                    cl = value3["cl"]
                    cd = value3["cd"]
                    #print(Re,ncrit,alpha,cl)
                    _alpha.append(alpha)
                    _cl.append(cl)
                    _cd.append(cd)
        
        polar = POLAR_CLASS(x=airfoil_x,y=airfoil_y,alpha=_alpha,Cl=_cl,Cd=_cd)
        M = Montgomerie(polar,reynolds=Re)

        #m_Alpha,m_Cl,m_Cd = M.calculate_extrapolation([],[],[])
        m_Alpha,m_Cl,m_Cd = draw_matplotlib(polar,M)
        f_cl = interp1d(_alpha,_cl,bounds_error=True)
        f_cd = interp1d(_alpha,_cd,bounds_error=True)
        for i in range(len(m_Alpha)):
            x.append(Re)
            y.append(m_Alpha[i])
            try:
                cl = f_cl(m_Alpha[i])
            except ValueError:
                cl = m_Cl[i]
            try:
                cd = f_cd(m_Alpha[i])
            except ValueError:
                cd = m_Cd[i]
            z_cl.append(cl)
            z_cd.append(cd)

    #for i in range(len(x)):
    #    print(x[i],y[i],z_cl[i],z_cd[i])

    #print("test2")

    fig = plt.figure()
    ax = Axes3D(fig)
    x,y,z_cl,z_cd = np.array(x),np.array(y),np.array(z_cl),np.array(z_cd)
    #ax.scatter(x,y,z_cl)
    #plt.show()

    ax.set_xlabel('re')
    ax.set_ylabel('alpha')
    ax.set_zlabel('cl')
    xi, yi = np.linspace(x.min(), x.max(), 10), np.linspace(y.min(), y.max(), 361)
    num_elements = len(xi)*len(yi)

    xi, yi = np.meshgrid(xi, yi)


    print("interp1")
    zi_cl = interp_at(x,y,z_cl,xi.ravel(),yi.ravel(),algorithm="linear",extrapolate=False)
    #print(len(xi.flatten()),len(yi.flatten()),len(zi_cl.flatten()))
    
    rbf_cl = scipy.interpolate.interp2d(xi, yi, zi_cl, kind='linear') #interpolacija je lahko linear, cubic, itd.
    
    #xi2,yi2 = np.linspace(50000,500000,10),np.linspace(-180,180,361)
    #xi2,yi2 = np.meshgrid(xi2,yi2)
    #zi_cl_2 = rbf_cl(xi2.flatten(), yi2.flatten())
    



    print("interp2")
    zi_cd = interp_at(x,y,z_cd,xi.ravel(),yi.ravel(),algorithm="linear",extrapolate=False)
    rbf_cd = scipy.interpolate.Rbf(xi,yi, zi_cd, function='linear')
    zi_cd_2 = rbf_cd(xi,yi)
    


    #print(len(xi2.flatten()),len(yi2.flatten()),len(zi_cl_2.flatten()))
    #ax.scatter(xi,yi,zi_cl)
    #ax.scatter(xi2.flatten(),yi2.flatten(),zi_cl_2.flatten())
    #ax.scatter(xi,yi,zi_cd)
    #plt.show()
    return rbf_cl,rbf_cd

#f,f2 = get_cl_cd_interpolation_function("http://airfoiltools.com/airfoil/details?airfoil=s826-nr")
#print(f(100000,15))
#print(f(200000,5))
#print(f(50000,5.0))