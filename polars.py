from scrape_polars import get_polars
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from pprint import pprint
from scipy.interpolate import Rbf
from montgomerie import Montgomerie
from xfoil import generate_polars_data
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
    

def get_interpolation_function(x,y,z,num_x=10,num_y=361,min_x=0,max_x=1e6):
    x,y,z = np.array(x),np.array(y),np.array(z)
    xi, yi = np.linspace(x.min(), x.max(), num_x), np.linspace(y.min(), y.max(), num_y)
    xi, yi = np.meshgrid(xi, yi)
    zi = interp_at(x,y,z,xi.ravel(),yi.ravel(),algorithm="linear",extrapolate=False)   
    fun = scipy.interpolate.interp2d(xi, yi, zi, kind='linear')
    return fun

def scrape_data(link):
    out = []
    data = get_polars(link)
    for Re,value in data.items():
        for ncrit,value2 in value.items():
            for alpha,value3 in value2.items():
                cl = value3["cl"]
                cd = value3["cd"]
                out.append([Re,ncrit,alpha,cl,cd])
    out = np.array(out)
    return out

#print(scrape_data("http://airfoiltools.com/airfoil/details?airfoil=s826-nr"))

def get_cl_cd_from_link(link,airfoil_x=[],airfoil_y=[]):
    data = scrape_data(link)
    return get_cl_cd_interpolation_function(data,airfoil_x,airfoil_y)

def get_cl_cd_from_xfoil(foil,airfoil_x=[],airfoil_y=[]):
    data = generate_polars_data(foil)
    return get_cl_cd_interpolation_function(data,airfoil_x,airfoil_y)

def get_cl_cd_interpolation_function(data,airfoil_x=[],airfoil_y=[]):
    #imp_polar = np.loadtxt(open("foils/NACA_0015_polar.csv", "rb"), delimiter=",", skiprows=1)
    print("Getting inteprolation function")
    
    x,y,z_cl,z_cd = [],[],[],[]

    Re_list = np.unique(data[:,0])
    ncrit_list = np.unique(data[:,1])
    ncrit_selected = ncrit_list[0]

    for Re in Re_list:
        rows_with_Re = data[np.in1d(data[:,0],Re)]
        rows_with_Re = rows_with_Re[np.in1d(rows_with_Re[:,1],ncrit_selected)]

        _alpha = rows_with_Re[:,2].flatten()
        _cl = rows_with_Re[:,3].flatten()
        _cd = rows_with_Re[:,4].flatten()

        M = Montgomerie(x=airfoil_x,y=airfoil_y,alpha=_alpha,Cl=_cl,Cd=_cd,Re=Re)

        m_Alpha,m_Cl,m_Cd = M.calculate_extrapolation()
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

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_xlabel('re')
    ax.set_ylabel('alpha')
    ax.set_zlabel('cl')

    x,y,z_cl,z_cd = np.array(x),np.array(y),np.array(z_cl),np.array(z_cd)
    ax.scatter(x,y,z_cl)
    
    print("interp1")
    fun_cl = get_interpolation_function(x,y,z_cl)

    print("interp2")
    fun_cd = get_interpolation_function(x,y,z_cd)

    #plt.show()

    print("Done with interp!")

    return fun_cl,fun_cd



#f,f2 = get_cl_cd_from_link("http://airfoiltools.com/airfoil/details?airfoil=s826-nr")
#f,f2 = get_cl_cd_from_xfoil("s826.dat")
#print(f(200000,13))

