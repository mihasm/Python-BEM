import numpy as np
from polars import scrape_data,get_extrapolated_data
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
import scipy.interpolate


def interp(re_in,alpha_in,re,alpha,cl):
	"""
	data has to be sorted, first by 0th column, then by 2nd column
	"""

	re_list,alpha_list,cl_list = np.unique(re),np.unique(alpha),np.unique(cl)

	if re_in >= re_list.max():
		indexes = np.where(re == re_list.max())
		alpha_selected = alpha[indexes]
		cl_selected = cl[indexes]
		return np.interp(alpha_in,alpha_selected,cl_selected)

	if re_in <= re_list.min():
		indexes = np.where(re==re_list.min())
		alpha_selected = alpha[indexes]
		cl_selected = cl[indexes]
		return np.interp(alpha_in,alpha_selected,cl_selected)

	re_bottom_index = np.where(re_list<re_in)[0][-1]
	re_bottom = re_list[re_bottom_index]
	re_top = re_list[re_bottom_index+1]

	indexes_bottom = np.where(re == re_bottom)
	indexes_top = np.where(re == re_top)

	alpha_bottom = alpha[indexes_bottom]
	cl_bottom = cl[indexes_bottom]

	alpha_top = alpha[indexes_top]
	cl_top = cl[indexes_top]

	_cl_1 = np.interp(alpha_in,alpha_bottom,cl_bottom)
	_cl_2 = np.interp(alpha_in,alpha_top,cl_top)

	cL = (_cl_2-_cl_1)/(re_top-re_bottom)*(re_in-re_bottom)+_cl_1
	return cL


def interp_at(x, y, v, xp, yp, algorithm='linear', extrapolate=False):
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
    grid = scipy.interpolate.griddata((x, y), v, (xp, yp), method=algorithm).ravel()
    if extrapolate and algorithm != 'nearest' and numpy.any(numpy.isnan(grid)):
        if xp.size > 2:
            grid = extrapolate_nans(xp, yp, grid)
        else:
            return scipy.interpolate.griddata((x,y),v,(xp,yp),method="nearest")
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
    v[nans] = scipy.interpolate.griddata((x[notnans], y[notnans]), v[notnans], (x[nans], y[nans]),
                                         method='nearest').ravel()
    return v


def get_interpolation_function(x, y, z, num_x=10, num_y=360):
    x, y, z = np.array(x), np.array(y), np.array(z)
    xi, yi = np.linspace(x.min(), x.max(), num_x), np.linspace(y.min(), y.max(), num_y)
    xi, yi = np.meshgrid(xi, yi)
    zi = interp_at(x, y, z, xi.ravel(), yi.ravel(), algorithm="linear", extrapolate=True)
    fun = scipy.interpolate.interp2d(xi, yi, zi, kind='linear')
    return fun


"""
fig = pyplot.figure()
ax = Axes3D(fig)

data = scrape_data("http://airfoiltools.com/airfoil/details?airfoil=s826-nr")

#print(data)
#print(data)
re,alpha,z_cl,z_cd = get_extrapolated_data(data)

re_array = np.array(re)
alpha_array = np.array(alpha)
z_cl_array = np.array(z_cl)
out = []
for i in range(len(re_array)):
	out.append([re_array[i],alpha_array[i],z_cl_array[i]])
data = np.array(out)
data = data[data[:,0].argsort()]
data = data[data[:,1].argsort(kind="mergesort")]

#print(data)
x,y = np.linspace(10000,1e6,10),np.linspace(-180,180,100)
#xi,yi = np.meshgrid(x,y)
a,b,c = [],[],[]
for _x in x:
	for _y in y:
		_z = interp(_x,_y,data[:,0].flatten(),data[:,1].flatten(),data[:,2].flatten())
		a.append(_x)
		b.append(_y)
		c.append(_z)

ax.scatter(a,b,c)
pyplot.show()
"""