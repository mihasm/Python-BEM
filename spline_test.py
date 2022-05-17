from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from math import sqrt
from numpy import ones,vstack
from numpy.linalg import lstsq


"""
Start of a method to determine the interpolated aerodynamic profile based on the smaller and the bigger profile.

1. Define a spline, and redraw the profile using n points (100 or more).
2. Find the nearest point.
3. Determine the equation of the line passing through the nearest point on the spline and the point on the other profile.
4. Create the interpolated point on the line, on a distance, which is determined by the interpolation factor.
5. Repeat for every point on the other profile.
"""

array_a = np.array( [ [1, 2], [4,8], [6, 15], [10, 6], [10, 3] ] )
x = array_a[:,0]
y = array_a[:,1]
tck,u     = interpolate.splprep( [x,y] ,s = 0 )
xnew,ynew = interpolate.splev( np.linspace( 0, 1, 200 ), tck, der=0)

def closest_point(x_point,y_point, x_list,y_list):
    dist = [0]*len(x_list)
    for i in range(len(x_list)):
        dist_v = np.sqrt((x_list[i]-x_point)**2+(y_list[i]-y_point)**2)
        dist[i] = dist_v
    min_index = np.argmin(dist)
    return min_index,x_list[min_index],y_list[min_index]
 

x_point = 2
y_point = 10

index,x_found,y_found = closest_point(x_point,y_point,xnew,ynew)

def line_coefficients_two_points(x1,y1,x2,y2):
    A = vstack([(x1,x2),ones(len((x1,x2)))]).T
    k, c = lstsq(A, (y1,y2),rcond=None)[0]
    return k,c

def distance_between_two_points()

k,c = line_coefficients_two_points(0,0,1,1)


plt.plot( x,y,'o' , xnew ,ynew,".")
plt.plot( x_point,y_point,"r*")
plt.plot( x_found,y_found,"r.")
plt.legend( [ 'data' , 'spline'] )
plt.axis([0,20,0,20])
plt.show()
