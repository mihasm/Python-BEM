import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d as mp3d
from turbine_data import SET_INIT
import numpy as np
"""
bot = [(0, 0, 0),
       (1, 0, 0),
       (1, 1, 0),
       (0, 1, 0),
       ]

top = [(0, 0, 1),
       (1, 0, 1),
       (1, 1, 1),
       (0, 1, 1),
       ]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
face1 = mp3d.art3d.Poly3DCollection([bot], alpha=0.5, linewidth=1)
face2 = mp3d.art3d.Poly3DCollection([top], alpha=0.5, linewidth=1)

ax.add_collection3d(face1)
ax.add_collection3d(face2)
plt.show()
"""
import math

def rotate(origin, point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return qx, qy

def rotate_array(x,y,xy,angle):
	angle = math.radians(angle)
	x_out, y_out = [],[]
	for _x,_y in zip(x,y):
		#print(_x,_y)
		_xq,_yq = rotate(xy,(_x,_y),angle)
		x_out.append(_xq)
		y_out.append(_yq)
	return x_out,y_out

def get_centroid(foil_x,foil_y):
	centroid = (np.sum(foil_x) / len(foil_x), np.sum(foil_y) / len(foil_y))
	return centroid



theta = SET_INIT["theta_in"]
r = SET_INIT["r_in"]
c = SET_INIT["c_in"]
foil = SET_INIT["foils_in"]
airfoils = SET_INIT["airfoils"]
foil_x,foil_y = SET_INIT["airfoils"]["s826"]["x"],SET_INIT["airfoils"]["s826"]["y"]

def scale_and_normalize(foil_x,foil_y,scale):
	foil_x,foil_y = np.array(foil_x),np.array(foil_y)
	centroid = get_centroid(foil_x,foil_y)
	foil_x = foil_x-centroid[0]
	foil_x,foil_y = foil_x*scale,foil_y*scale
	return foil_x,foil_y

#foil_x,foil_y = scale_and_normalize(foil_x,foil_y,5)
rot_x,rot_y = rotate_array(foil_x,foil_y,(0,0),45)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
out_x,out_y,out_z = [],[],[]
for i in range(len(r)):
	z = r[i]
	_c = c[i]
	_foil = foil[i]
	_theta = theta[i]
	_foil_x,_foil_y = airfoils[_foil]["x"],airfoils[_foil]["y"]
	_foil_x,_foil_y = scale_and_normalize(_foil_x,_foil_y,_c)
	_foil_x,_foil_y = rotate_array(_foil_x,_foil_y,(0,0),_theta)
	for _x,_y in zip(_foil_x,_foil_y):
		out_x.append(_x)
		out_y.append(_y)
		out_z.append(z)


#plt.plot(foil_x,foil_y)
#plt.plot(rot_x,rot_y,"-r")
ax.scatter(out_x,out_y,out_z)
X,Y,Z = np.array(out_x),np.array(out_y),np.array(out_z)
max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max() / 2.0

mid_x = (X.max()+X.min()) * 0.5
mid_y = (Y.max()+Y.min()) * 0.5
mid_z = (Z.max()+Z.min()) * 0.5
ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)
#ax.set_aspect("equal")
#plt.axis("equal")
plt.show()