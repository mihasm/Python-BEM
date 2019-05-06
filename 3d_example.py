import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d as mp3d
from turbine_data import SET_INIT
import numpy as np
from matplotlib import cm
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
facecolors = [cm.jet(x) for x in np.random.rand(20)]


def create_face(p1,p2,p3,p4,*args,**kwargs):
	coords = [p1,p2,p3,p4]
	face = mp3d.art3d.Poly3DCollection([coords],facecolors=facecolors, alpha=1.0, linewidth=0.0, *args,**kwargs)
	return face

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
data_faces = []
for i in range(len(r)):
	z = r[i]
	data_faces.append([z])
	_c = c[i]
	_foil = foil[i]
	_theta = theta[i]
	_foil_x,_foil_y = airfoils[_foil]["x"],airfoils[_foil]["y"]
	_foil_x,_foil_y = scale_and_normalize(_foil_x,_foil_y,_c)
	_foil_x,_foil_y = rotate_array(_foil_x,_foil_y,(0,0),_theta)
	list_x,list_y = [],[]
	for _x,_y in zip(_foil_x,_foil_y):
		out_x.append(_x)
		out_y.append(_y)
		out_z.append(z)
		list_x.append(_x)
		list_y.append(_y)
	data_faces[-1].append(list_x)
	data_faces[-1].append(list_y)


#plt.plot(foil_x,foil_y)
#plt.plot(rot_x,rot_y,"-r")
#ax.scatter(out_x,out_y,out_z)
X,Y,Z = np.array(out_x),np.array(out_y),np.array(out_z)
max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max() / 2.0

mid_x = (X.max()+X.min()) * 0.5
mid_y = (Y.max()+Y.min()) * 0.5
mid_z = (Z.max()+Z.min()) * 0.5
ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)
ax.set_aspect("equal")
Axes3D.mouse_init(ax,rotate_btn=1,zoom_btn=2)


from matplotlib.colors import LightSource

#ls = LightSource(azdeg=0, altdeg=65)
# shade data, creating an rgb array.
#rgb = ls.shade(z, plt.cm.RdYlBu)


for i in range(len(data_faces)-1):
	z = data_faces[i][0]
	z_up = data_faces[i+1][0]
	for j in range(len(data_faces[i][1])-1):
		x = data_faces[i][1][j]
		y = data_faces[i][2][j]
		x_up = data_faces[i+1][1][j]
		y_up = data_faces[i+1][2][j]
		x_next = data_faces[i][1][j+1]
		y_next = data_faces[i][2][j+1]
		x_up_next = data_faces[i+1][1][j+1]
		y_up_next = data_faces[i+1][2][j+1]
		face = create_face((x,y,z),(x_next,y_next,z),(x_up_next,y_up_next,z_up),(x_up,y_up,z_up))
		ax.add_collection3d(face)


#face = create_face((0,0,0),(0,0,1),(1,0,0),(1,1,0))
#ax.add_collection3d(face)
#plt.axis("equal")
plt.show()

