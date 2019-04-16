__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.3.0"
__maintainer__ = "Miha Smrekar"
__email__ = "miha.smrekar9@gmail.com"
__status__ = "Development"

import numpy
from scipy import interpolate
from numpy import array
import os


def transpose(a):
    """
    Transposes matrix-like array a.

    Input:
    [[1,2],
     [3,4],
     [5,6]]

    Output:
    [[1,3,5],
     [2,4,6]]

    :param a: input array
    :return: a^T
    """
    o = []
    for i in range(len(a[0])):
        o.append([])
        for _r in a:
            o[i].append(_r[i])
    return o


def dict_to_ar(inp_dict):
    prep = []
    i = 0
    for k, v in inp_dict.items():
        prep.append([k])
        for j in v:
            if isinstance(j, numpy.ndarray):
                j = numpy.array2string(j, max_line_width=10000000)
            prep[i].append(str(j))
        i += 1
    prep = transpose(prep)
    return prep


def dict_to_csv(inp_dict, delimiter=";"):
    """
    Creates string from input dict that is formatted like a csv file.

    E.g.:
        inp_dict={"a":[1,2,3],"b":4,"c":[5,6,7]}

        out:
        "a,b,c,
         1,4,5,
         2,4,6,
         3,4,7,"

    :param delimiter: delimiter string (default ;)
    :param inp_dict: input dictionary
    :return: csv-like string
    """
    prep = dict_to_ar(inp_dict)
    out = ""
    for _r in prep:
        for e in _r:
            out += str(e) + delimiter
        out += "\n"
    return out


def array_to_csv(in_ar, delimiter="\t"):
    out = ""
    for r in range(len(in_ar)):
        for c in in_ar[r]:
            if c == None:
                c = ""
            out += c + delimiter
        out = out[0:-1]
        out += "\n"
    return out


def sort_xy(array_x, array_y):
    """
    Sorts two arrays by values in the first array.
    Useful for sorting pairs of x,y points.

    Returns sorted arrays.

    Arrays have to be the same length.

    :param array_x: x values
    :param array_y: y values
    :return: x,y (sorted)
    """
    out = []
    if len(array_x) == len(array_y):
        for i in range(len(array_x)):
            out.append((array_x[i], array_y[i]))
        out = sorted(out, key=lambda k: k[0])
        out_x, out_y = [], []
        for n in range(len(out)):
            out_x.append(out[n][0])
            out_y.append(out[n][1])
        return out_x, out_y
    else:
        raise Exception(
            "Cannot create XY pairs with arrays with different num of elements"
        )


# noinspection PyUnboundLocalVariable
def transitions_to_nearest_profiles(r, foils):
    # print("transitioning foils")
    # print("transition foils before",foils)
    r_orig = r.copy()
    foils_orig = foils.copy()
    if len(r_orig) != len(foils_orig):
        raise Exception("Lengths of arrays r and foils must match!")
    for i in range(len(r)):
        # print("currently checking i",i)
        # if there is a transition between two foils,
        # select the closest foil to the given transition
        if "transition" in foils_orig[i].lower():
            _r = r[i]
            nearest_i_down = int(str(i))
            while True:
                nearest_i_down -= 1
                if nearest_i_down < 0:
                    nearest_i_down = None
                    break
                if "transition" in foils_orig[nearest_i_down].lower():
                    pass
                else:
                    break
            # print("##nearest_i_down",nearest_i_down)
            nearest_i_up = int(str(i))
            while True:
                nearest_i_up += 1
                if nearest_i_up >= len(r):
                    nearest_i_up = None
                    break
                if "transition" in foils_orig[nearest_i_up].lower():
                    pass
                else:
                    break
            # print("##nearest_i_up",nearest_i_up)

            if nearest_i_up == None and nearest_i_down == None:
                raise Exception("Transition %s doesn't end with any profile" % i)

            dr_down,dr_up = None, None

            if nearest_i_down != None:
                dr_down = abs(_r - r[nearest_i_down])

            if nearest_i_up != None:
                dr_up = abs(_r - r[nearest_i_up])

            if nearest_i_up != None and nearest_i_down != None:
                # noinspection PyUnboundLocalVariable
                if dr_down < dr_up:
                    best_i = nearest_i_down
                else:
                    best_i = nearest_i_up
            else:
                if nearest_i_up:
                    best_i = nearest_i_up
                else:
                    best_i = nearest_i_down

            foils[i] = foils_orig[best_i]
    # print("transition foils after",foils)
    return foils


def interpolate_geom(r, c, theta, foils, num=None, linspace_interp=False):
    """
    interpolates c,r,theta with num elements:
    """
    # print("interpolating")
    # print("foils before",foils)
    c_interpolator = interpolate.interp1d(r, c)
    theta_interpolator = interpolate.interp1d(r, theta)
    r_orig = r.copy()
    foils_orig = foils.copy()
    foils_orig = transitions_to_nearest_profiles(r_orig, foils_orig)
    if linspace_interp:
        r = numpy.linspace(start=r[0], stop=r[-1], num=num + 1)
        c = c_interpolator(r)
        theta = theta_interpolator(r)
        foils = []
        for _r in r:
            closest_index = find_nearest(r_orig, _r)
            foils.append(foils_orig[closest_index])
    else:
        foils = foils_orig

    # calculate dr
    r_shifted = [r[0]]
    for _r in r:
        r_shifted.append(_r)
    r_shifted = array(r_shifted[:-1])
    dr = r - r_shifted
    # print("foils after",foils)
    return r, c, theta, foils, dr


def find_nearest(_array, value):
    _array = numpy.asarray(_array)
    idx = (numpy.abs(_array - value)).argmin()
    return idx


def to_float(inp):
    if isinstance(inp, str):
        inp = inp.replace(",", ".")
    return float(inp)


class Printer:
    def __init__(self, arr):
        self.out = arr

    def print(self, *args, add_newline=True):
        out_str = ""
        i = 0
        for a in args:
            if i > 0:
                out_str += " "
            if isinstance(a, float):
                a = "%.3f" % round(a, 3)
            out_str += str(a)
            i += 1
        print(out_str)
        if add_newline:
            out_str += "\n"
        self.out.append(out_str)
        return out_str


import copy


def fltr(node, vals):
    print(node)
    if isinstance(node, dict):
        retVal = {}
        for key, value in node.items():
            if isinstance(value, numpy.ndarray):
                node[key] = value.tolist()
            if isinstance(key, vals) and isinstance(value, vals):
                retVal[key] = copy.deepcopy(node[key])
            elif isinstance(node[key], list) or isinstance(node[key], dict):
                child = fltr(node[key], vals)
                if child:
                    retVal[key] = child
        if retVal:
            return retVal
        else:
            return None
    elif isinstance(node, list):
        retVal = []
        for entry in node:
            child = fltr(entry, vals)
            if child:
                retVal.append(child)
        if retVal:
            return retVal
        else:
            return None


def generate_dat(name, x, y):
    print("generating dat")
    out = ""
    out += name + "\n"
    for i in range(len(x)):
        _x = float(x[i])
        _y = float(y[i])
        if _y >= 0:
            out += "%.6f   %.6f\n" % (_x, _y)
        else:
            out += "%.6f  %.6f\n" % (_x, _y)

    f = open(os.path.join("foils", name + ".dat"), "w")
    f.write(out)
    f.close()
    return out

def sort_data(data,columns=[0,2]):
    if len(columns) == 0:
        raise Exception("Sorting must be done for more than zero columns.")
    first = False
    for i in columns:
        if first == False:
            data = data[data[:,i].argsort()] #sort by reynolds
            first = True
        else:
            data = data[data[:,i].argsort(kind="mergesort")] #sort by alpha
    return data
