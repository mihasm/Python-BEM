__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.2.4"
__maintainer__ = "Miha Smrekar"
__email__ = "miha.smrekar9@gmail.com"
__status__ = "Development"

import numpy


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
