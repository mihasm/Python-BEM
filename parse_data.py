__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Miha Smrekar"
__email__ = "miha.smrekar9@gmail.com"
__status__ = "Development"

from openpyxl import load_workbook
import numpy
from scipy import interpolate
from scipy.ndimage.interpolation import shift

try:
    wb = load_workbook('program_v3.xlsm')
    ws = wb['program']
except FileNotFoundError:
    print(
        "File 'program_v3.xlsm' was not found, please rename your excel file to 'program_v3.xlsm' and move it to the "
        "folder of this file.")


def parse_podatki():
    """
    Parses turbine geometry data from excel file program_v3.xlsm.
    Extracts cL and cD curve data.
    :return: Returns dictionary with resulting points as arrays.
    """
    ws2 = wb['podatki']

    COLUMN_CL_X = 1
    COLUMN_CL_Y = 2
    COLUMN_CD_X = 4
    COLUMN_CD_Y = 5

    start_row = 2

    cL_x = []
    cL_y = []
    cD_x = []
    cD_y = []

    i = start_row

    while True:
        value_x = ws2.cell(row=i, column=COLUMN_CL_X).value
        value_y = ws2.cell(row=i, column=COLUMN_CL_Y).value

        if value_x == None and value_y == None:
            break
        else:
            cL_x.append(float(value_x))
            cL_y.append(float(value_y))
            i += 1

    i = start_row

    while True:
        value_x = ws2.cell(row=i, column=COLUMN_CD_X).value
        value_y = ws2.cell(row=i, column=COLUMN_CD_Y).value

        if value_x == None and value_y == None:
            break
        else:
            cD_x.append(float(value_x))
            cD_y.append(float(value_y))
            i += 1

    return {
        'cL_x': numpy.array(cL_x),
        'cL_y': numpy.array(cL_y),
        'cD_x': numpy.array(cD_x),
        'cD_y': numpy.array(cD_y)
    }


def parse_sections():
    """
    Parses turbine geometry data from excel file program_v3.xlsm.
    :return: Returns section radiuses, chord lengths, chord angles and section heights.
    """
    sections_radius = []
    chord_lengths = []
    chord_angles = []

    for row in range(9, 31):
        sections_radius.append(ws.cell(row=row, column=1).value)  # mm
        chord_lengths.append(ws.cell(row=row, column=2).value)  # mm
        chord_angles.append(ws.cell(row=row, column=3).value)  # degrees

    sections_radius = numpy.array(sections_radius)
    sections_radius_shifted = shift(sections_radius, +1, cval=sections_radius[0] - 30)
    dr = sections_radius - sections_radius_shifted
    dr = dr * 1e-3  # m
    sections_radius = sections_radius * 1e-3

    chord_lengths = numpy.array(chord_lengths) * 1e-3  # m
    chord_angles = numpy.array(chord_angles)
    chord_angles = 90.0 - chord_angles
    return sections_radius, chord_lengths, chord_angles, dr


podatki = parse_podatki()
f_c_L = interpolate.interp1d(podatki['cL_x'], podatki['cL_y'], fill_value=(podatki['cL_y'][0], podatki['cL_y'][-1]),
                             bounds_error=False)
f_c_D = interpolate.interp1d(podatki['cD_x'], podatki['cD_y'], fill_value=(podatki['cD_y'][0], podatki['cD_y'][-1]),
                             bounds_error=False)
