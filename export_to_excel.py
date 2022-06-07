from openpyxl import Workbook
from openpyxl.drawing.image import Image


from test_results import res
from utils import transpose
import json
from matplotlib import pyplot as plt

SET_INIT = json.loads(open("karlsen.bem").read())

def generate_excel(settings, results):
    """

    :param settings:
    :param results:
    """
    wb = Workbook()
    ws = wb.active
    ws.title = "Turbine info"

    ### ADD BASIC DATA
    ws.append(["Turbine name", settings['turbine_name']])
    ws.append(["Tip radius", settings['R'], 'm'])
    ws.append(["Hub radius", settings['Rhub'], 'm'])
    ws.append(["Number of blades", settings['B']])
    ws.append(["Fluid density", settings['rho'], 'kg/m3'])
    ws.append(["Kinematic viscosity", settings['kin_viscosity'], 'm2/s'])

    ### ADD GEOMETRY DATA
    data_start = 10
    num_sections = len(settings['r'])

    row = data_start
    column = 1
    data_turbine = transpose([settings['r'], settings['c'], settings['theta'], settings['foils']])
    data_turbine = [
                       ["r [m]", 'c [m]', 'theta [m]', 'profil [m]']
                   ] + data_turbine
    for r in data_turbine:
        column = 1
        for v in r:
            ws.cell(row=row, column=column, value=v)
            column += 1
        row += 1
    
    ### ADD R,C,THETA CHART

    fig, ax1 = plt.subplots()
    color = 'tab:blue'

    ax1.set_xlabel('r [m]')
    ax1.set_ylabel('c [m]', color=color)
    ax1.plot(settings['r'],settings['c'], color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:red'
    ax2.set_ylabel('θ [°]', color=color)  # we already handled the x-label with ax1
    ax2.plot(settings['r'],settings['theta'], color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()  # otherwise the right y-label is slightly clipped

    plt.savefig("graph.png")

    image = Image('graph.png')

    ws.add_image(image,'E2')


    wb.save('test.xlsx')


generate_excel(SET_INIT, res)
