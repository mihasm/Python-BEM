from pprint import pprint

import matplotlib.cm as cm
import matplotlib.pyplot as plt
from montgomerie import Montgomerie
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.interpolate import Rbf,griddata,interp1d
from xfoil import generate_polars_data

def scrape_data(link):
    out = []
    data = get_polars(link)
    for Re, value in data.items():
        for ncrit, value2 in value.items():
            for alpha, value3 in value2.items():
                cl = value3["cl"]
                cd = value3["cd"]
                out.append([Re, ncrit, alpha, cl, cd])
    out = np.array(out)
    return out


def get_extrapolated_data(data, airfoil_x=[], airfoil_y=[]):
    # imp_polar = np.loadtxt(open("foils/NACA_0015_polar.csv", "rb"), delimiter=",", skiprows=1)
    print("Getting inteprolation function")

    out = []

    x, y, z_cl, z_cd = [], [], [], []

    Re_list = np.unique(data[:, 0])
    ncrit_list = np.unique(data[:, 1])
    ncrit_selected = ncrit_list[0]

    for Re in Re_list:
        rows_with_Re = data[np.in1d(data[:, 0], Re)]
        rows_with_Re = rows_with_Re[np.in1d(
            rows_with_Re[:, 1], ncrit_selected)]

        _alpha = rows_with_Re[:, 2].flatten()
        _cl = rows_with_Re[:, 3].flatten()
        _cd = rows_with_Re[:, 4].flatten()

        M = Montgomerie(x=airfoil_x, y=airfoil_y,
                        alpha=_alpha, Cl=_cl, Cd=_cd, Re=Re)

        m_Alpha, m_Cl, m_Cd = M.calculate_extrapolation()
        f_cl = interp1d(_alpha, _cl, bounds_error=True)
        f_cd = interp1d(_alpha, _cd, bounds_error=True)
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
            out.append([Re, ncrit_selected, m_Alpha[i], cl, cd])
    return np.array(out)

import requests
from bs4 import BeautifulSoup as bs
from urllib.parse import parse_qsl, urlparse
from urllib import parse


def get_polars(link):
    results = {}
    print("Getting data from",link)
    r = requests.get(link)
    txt = r.text
    soup = bs(txt, features='html.parser')
    table_of_polars = soup.find_all("table", {"class": "polar"})[0]
    rows = table_of_polars.find_all("tr")[1:-1]
    links = []
    for r in rows:
        data = r.find_all("td")
        checkbox, name, reynolds, ncrit, maxclcd, description, source, details = [
            d for d in data]
        details_link = "http://airfoiltools.com" + \
            details.find_all("a")[0].get('href')
        links.append(details_link)

    print("List of links:",links)

    csv_links = []

    for l in links:
        print("getting csv links from",l)
        r = requests.get(l)
        txt = r.text
        soup = bs(txt, features='html.parser')
        details_table = soup.find_all("table", {"class": "details"})[0]
        td1 = details_table.find_all("td", {"class": "cell1"})
        _links = details_table.find_all("a")
        for _l in _links:
            if "csv?polar" in _l.get("href"):
                csv_links.append("http://airfoiltools.com" + _l.get("href"))
                break
    
    print("CSV links:",csv_links)

    for l in csv_links:
        print("getting data from csv link",l)
        lines = None
        r = requests.get(l)
        text = r.text
        lines = text.splitlines()
        start = lines.index("Alpha,Cl,Cd,Cdp,Cm,Top_Xtr,Bot_Xtr")
        reynolds_number = [float(i.split(",")[1]) for i in lines if "Reynolds number," in i][0]
        ncrit = [float(i.split(",")[1]) for i in lines if "Ncrit," in i][0]
        mach = [float(i.split(",")[1]) for i in lines if "Mach," in i][0]

        for line_num in range(start + 1, len(lines)):
            alpha, cl, cd = None, None, None
            alpha, cl, cd = lines[line_num].split(",")[:3]
            alpha, cl, cd = float(alpha), float(cl), float(cd)
            if not reynolds_number in results.keys():
                results[reynolds_number] = {}
            if not ncrit in results[reynolds_number].keys():
                results[reynolds_number][ncrit] = {}
            results[reynolds_number][ncrit][alpha] = {"cl": cl, "cd": cd}

    return results


def get_x_y_from_link(link):
    print("Getting x-y airfoil data from",link)

    params = parse_qsl(urlparse(link.strip()).query, keep_blank_values=True)
    selig_link = "http://airfoiltools.com/airfoil/seligdatfile?airfoil="+params[0][1]
    print("Presumed Selig dat link:")
    print(selig_link)

    r = requests.get(selig_link)
    text = r.text
    print("Got text:")
    print(text)
    print("Parsing...")
    lines = text.splitlines()

    x,y=[],[]

    for l in lines[1:]:
        stripped_line = l.strip().replace("  "," ")
        x_l,y_l = stripped_line.split(" ")
        x.append(float(x_l))
        y.append(float(y_l))

    return x,y