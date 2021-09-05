import os

import numpy as np
from PyQt5.QtWidgets import QWidget, QGridLayout, QFormLayout, QLabel, QLineEdit, QComboBox, QCheckBox, QPushButton, \
                            QScrollArea, QVBoxLayout
from numpy.core._multiarray_umath import array
from scipy import interpolate

from main import application_path
from UI.helpers import MatplotlibWindow, PrintoutWindow, AdkinsThread
from turbine_data import SET_INIT
from utils import to_float, interpolate_geom, generate_chord_lengths_betz, generate_twists_betz, \
    generate_chord_lengths_schmitz, generate_twists_schmitz, create_folder, create_macro_text, \
    generate_propeller_larabee, generate_propeller_adkins
from UI.Table import Table
from visualize import create_3d_blade


class WindTurbineProperties(QWidget):
    """
    Class used for storing the main wind turbine information, such as its name, number of blades, radiuses, and
    its geometry information (radius, chord length, twist angle and airfoil name for each wind blade section).
    """

    def __init__(self, parent=None):
        super(WindTurbineProperties, self).__init__(parent)

        self.main = self.parent()

        self.grid = QGridLayout()
        self.setLayout(self.grid)

        self.left = QWidget()
        self.fbox = QFormLayout()
        self.left.setLayout(self.fbox)

        self.scroll_area = QScrollArea()
        self.scroll_widget = QWidget()
        self.scroll_widget_layout = QVBoxLayout()

        self.scroll_widget.setLayout(self.scroll_widget_layout)
        self.scroll_area.setWidget(self.left)
        self.scroll_area.setWidgetResizable(True)

        self.grid.addWidget(self.scroll_area, 1, 1)

        self.table_properties = Table()
        self.table_properties.createEmpty(4, 30)
        self.table_properties.set_labels(["r [m]", "c [m]", "theta [deg]", "airfoil"])

        #self.grid.addWidget(self.left, 1, 1)
        self.grid.addWidget(self.table_properties, 1, 2)

        _name = QLabel("Turbine Name")
        self.name = QLineEdit()
        self.fbox.addRow(_name, self.name)
        self.name.textEdited.connect(self.main.set_title)

        _Rhub = QLabel("Hub radius [m]")
        self.Rhub = QLineEdit()
        self.Rhub.setText("0.1")
        self.fbox.addRow(_Rhub, self.Rhub)
        self.Rhub.setToolTip("Radij pesta. Vpliv ima le pri Hub-Loss popravku.")

        _R = QLabel("Tip radius [m]")
        self.R = QLineEdit()
        self.R.setText("0.776")
        self.fbox.addRow(_R, self.R)
        self.R.setToolTip("Radij turbine. Vpliv ima na izračun moči ter pri Tip-Loss popravku.")

        _B = QLabel("Number of blades")
        self.B = QLineEdit()
        self.B.setText("5")
        self.fbox.addRow(_B, self.B)
        self.B.setToolTip("Število lopatic.")

        _blade_design = QLabel("Blade design")
        self.blade_design = QComboBox()
        self.blade_design.addItems(["Filled", "Hollow", "Spar"])
        self.fbox.addRow(_blade_design, self.blade_design)
        self.B.setToolTip("Pomembno za statični izračun.")

        _blade_thickness = QLabel("Blade thickness [m]")
        self.blade_thickness = QLineEdit()
        self.blade_thickness.setText("0.001")
        self.fbox.addRow(_blade_thickness, self.blade_thickness)
        self.blade_thickness.setToolTip("Debelina lopatice / Spar Cap-a.")

        _mass_density = QLabel("Mass density [kg/m^3]")
        self.mass_density = QLineEdit()
        self.mass_density.setText("1060")
        self.fbox.addRow(_mass_density, self.mass_density)
        self.mass_density.setToolTip("Debelina lopatice / Spar Cap-a.")

        self.fbox.addRow(QLabel("————— Scale and interpolation —————"))

        _geometry_scale = QLabel("Scale factor")
        self.geometry_scale = QLineEdit()
        self.geometry_scale.setText("1.0")
        self.fbox.addRow(_geometry_scale, self.geometry_scale)
        self.geometry_scale.setToolTip("Scale factor")

        _linspace_interp = QLabel("Interpolate geometry")
        self.linspace_interp = QCheckBox()
        self.fbox.addRow(_linspace_interp, self.linspace_interp)
        self.linspace_interp.setToolTip("Interpolate_geom")

        _num_interp = QLabel("Number of interpolation points")
        self.num_interp = QLineEdit()
        self.num_interp.setText("25")
        self.fbox.addRow(_num_interp, self.num_interp)
        self.num_interp.setToolTip("Number of interpolation points")

        self.fbox.addRow(QLabel("————— Generate geometry —————"))

        _num_gen_sections = QLabel("Number of gen. sections")
        self.num_gen_sections = QLineEdit()
        self.num_gen_sections.setText("10")
        self.fbox.addRow(_num_gen_sections, self.num_gen_sections)
        self.num_gen_sections.setToolTip("Število odsekov za generiranje.")

        _design_tsr = QLabel("Design TSR")
        self.design_tsr = QLineEdit()
        self.design_tsr.setText("7")
        self.fbox.addRow(_design_tsr, self.design_tsr)
        self.design_tsr.setToolTip("Željeni TSR za generiranje.")

        _design_aoa = QLabel("Design AoA [°]")
        self.design_aoa = QLineEdit()
        self.design_aoa.setText("7")
        self.fbox.addRow(_design_aoa, self.design_aoa)
        self.design_aoa.setToolTip("Željeni AoA za generiranje.")

        _design_cl = QLabel("Cl @ Design AoA")
        self.design_cl = QLineEdit()
        self.design_cl.setText("1.4")
        self.fbox.addRow(_design_cl, self.design_cl)
        self.design_cl.setToolTip("Koeficient vzgona pri željenem AoA.")

        _design_airfoil = QLabel("Airfoil")
        self.design_airfoil = QLineEdit()
        self.design_airfoil.setText("s826")
        self.fbox.addRow(_design_airfoil, self.design_airfoil)
        self.design_airfoil.setToolTip("Tekst za v stolpec airfoil.")

        _design_RPM = QLabel("RPM (prop)")
        self.design_RPM = QLineEdit()
        self.design_RPM.setText("5000")
        self.fbox.addRow(_design_RPM, self.design_RPM)
        self.design_RPM.textChanged.connect(self.update_j)
        self.design_RPM.setToolTip("RPM propelerja")

        _design_velocity = QLabel("Velocity (prop) [m/s]")
        self.design_velocity = QLineEdit()
        self.design_velocity.setText("15")
        self.fbox.addRow(_design_velocity, self.design_velocity)
        self.design_velocity.textChanged.connect(self.update_j)
        self.design_velocity.setToolTip("Hitrost letala")

        _design_thrust = QLabel("Thrust (prop) [N]")
        self.design_thrust = QLineEdit()
        self.design_thrust.setText("50")
        self.design_thrust.textChanged.connect(self.update_j)
        self.fbox.addRow(_design_thrust, self.design_thrust)
        self.design_thrust.setToolTip("Potisk")

        _design_drag_lift_ratio = QLabel("Cd/Cl @ Design AoA (prop)")
        self.design_drag_lift_ratio = QLineEdit()
        self.design_drag_lift_ratio.setText("0.02")
        self.fbox.addRow(_design_drag_lift_ratio, self.design_drag_lift_ratio)
        self.design_drag_lift_ratio.setToolTip("Razmerje Drag/Lift (samo pri Larrabee)")

        _design_rho = QLabel("Air density [kg/m3]")
        self.design_rho = QLineEdit()
        self.design_rho.setText("1.225")
        self.fbox.addRow(_design_rho, self.design_rho)
        self.design_rho.setToolTip("Gostota zraka")

        _design_iters = QLabel("Num. Iterations")
        self.design_iters = QLineEdit()
        self.design_iters.setText("15")
        self.fbox.addRow(_design_iters, self.design_iters)
        self.design_iters.setToolTip("St. iteracij (Adkins)")

        _design_relaxation = QLabel("Relax. factor (Adkins)")
        self.design_relaxation = QLineEdit()
        self.design_relaxation.setText("0.3")
        self.fbox.addRow(_design_relaxation, self.design_relaxation)
        self.design_relaxation.setToolTip("Relaksacijski faktor (Adkins)")

        _convergence_criterion_adkins = QLabel("Convergence crit.")
        self.convergence_criterion_adkins = QLineEdit()
        self.convergence_criterion_adkins.setText("0.05")
        self.fbox.addRow(_convergence_criterion_adkins, self.convergence_criterion_adkins)
        self.convergence_criterion_adkins.setToolTip("Konvergenčni kriterij")

        _design_method = QLabel("Design method.")
        self.design_method = QComboBox()
        self.design_method.addItems(["Betz", "Schmitz", "Larrabee (prop)", "Adkins (prop)"])
        self.fbox.addRow(_design_method, self.design_method)
        self.design_method.setToolTip("Metoda dizajniranja.")

        _button_generate_geometry = QLabel("Generate geometry.")
        self.button_generate_geometry = QPushButton("Generate")
        self.fbox.addRow(_button_generate_geometry, self.button_generate_geometry)
        self.button_generate_geometry.clicked.connect(self.generate_geometry)
        self.button_generate_geometry.setToolTip("Generiraj geometrijo (povozi predhodno!).")

        self.fbox.addRow(QLabel("————— Export to Solidworks —————"))

        self.export_button = QPushButton("Export curve data")
        self.export_button.clicked.connect(self.export)
        self.fbox.addRow("Export:", self.export_button)
        self.export_button.setToolTip(
            "Krivulje na vseh radijih lopatice se shranijo v posamezne datoteke. Solidworks makro se nato zgenerira v Python konzoli.")

        self.flip_turning_direction = QCheckBox()
        self.fbox.addRow("Flip turning direction", self.flip_turning_direction)

        self.propeller_geom = QCheckBox()
        self.fbox.addRow("Propeller", self.propeller_geom)

        self.fbox.addRow(QLabel("—————————————————————————"))

        _button_calculate_pitch = QLabel("Calculate pitch")
        self.button_calculate_pitch = QPushButton("Calculate pitch")
        self.fbox.addRow(_button_calculate_pitch, self.button_calculate_pitch)
        self.button_calculate_pitch.clicked.connect(self.calculate_pitch)
        self.button_calculate_pitch.setToolTip("Izračun zasuka propelerja v inčih")

        

        _button_create_geometry_graph = QLabel("Create R,C,θ graph.")
        self.button_create_geometry_graph = QPushButton("Create R,C,θ graph.")
        self.fbox.addRow(_button_create_geometry_graph, self.button_create_geometry_graph)
        self.button_create_geometry_graph.clicked.connect(self.create_geometry_graph)
        self.button_create_geometry_graph.setToolTip("Izris grafa R,C,θ.")

        _button_create_pitch_graph = QLabel("Create pitch graph.")
        self.button_create_pitch_graph = QPushButton("Create pitch graph.")
        self.fbox.addRow(_button_create_pitch_graph, self.button_create_pitch_graph)
        self.button_create_pitch_graph.clicked.connect(self.create_pitch_graph)
        self.button_create_pitch_graph.setToolTip("Izris grafa koraka propelerja.")

        self.pitch_inches = QLabel("0")
        self.fbox.addRow("Pitch @ 0.75R [in.]", self.pitch_inches)

        self.J_string = QLabel("0")
        self.fbox.addRow("J (prop sign):", self.J_string)

        self.c_T_string = QLabel("0")
        self.fbox.addRow("c_T@J (prop design):", self.c_T_string)

        self.window = None
        self.thread = None

    def get_settings(self):
        """
        Used to get the basic wind turbine settings in a dictionary format, namely the:
        -hub radius (Rhub)
        -tip radius (R)
        -number of blades (B)
        -turbine name (turbine_name)
        -geometry (r,c,theta,foils) as four separate list objects
        :return: dict: Settings dictionary (Basic wind turbine information)
        """
        out = {"Rhub": to_float(self.Rhub.text()),
               "R": to_float(self.R.text()),
               "B": int(self.B.text()),
               "turbine_name": self.name.text(),
               "geometry_scale": to_float(self.geometry_scale.text()),
               "linspace_interp": self.linspace_interp.isChecked(),
               "num_interp": int(self.num_interp.text()),
               "blade_thickness": to_float(self.blade_thickness.text()),
               "blade_design": self.blade_design.currentIndex(),
               "mass_density": to_float(self.mass_density.text()),
               "design_RPM":to_float(self.design_RPM.text()),
               "design_velocity":to_float(self.design_velocity.text()),
               "design_thrust":to_float(self.design_thrust.text()),
               "design_drag_lift_ratio":to_float(self.design_drag_lift_ratio.text()),
               "design_rho":to_float(self.design_rho.text()),
               "num_gen_sections":to_float(self.num_gen_sections.text()),
               "design_tsr":to_float(self.design_tsr.text()),
               "design_aoa":to_float(self.design_aoa.text()),
               "design_cl":to_float(self.design_cl.text()),
               "design_iters":to_float(self.design_iters.text()),
               "design_relaxation":to_float(self.design_relaxation.text()),
               "design_airfoil":self.design_airfoil.text(),
               "design_method":to_float(self.design_method.currentIndex()),
               "convergence_criterion_adkins":to_float(self.convergence_criterion_adkins.text())}
        
        geom_array = self.table_properties.get_values()
        r, c, theta, foils = [], [], [], []
        for row in geom_array:
            if row[0] != "" and row[1] != "" and row[2] != "":
                r.append(to_float(row[0]))
                c.append(to_float(row[1]))
                theta.append(to_float(row[2]))
                foils.append(row[3])
        out["r"] = array(r)
        out["c"] = array(c)
        out["theta"] = array(theta)
        out["foils"] = foils
        _r = out["r"]
        _c = out["c"]
        _theta = out["theta"]
        _foils = out["foils"]
        out["R"] = out["R"]
        out["Rhub"] = out["Rhub"]
        r, c, theta, foils, dr = interpolate_geom(_r, _c, _theta, _foils, out["R"], out["Rhub"], out["num_interp"], out["linspace_interp"],
                                                  out["geometry_scale"])
        out["r"], out["c"], out["theta"], out["foils"], out["dr"] = r, c, theta, foils, dr
        out["r_in"], out["c_in"], out["theta_in"], out["foils_in"] = _r, _c, _theta, _foils

        out["propeller_pitch"] = 2*np.pi*np.array(out["r"])*np.tan(np.radians(out["theta"]))*39.3701

        return out

    def update_j(self):
        try:
            J = float(self.design_velocity.text()) / (float(self.design_RPM.text()) / 60 * 2 * float(self.R.text()))
            self.J_string.setText("%.2f" % (J))
            c_T=float(self.design_thrust.text())/(float(self.design_rho.text())*(float(self.design_RPM.text())/60)**2*(float(self.R.text())*2)**4)
            self.c_T_string.setText("%.2f" % (c_T))
        except:
            print("couldnt update J/cT")

    def calculate_pitch(self):
        out = self.get_settings()
        r,theta = np.array(out["r_in"]), np.array(out["theta_in"])
        R = out["R"]
        r_R = r/R
        f = interpolate.interp1d(r_R,theta)
        beta = np.deg2rad(f(0.75)) # theta at 0.75 R
        D = 2*R*39.3701 # inches
        pitch_inches = np.pi*0.75*D*np.tan(beta)
        self.pitch_inches.setText(str(pitch_inches))
        return


    def create_geometry_graph(self):
        out = self.get_settings()
        self.gw = MatplotlibWindow()
        self.gw.setWindowTitle("r,c,θ graph")

        self.gw.ax = self.gw.figure.add_subplot(111)
        self.gw.ax.set_title("c(r) and θ(r)")

        self.gw.ax2 = self.gw.ax.twinx()
        self.gw.ax.plot(out["r"], out["c"], color="b")
        self.gw.ax2.plot(out["r"], out["theta"], color="r")

        self.gw.ax.set_xlabel("Radius r [m]")
        self.gw.ax.set_ylabel("Chord c [m]", color="tab:blue")
        self.gw.ax2.set_ylabel("Twist θ [°]", color="tab:red")
        self.gw.ax.tick_params(axis='y', labelcolor="tab:blue")
        self.gw.ax2.tick_params(axis='y', labelcolor="tab:red")

    def create_pitch_graph(self):
        out = self.get_settings()
        self.gw_pitch = MatplotlibWindow()
        self.gw_pitch.setWindowTitle("Pitch graph")

        self.gw_pitch.ax = self.gw_pitch.figure.add_subplot(111)
        self.gw_pitch.ax.set_title("p(r)")
        #print(out["propeller_pitch"])
        self.gw_pitch.ax.plot(np.array(out["r"])/out["R"], out["propeller_pitch"])
        self.gw_pitch.ax.set_xlabel("r/R")
        self.gw_pitch.ax.set_ylabel("Pitch p [in]")

    def adkins_completion(self):
        array_out = []
        chords,thetas = self.adkins_return
        R = float(self.R.text())
        Rhub = float(self.Rhub.text())
        num_gen_sections = int(float(self.num_gen_sections.text()))
        radiuses = np.linspace(Rhub, R, num_gen_sections)
        airfoil = self.design_airfoil.text()
        for r in range(num_gen_sections):
            array_out.append([round(radiuses[r], 4), round(chords[r], 4), round(thetas[r], 4), airfoil])

        self.table_properties.createTable(array_out)

        self.calculate_pitch()

    def generate_geometry(self):
        array_out = []
        R = float(self.R.text())
        Rhub = float(self.Rhub.text())
        num_gen_sections = int(float(self.num_gen_sections.text()))
        radiuses = np.linspace(Rhub, R, num_gen_sections)
        Cl_max = float(self.design_cl.text())
        B = float(self.B.text())
        TSR = float(self.design_tsr.text())
        method = self.design_method.currentIndex()
        airfoil = self.design_airfoil.text()
        design_aoa = float(self.design_aoa.text())
        RPM = float(self.design_RPM.text())
        v = float(self.design_velocity.text())
        T = float(self.design_thrust.text())
        drag_lift_ratio = float(self.design_drag_lift_ratio.text())
        rho = float(self.design_rho.text())

        if method == 0:
            chords = generate_chord_lengths_betz(radiuses=radiuses, R=R, Cl_max=Cl_max, B=B, TSR=TSR)
            thetas = generate_twists_betz(radiuses=radiuses, R=R, TSR=TSR, alpha_d=design_aoa)

        elif method == 1:
            chords = generate_chord_lengths_schmitz(radiuses=radiuses, R=R, Cl_max=Cl_max, B=B, TSR=TSR)
            thetas = generate_twists_schmitz(radiuses=radiuses, R=R, TSR=TSR, alpha_d=design_aoa)

        elif method == 2:
            chords, thetas = generate_propeller_larabee(radiuses=radiuses, R=R, B=B, RPM=RPM, drag_lift_ratio=drag_lift_ratio, v=v, T=T, rho=rho, cl=Cl_max)

        elif method == 3:
            input_data = self.main.get_all_settings()
            if self.window != None:
                self.window.close()
            self.window = PrintoutWindow(self)
            if self.thread != None:
                if self.thread.isRunning():
                    self.thread.terminate()
            self.thread = AdkinsThread(self)
            self.thread.set_params(input_data)
            self.thread.completeSignal.connect(self.adkins_completion)
            self.thread.start()
            print("Done")
            return

        for r in range(num_gen_sections):
            array_out.append([round(radiuses[r], 4), round(chords[r], 4), round(thetas[r], 4), airfoil])

        self.table_properties.createTable(array_out)

        self.calculate_pitch()

    def set_settings(self, dict_settings):
        """
        Reads the settings from the input dictionary and sets the values in the UI.
        :param dict_settings: dict: Settings dictionary.
        """
        if "Rhub" in dict_settings:
            t = str(dict_settings["Rhub"])
            self.Rhub.setText(t)
        if "R" in dict_settings:
            t = str(dict_settings["R"])
            self.R.setText(t)
        if "geometry_scale" in dict_settings:
            t = str(dict_settings["geometry_scale"])
            self.geometry_scale.setText(t)
        if "linspace_interp" in dict_settings:
            t = dict_settings["linspace_interp"]
            self.linspace_interp.setChecked(t)
        if "num_interp" in dict_settings:
            t = str(dict_settings["num_interp"])
            self.num_interp.setText(t)
        if "B" in dict_settings:
            t = str(dict_settings["B"])
            self.B.setText(t)
        if "r_in" in dict_settings and "c_in" in dict_settings and "theta_in" in dict_settings and "foils_in" in dict_settings:
            _array = []
            for r in range(len(dict_settings["r_in"])):
                _r = dict_settings["r_in"][r]
                _c = dict_settings["c_in"][r]
                _theta = dict_settings["theta_in"][r]
                _f = dict_settings["foils_in"][r]
                _array.append([_r, _c, _theta, _f])
            self.table_properties.createTable(_array)
        if "turbine_name" in dict_settings:
            t = str(dict_settings["turbine_name"])
            self.name.setText(t)
        else:
            self.name.setText("")
        if "blade_thickness" in dict_settings:
            t = str(dict_settings["blade_thickness"])
            self.blade_thickness.setText(t)
        if "blade_design" in dict_settings:
            self.blade_design.setCurrentIndex(dict_settings["blade_design"])
        if "mass_density" in dict_settings:
            t = str(dict_settings["mass_density"])
            self.mass_density.setText(t)
        
        if "num_gen_sections" in dict_settings:
            t = str(dict_settings["num_gen_sections"])
            self.num_gen_sections.setText(t)
        if "design_tsr" in dict_settings:
            t = str(dict_settings["design_tsr"])
            self.design_tsr.setText(t)
        if "design_aoa" in dict_settings:
            t = str(dict_settings["design_aoa"])
            self.design_aoa.setText(t)
        if "design_cl" in dict_settings:
            t = str(dict_settings["design_cl"])
            self.design_cl.setText(t)
        if "design_airfoil" in dict_settings:
            t = str(dict_settings["design_airfoil"])
            self.design_airfoil.setText(t)
        if "design_RPM" in dict_settings:
            t = str(dict_settings["design_RPM"])
            self.design_RPM.setText(t)
        if "design_velocity" in dict_settings:
            t = str(dict_settings["design_velocity"])
            self.design_velocity.setText(t)
        if "design_thrust" in dict_settings:
            t = str(dict_settings["design_thrust"])
            self.design_thrust.setText(t)
        if "design_drag_lift_ratio" in dict_settings:
            t = str(dict_settings["design_drag_lift_ratio"])
            self.design_drag_lift_ratio.setText(t)
        if "design_rho" in dict_settings:
            t = str(dict_settings["design_rho"])
            self.design_rho.setText(t)
        if "design_iters" in dict_settings:
            t = str(dict_settings["design_iters"])
            self.design_iters.setText(t)
        if "design_relaxation" in dict_settings:
            t = str(dict_settings["design_relaxation"])
            self.design_relaxation.setText(t)
        if "convergence_criterion_adkins" in dict_settings:
            t = str(dict_settings["convergence_criterion_adkins"])
            self.convergence_criterion_adkins.setText(t)
        if "design_method" in dict_settings:
            t = dict_settings["design_method"]
            self.design_method.setCurrentIndex(t)

        self.calculate_pitch()
        self.update_j()

    def export(self):
        """
        Exports the wind turbine geometry data as spatial data (XYZ points), and creates a VB macro, which
        you can use to import the geometry into the Solidworks 3D modeller.
        :return:
        """

        if self.window != None:
            self.window.close()
        self.window = PrintoutWindow(self)
        self.window.setWindowTitle("Solidworks Export Macro")

        settings_fetched = self.parent().parent().parent().get_all_settings()
        if settings_fetched == None:
            return
        data = create_3d_blade(settings_fetched, self.flip_turning_direction.isChecked(),
                               self.propeller_geom.isChecked())

        # DRAW MATPLOTLIB ############################
        self.w = MatplotlibWindow()
        self.w.setWindowTitle("Export 3D preview")
        self.w.ax = self.w.figure.add_subplot(111, projection="3d")
        self.w.ax.scatter(data["X"], data["Y"], data["Z"])
        X, Y, Z = array(data["X"]), array(data["Y"]), array(data["Z"])

        # Create cubic bounding box to simulate equal aspect ratio
        max_range = np.array([X.max() - X.min(), Y.max() - Y.min(), Z.max() - Z.min()]).max()
        Xb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][0].flatten() + 0.5 * (X.max() + X.min())
        Yb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][1].flatten() + 0.5 * (Y.max() + Y.min())
        Zb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][2].flatten() + 0.5 * (Z.max() + Z.min())
        # Comment or uncomment following both lines to test the fake bounding box:
        for xb, yb, zb in zip(Xb, Yb, Zb):
            self.w.ax.plot([xb], [yb], [zb], 'w')
        # self.w.ax.set_aspect("equal")
        ##############################################

        # Create necessary folders
        create_folder(os.path.join(application_path, "export"))
        folder_path = os.path.join(application_path, "export", SET_INIT["turbine_name"])
        create_folder(folder_path)


        #print("Exporting... (and converting m to mm)")

        filenames = []
        data_out = []
        for z, x_data, y_data in data["data"]:
            #print("Exporting z=" + str(z), "[m]")
            #z = z * 1e3  # in mm
            file_name = os.path.join(folder_path, "z_%s.txt" % z)
            filenames.append(os.path.join(os.getcwd(), file_name))
            # print(file_name)
            f = open(os.path.join(folder_path, "z_%s.txt" % z), "w")
            for x, y in zip(x_data, y_data):
                x_out, y_out, z_out = x * 1e3, y * 1e3, z * 1e3  # in mm
                f.write("%s\t%s\t%s\n" % (y_out, z_out, x_out)) # Change indexes for SW
            f.close()

            data_out.append([y_data,np.zeros(len(y_data))+z,x_data])

        #print("Filenames:", filenames)
        macro_text = create_macro_text(filenames,data_out)

        print("'===============MACRO START==================")
        print(macro_text)
        print("'===============MACRO   END==================")