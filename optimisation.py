__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Miha Smrekar"
__email__ = "miha.smrekar9@gmail.com"
__status__ = "Development"

import numbers

import numpy

from cp_curve import calculate_power
from parse_data import parse_sections
from parse_data import f_c_L, f_c_D


class Optimizer:
    """
    Class used for optimizing blade geometry.
    TODO: add input arguments to init, so that user will have to pass geometry data and functions for initialisation.
    """
    def __init__(self):
        sections_radius, chord_lengths, chord_angles, dr = parse_sections()

        self.sections_radius = sections_radius
        self.chord_lengths = chord_lengths
        self.chord_angles = chord_angles
        self.dr = dr

        self.chord_angles_orig = numpy.copy(chord_angles)
        self.target_speed = 7
        self.target_rpm = 500
        self.R = .776
        self.Rhub = 0.1
        self.B = 5

        self.f_c_L = f_c_L
        self.f_c_D = f_c_D

    def _power(self, _r=None, _delta=None):
        """
        Helper function that calculates power.

        If _r and _delta are provided, returns power with changing the blade section number [_r] for [_delta] degrees.

        If none of the arguments are provided, just calculates power.

        :param _r: section number (int)
        :param _delta: change of twist - theta - for given section [deg]
        :return: power [W]
        """
        if _r and _delta:
            chord_angles_copy = numpy.copy(self.chord_angles)
            chord_angles_copy[_r] += _delta
            return calculate_power(
                speed_wind=self.target_speed, rpm=self.target_rpm, sections_radius=self.sections_radius,
                chord_lengths=self.chord_lengths,
                chord_angles=chord_angles_copy, dr=self.dr, R=self.R, B=self.B, f_c_L=self.f_c_L,
                f_c_D=self.f_c_D,
                Rhub=self.Rhub)["power"]
        return calculate_power(
            speed_wind=self.target_speed, rpm=self.target_rpm, sections_radius=self.sections_radius,
            chord_lengths=self.chord_lengths,
            chord_angles=self.chord_angles, dr=self.dr, R=self.R, B=self.B, f_c_L=self.f_c_L,
            f_c_D=self.f_c_D,
            Rhub=self.Rhub)["power"]

    def direction(self, _r=None, _delta=None):
        """
        Determines direction of changing the blade section number [_r] twist for an angle [_delta],
        that will provide better power than the initial twist angle.

        If better power is gained through increasing twist (twist=twist+_delta), returns True.
        If better power is gained through decreasing twist (twist=twist-_delta), returns False.
        If better power cannot be provided through increasing or decreasing, returns None.

        :param _r: blade section number (int)
        :param _delta: change of twist - theta - for given section [deg]
        :return: True, False or None
        """
        power_orig = self._power()
        if _delta <= 0.0:
            raise Exception("Delta has to be a real number greater than 0.")
        if isinstance(_r, int) and isinstance(_delta, numbers.Real):
            power_up = self._power(_r, _delta)
            power_down = self._power(_r, -_delta)
            if power_up > power_down:
                return True
            if power_down > power_up:
                return False
            if power_up == power_down:
                if power_up == power_orig and power_down == power_orig:
                    return None

        else:
            raise Exception("r and delta have to be provided.")

    def optimize_angles(self):
        """
        Optimizes angles of attack by changing twist angles.
        Prints results.
        """
        power_preliminary = self._power()

        print("Optimizing angles for target speed", self.target_speed,
              "m/s, target rpm", self.target_rpm,
              ", number of blades:", 5,
              ", tip radius:", self.R)

        MINDELTA = 0.01

        print("power at start is", power_preliminary)

        for r in range(len(self.sections_radius)):
            print("Currently optimizing section", r)
            power = self._power()
            delta = 5.0
            current_delta = 0.0

            while True:
                print("    Current delta", delta)
                if delta < MINDELTA:
                    print("    Delta limit reached")
                    break

                d = self.direction(r, delta)
                run_optimize = False

                if d == True:
                    print("    Direction is up")
                    current_delta = delta
                    run_optimize = True
                if d == False:
                    print("    Direciton is down")
                    current_delta = -delta
                    run_optimize = True
                if d == None:
                    print("    No direction provides better results.")

                if run_optimize:
                    power_old = power
                    while True:
                        print("        Changing blade section", r, "angle from", self.chord_angles[r], "to",
                              self.chord_angles[r] + current_delta)
                        self.chord_angles[r] += current_delta
                        power = self._power()
                        if power <= power_old:
                            print("        New power", power, "is less than or equal to old power", power_old)
                            break
                        print("        New power", power, "is greater than old power", power_old)
                        power_old = power
                    print("    Changing delta from", delta, "to", delta * 0.1)

                delta = delta * 0.1

        power_new = self._power()

        print("power at start was", power_preliminary)
        print("power at end is", power_new)
        print("Power increase", power_new - power_preliminary)
        percentage_increase = power_new / power_preliminary * 100
        print("Percentage increase", percentage_increase)
        print("Old angles", self.chord_angles_orig)
        print("New angles", self.chord_angles)

    
    def optimize_pitch(self, min_add_angle=-30, max_add_angle=30, step=0.5):
        """
        This function calculates the optimum pitch angle of the blade for the given wind speed and rotational velocity.

        :param min_add_angle: lower limit [degrees]. Default: -30
        :param max_add_angle: upper limit [degrees]. Default: +30
        :param step: step for changing angle [degrees]
        :return: best angle change [degrees]
        """
        print("\nOptimizing for wind speed of ", self.target_speed, "m/s... and rpm ", self.target_rpm)
        power_orig = \
            calculate_power(speed_wind=self.target_speed, rpm=self.target_rpm, sections_radius=self.sections_radius,
                chord_lengths=self.chord_lengths,
                chord_angles=self.chord_angles_orig, dr=self.dr, R=self.R, B=self.B, f_c_L=self.f_c_L,
                f_c_D=self.f_c_D, Rhub=self.Rhub)["power"]
        print("Without turning the blade, the power is:", "%.2f" % round(power_orig), "Watts at wind speed", self.target_speed,
              "m/s and rotational velocity", self.target_rpm, "rpm")
        current_add_angle = min_add_angle
        results = []
        while current_add_angle <= max_add_angle:
            print("Testing pitch:", current_add_angle, "° at rpm", self.target_rpm)
            power = \
                calculate_power(speed_wind=self.target_speed, rpm=self.target_rpm, sections_radius=self.sections_radius,
                chord_lengths=self.chord_lengths,
                chord_angles=self.chord_angles_orig, dr=self.dr, R=self.R, B=self.B, f_c_L=self.f_c_L,
                f_c_D=self.f_c_D, Rhub=self.Rhub,add_angle=current_add_angle)["power"]
            results.append((current_add_angle, power))
            current_add_angle += step
            print("---")
        best_angle, best_power = sorted(results, key=lambda x: x[1], reverse=True)[0]
        print("The best angle to turn would be:", "%.3f" % round(best_angle, 2), "°, at this pitch angle the power is:",
              "%.2f" % round(best_power), "Watts")
        print("That is a power gain of ", "%.2f" % round(best_power - power_orig, 2), "Watts, or",
              round(((best_power - power_orig) / power_orig) * 100, 2), "% better than the original version.")
        return best_angle


#Optimizer().optimize_pitch()
Optimizer().optimize_angles()