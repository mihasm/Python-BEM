__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Miha Smrekar"
__email__ = "miha.smrekar9@gmail.com"
__status__ = "Development"

import numpy

from induction_factors import Calculator


# from cp_curve import calculate_power
def calculate_power(
    speed_wind,
    rpm,
    sections_radius,
    chord_lengths,
    chord_angles,
    dr,
    R,
    Rhub,
    B,
    f_c_L,
    f_c_D,
    add_angle=None,
):
    """
    Returns calculated power using BEM analysis.

    Inputs are wind speed, rotational velocity, blade geometry, number of blades, and
    functions for calculating lift and drag coefficients.

    Output is a dictionary with all results.

    :param speed_wind: wind speed [m/s]
    :param rpm: rotational velocity [rpm]
    :param sections_radius: np array of section sections_radius [m]
    :param chord_lengths: np array of section chord lengths [m]
    :param chord_angles: np array of chord angles [degrees]
    :param dr: np array of section heights [m]
    :param R: outer (tip) radius [m]
    :param Rhub: hub radius [m]
    :param B: number of blades
    :param f_c_L: function for calculating lift coefficient
    :param f_c_D: function for calculating drag coefficient
    :param add_angle: [degrees]
    :return: dict with results
    """
    if add_angle != None:
        chord_angles = chord_angles + add_angle
    results = Calculator(f_c_L, f_c_D).run_array(
        chord_angle=chord_angles,
        B=B,
        c=chord_lengths,
        r=sections_radius,
        dr=dr,
        rpm=rpm,
        v=speed_wind,
        R=R,
        Rhub=Rhub,
    )
    return results


class Optimizer:
    """
    Class used for optimizing blade geometry.
    TODO: add input arguments to init, so that user will have to pass geometry data and functions for initialisation.
    """

    def __init__(self, r, c, theta, dr, B, R, Rhub, f_c_L, f_c_D):
        # sections_radius, chord_lengths, chord_angles, dr = parse_sections()

        self.sections_radius = r
        self.chord_lengths = c
        self.chord_angles = theta
        self.dr = dr

        self.chord_angles_orig = numpy.copy(self.chord_angles)

        self.target_speed = 7
        self.target_rpm = 500

        self.R = R
        self.Rhub = Rhub
        self.B = B

        self.f_c_L = f_c_L
        self.f_c_D = f_c_D

    def _power(
        self,
        _r=None,
        _delta=None,
        add_angle=None,
        wind_speed=None,
        rpm=None,
        power_only=True,
    ):
        """
        Helper function that calculates power.

        If _r and _delta are provided, returns power with changing the blade section number [_r] for [_delta] degrees.

        If _r and _delta are not provided, just calculates power.

        :param _r: section number (int)
        :param _delta: change of twist - theta - for given section [deg]
        :return: power [W]
        """
        if _r and _delta:
            chord_angles_copy = numpy.copy(self.chord_angles)
            chord_angles_copy[_r] += _delta
            res = calculate_power(
                speed_wind=wind_speed,
                rpm=rpm,
                sections_radius=self.sections_radius,
                chord_lengths=self.chord_lengths,
                chord_angles=chord_angles_copy,
                dr=self.dr,
                R=self.R,
                B=self.B,
                f_c_L=self.f_c_L,
                f_c_D=self.f_c_D,
                add_angle=add_angle,
                Rhub=self.Rhub,
            )
        else:
            res = calculate_power(
                speed_wind=wind_speed,
                rpm=rpm,
                sections_radius=self.sections_radius,
                chord_lengths=self.chord_lengths,
                chord_angles=self.chord_angles,
                dr=self.dr,
                R=self.R,
                B=self.B,
                f_c_L=self.f_c_L,
                f_c_D=self.f_c_D,
                add_angle=add_angle,
                Rhub=self.Rhub,
            )

        if power_only:
            return res["power"]
        else:
            return res

    def direction(self, _r, _delta, wind_speed, rpm):
        """
        Determines direction of changing the blade section number [_r] twist for an angle [_delta],
        that will provide better power than the initial twist angle.

        If better power is gained through increasing twist (twist=twist+_delta), returns True.
        If better power is gained through decreasing twist (twist=twist-_delta), returns False.
        If better power cannot be provided through increasing or decreasing, returns None.

        :param rpm: rotational velocity [RPM]
        :param wind_speed: wind speed [m/s]
        :param _r: blade section number (int)
        :param _delta: change of twist - theta - for given section [deg]
        :return: True, False or None
        """
        power_orig = self._power(wind_speed=wind_speed, rpm=rpm)
        power_up = self._power(_r, _delta, wind_speed=wind_speed, rpm=rpm)
        power_down = self._power(_r, -_delta, wind_speed=wind_speed, rpm=rpm)
        if power_up > power_down:
            return True
        if power_down > power_up:
            return False
        if power_up == power_down:
            if power_up == power_orig and power_down == power_orig:
                return None

    def optimize_angles(self, target_speed, target_rpm):
        """
        Optimizes angles of attack by changing twist angles.
        Prints results.
        """
        power_preliminary = self._power(wind_speed=target_speed, rpm=target_rpm)

        print(
            "Optimizing angles for target speed",
            target_speed,
            "m/s, target rpm",
            target_rpm,
            ", number of blades:",
            5,
            ", tip radius:",
            self.R,
        )

        MINDELTA = 0.01

        print("power at start is", power_preliminary)

        for r in range(len(self.sections_radius)):
            print("Currently optimizing section", r)
            power = self._power(wind_speed=target_speed, rpm=target_rpm)
            delta = 5.0
            current_delta = 0.0

            while True:
                print("    Current delta", delta)
                if delta < MINDELTA:
                    print("    Delta limit reached")
                    break

                d = self.direction(r, delta, wind_speed=target_speed, rpm=target_rpm)
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
                        print(
                            "        Changing blade section",
                            r,
                            "angle from",
                            self.chord_angles[r],
                            "to",
                            self.chord_angles[r] + current_delta,
                        )
                        self.chord_angles[r] += current_delta
                        power = self._power(wind_speed=target_speed, rpm=target_rpm)
                        if power <= power_old:
                            print(
                                "        New power",
                                power,
                                "is less than or equal to old power",
                                power_old,
                            )
                            break
                        print(
                            "        New power",
                            power,
                            "is greater than old power",
                            power_old,
                        )
                        power_old = power
                    print("    Changing delta from", delta, "to", delta * 0.1)

                delta = delta * 0.1

        power_new = self._power(wind_speed=target_speed, rpm=target_rpm)

        print("power at start was", power_preliminary)
        print("power at end is", power_new)
        print("Power increase", power_new - power_preliminary)
        percentage_increase = power_new / power_preliminary * 100
        print("Percentage increase", percentage_increase)
        old_angles = self.chord_angles_orig
        new_angles = self.chord_angles
        print("Old angles", old_angles)
        print("New angles", new_angles)

        # reset chord angles so other functions work properly
        self.chord_angles = numpy.copy(self.chord_angles_orig)
        return new_angles

    def optimize_pitch(
        self,
        target_speed,
        target_rpm=None,
        min_add_angle=-10,
        max_add_angle=10,
        step=0.5,
    ):
        """
        This function calculates the optimum pitch angle of the blade for the given wind speed and rotational velocity.

        :param target_rpm: rotational velocity [RPM]
        :param target_speed: wind speed [m/s]
        :param min_add_angle: lower limit [degrees]. Default: -30
        :param max_add_angle: upper limit [degrees]. Default: +30
        :param step: step for changing angle [degrees]
        :return: best angle change [degrees]
        """
        if target_rpm != None:
            print(
                "\nOptimizing for wind speed of ",
                target_speed,
                "m/s... and rpm ",
                target_rpm,
            )

            power_orig = self._power(wind_speed=target_speed, rpm=target_rpm)
            print(
                "Without turning the blade, the power is:",
                "%.2f" % round(power_orig),
                "Watts at wind speed",
                target_speed,
                "m/s and rotational velocity",
                target_rpm,
                "rpm",
            )

            current_add_angle = min_add_angle
            results = []
            while current_add_angle <= max_add_angle:
                print("Testing pitch:", current_add_angle, "° at rpm", target_rpm)
                power = self._power(
                    add_angle=current_add_angle, wind_speed=target_speed, rpm=target_rpm
                )
                results.append((current_add_angle, power))
                current_add_angle += step
                print("---")

            best_angle, best_power = sorted(results, key=lambda x: x[1], reverse=True)[
                0
            ]

            print(
                "The best angle to turn would be:",
                "%.3f" % round(best_angle, 2),
                "°, at this pitch angle the power is:",
                "%.2f" % round(best_power),
                "Watts",
            )
            print(
                "That is a power gain of ",
                "%.2f" % round(best_power - power_orig, 2),
                "Watts, or",
                round(((best_power - power_orig) / power_orig) * 100, 2),
                "% better than the original version.",
            )

            return best_angle
        else:
            print("Getting best power with original pitch.")
            power_orig = (None, None)  # power, rpm
            for rpm in numpy.linspace(start=100, stop=3000, num=30):
                res = self._power(wind_speed=target_speed, rpm=rpm, power_only=False)
                if 0.0 < res["cp"] <= 0.6:
                    p = res["power"]
                    print("    rpm:", rpm, "p:", p)
                    if power_orig[0] == None and power_orig[1] == None:
                        power_orig = (p, rpm)
                    elif p > power_orig[0]:
                        power_orig = (p, rpm)
            print("Original power", power_orig)
            current_add_angle = min_add_angle
            results = []
            print("Getting power for different pitches")
            while current_add_angle <= max_add_angle:
                print("    Current add angle:", current_add_angle)
                best_result_angle = (None, None, None)
                for rpm in numpy.linspace(start=100, stop=3000, num=30):
                    res = self._power(
                        wind_speed=target_speed,
                        rpm=rpm,
                        add_angle=current_add_angle,
                        power_only=False,
                    )
                    if 0.0 < res["cp"] <= 0.6:
                        p = res["power"]
                        print("        rpm:", rpm, "p:", p)
                        if (
                            best_result_angle[0] == None
                            and best_result_angle[1] == None
                        ):
                            best_result_angle = (p, rpm, current_add_angle)
                        if p > best_result_angle[0]:
                            best_result_angle = (p, rpm, current_add_angle)
                    else:
                        p = res["power"]
                        print("        rpm:", rpm, "p:", p, "not possible / bad result")
                results.append(best_result_angle)
                current_add_angle += step
            best = sorted(results, key=lambda x: x[0], reverse=True)[0]
            print("Best power is", best[0], "at rpm", best[1], "using pitch", best[2])
            print("That is a power gain of", best[0] - power_orig[0], "Watts.")
        best_chord_angles = self.chord_angles_orig + best[2]
        return best_chord_angles
