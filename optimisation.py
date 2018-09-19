__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.2.6"
__maintainer__ = "Miha Smrekar"
__email__ = "miha.smrekar9@gmail.com"
__status__ = "Development"

import numpy

from induction_factors import Calculator


from cp_curve import calculate_power


class Optimizer:
    """
    Class used for optimizing blade geometry.
    TODO: add input arguments to init, so that user will have to pass geometry data and functions for initialisation.
    """

    def __init__(self,inp_args):
        # sections_radius, chord_lengths, chord_angles, dr = parse_sections()
        self.inp_args = inp_args

    def _power(
        self,
        wind_speed,
        rpm,
        power_only=True,
         _r=None,
        _delta=None,
        add_angle=None,
    ):
        """
        Helper function that calculates power.

        If _r and _delta are provided, returns power with changing the blade section number [_r] for [_delta] degrees.

        If _r and _delta are not provided, just calculates power.

        :param _r: section number (int)
        :param _delta: change of twist - theta - for given section [deg]
        :return: power [W]
        """
        self.inp_args["v"] = wind_speed
        self.inp_args["rpm"] = rpm
        if _r and _delta:
            chord_angles_copy = numpy.copy(self.inp_args["theta"])
            chord_angles_copy[_r] += _delta
            _inp_args = {**self.inp_args}
            _inp_args["theta"] = chord_angles_copy
            res = calculate_power(
                _inp_args
            )
        else:
            res = calculate_power(
                self.inp_args
            )

        if power_only:
            return res["power"]
        else:
            return res

    def direction(self, _r, _delta, _wind_speed, _rpm, *args, **kwargs):
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
        power_orig = self._power(wind_speed=_wind_speed, rpm=_rpm)
        power_up = self._power(_r=_r, _delta=_delta, wind_speed=_wind_speed, rpm=_rpm)
        power_down = self._power(_r=_r, _delta=-_delta, wind_speed=_wind_speed, rpm=_rpm)
        if power_up > power_down:
            return True
        if power_down > power_up:
            return False
        if power_up == power_down:
            if power_up == power_orig and power_down == power_orig:
                return None

    def optimize_angles(
        self,
        inp_args,
    ):
        """
        Optimizes angles of attack by changing twist angles.
        Prints results.
        """
        power_preliminary = self._power(wind_speed=inp_args["target_speed"], rpm=inp_args["target_rpm"])
        return_print = inp_args["return_print"]
        return_results = inp_args["return_results"]
        return_print.append(
            "Optimizing angles for target speed "
            + str(inp_args["target_speed"])
            + " [m/s], target rpm "
            + str(inp_args["target_rpm"])
            + " [RPM], number of blades:"
            + str(5)
            + ", tip radius:"
            + str(inp_args["R"])
            + " [m]\n"
            ", hub radius:" + str(inp_args["Rhub"]) + " [m]\n"
        )

        return_print.append("power at start is " + str(power_preliminary) + "\n")

        for r in range(len(inp_args["r"])):
            return_print.append("Currently optimizing section " + str(r) + "\n")
            power = self._power(wind_speed=inp_args["target_speed"], rpm=inp_args["target_rpm"])
            delta = inp_args["delta_start"]
            current_delta = 0.0

            while True:
                return_print.append("    Current delta " + str(delta) + "\n")
                if delta < inp_args["min_delta"]:
                    return_print.append("    Delta limit reached\n")
                    break

                d = self.direction(r, delta, _wind_speed=inp_args["target_speed"], _rpm=inp_args["target_rpm"])
                run_optimize = False

                if d == True:
                    return_print.append("    Direction is up\n")
                    current_delta = delta
                    run_optimize = True
                if d == False:
                    return_print.append("    Direction is down\n")
                    current_delta = -delta
                    run_optimize = True
                if d == None:
                    return_print.append("    No direction provides better results.\n")

                if run_optimize:
                    power_old = power
                    while True:
                        return_print.append(
                            "        Changing blade section "
                            + str(r)
                            + " angle from "
                            + str(inp_args["theta"][r])
                            + " to "
                            + str(inp_args["theta"][r]+current_delta)
                            + "\n"
                        )
                        inp_args["theta"][r] += current_delta
                        power = self._power(wind_speed=inp_args["target_speed"], rpm=inp_args["target_rpm"])
                        if power <= power_old:
                            return_print.append(
                                "        New power "
                                + str(power)
                                + " is less than or equal to old power "
                                + str(power_old)
                                + "\n"
                            )
                            break
                        return_print.append(
                            "        New power "
                            + str(power)
                            + " is greater than old power "
                            + str(power_old)
                            + "\n"
                        )
                        power_old = power
                    return_print.append(
                        "    Changing delta from "
                        + str(delta)
                        + " to "
                        + str(delta * inp_args["decrease_factor"])
                        + "\n"
                    )

                delta = delta * inp_args["decrease_factor"]

        power_new = self._power(wind_speed=inp_args["target_speed"], rpm=inp_args["target_rpm"])

        return_print.append("power at start was " + str(power_preliminary) + "\n")
        return_print.append("power at end is " + str(power_new) + "\n")
        return_print.append(
            "Power increase " + str(power_new - power_preliminary) + "\n"
        )
        percentage_increase = power_new / power_preliminary * 100
        return_print.append("Percentage " + str(percentage_increase) + "\n")
        old_angles = self.chord_angles_orig
        new_angles = self.chord_angles
        return_print.append("Old angles " + str(old_angles) + "\n")
        return_print.append("New angles " + str(new_angles) + "\n")
        return_print.append("----------------------" + "\n")

        # reset chord angles so other functions work properly
        self.chord_angles = numpy.copy(self.chord_angles_orig)

        return_print.append("!!!!EOF!!!!")
        return new_angles

    def optimize_pitch(
        self,
        target_speed,
        target_rpm=None,
        min_add_angle=-10,
        max_add_angle=10,
        step=0.5,
        return_print=None,
        return_results=None,
        *args,
        **kwargs,
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
        return_print = inp_args["return_print"]
        return_results = inp_args["return_results"]

        if inp_args["target_rpm"] != None:
            return_print.append(
                "\nOptimizing for wind speed of "
                + str(inp_args["target_speed"])
                + "m/s... and rpm "
                + str(inp_args["target_rpm"])
                + " RPM\n"
            )

            power_orig = self._power(wind_speed=inp_args["target_speed"], rpm=inp_args["target_rpm"])
            return_print.append(
                "Without turning the blade, the power is: "
                + "%.2f" % round(power_orig)
                + " [W] at wind speed "
                + str(inp_args["target_speed"])
                + "[m/s] and rotational velocity "
                + str(inp_args["target_rpm"])
                + "[RPM]"
                + "\n"
            )

            current_add_angle = min_add_angle
            results = []
            while current_add_angle <= max_add_angle:
                return_print.append(
                    "Testing pitch: "
                    + str(current_add_angle)
                    + "° at rpm "
                    + str(inp_args["target_rpm"])
                    + " [RPM]\n"
                )
                power = self._power(
                    add_angle=current_add_angle, wind_speed=inp_args["target_speed"], rpm=inp_args["target_rpm"]
                )
                results.append((current_add_angle, power))
                current_add_angle += step
                return_print.append("    power:" + str(power) + "\n")
                return_print.append("---\n")

            best_angle, best_power = sorted(results, key=lambda x: x[1], reverse=True)[
                0
            ]

            return_print.append(
                "The best angle to turn would be:"
                + "%.3f" % round(best_angle, 2)
                + "°, at this pitch angle the power is:"
                + "%.2f" % round(best_power)
                + "Watts"
                + "\n"
            )
            return_print.append(
                "That is a power gain of "
                + "%.2f" % round(best_power - power_orig, 2)
                + "Watts"
                + "\n"
            )
            return_print.append("!!!!EOF!!!!")
            return best_angle
        else:
            return_print.append("Getting best power with original pitch.\n")
            power_orig = (None, None)  # power, rpm
            for rpm in numpy.linspace(start=100, stop=3000, num=30):
                res = self._power(wind_speed=inp_args["target_speed"], rpm=rpm, power_only=False)
                if 0.0 < res["cp"] <= 0.6:
                    p = res["power"]
                    return_print.append("    rpm:" + str(rpm) + " p:" + str(p) + "\n")
                    if power_orig[0] == None and power_orig[1] == None:
                        power_orig = (p, rpm)
                    elif p > power_orig[0]:
                        power_orig = (p, rpm)
            return_print.append("Original power" + str(power_orig) + "\n")
            current_add_angle = min_add_angle
            results = []
            return_print.append("Getting power for different pitches" + "\n")
            while current_add_angle <= max_add_angle:
                return_print.append(
                    "    Current add angle:" + str(current_add_angle) + "\n"
                )
                best_result_angle = (None, None, None)
                for rpm in numpy.linspace(start=100, stop=3000, num=30):
                    res = self._power(
                        wind_speed=inp_args["target_speed"],
                        rpm=rpm,
                        add_angle=current_add_angle,
                        power_only=False,
                    )
                    if 0.0 < res["cp"] <= 0.6:
                        p = res["power"]
                        return_print.append(
                            "        rpm:" + str(rpm) + "p:" + str(p) + "\n"
                        )
                        if (
                            best_result_angle[0] == None
                            and best_result_angle[1] == None
                        ):
                            best_result_angle = (p, rpm, current_add_angle)
                        if p > best_result_angle[0]:
                            best_result_angle = (p, rpm, current_add_angle)
                    else:
                        p = res["power"]
                        return_print.append(
                            "        rpm:"
                            + str(rpm)
                            + " p:"
                            + str(p)
                            + " not possible / bad result"
                            + "\n"
                        )
                results.append(best_result_angle)
                current_add_angle += step
            best = sorted(results, key=lambda x: x[0], reverse=True)[0]
            return_print.append(
                "Best power is "
                + str(best[0])
                + " at rpm "
                + str(best[1])
                + " using pitch "
                + str(best[2])
                + "°\n"
            )
            return_print.append(
                "That is a power gain of "
                + str(best[0] - power_orig[0])
                + " Watts."
                + "\n"
            )
        best_chord_angles = self.chord_angles_orig + best[2]
        return_print.append("!!!!EOF!!!!")
        return best_chord_angles
