__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.2.8"
__maintainer__ = "Miha Smrekar"
__email__ = "miha.smrekar9@gmail.com"
__status__ = "Development"

import numpy

from induction_factors import Calculator
from cp_curve import calculate_power
from utils import Printer


class Optimizer:
    """
    Class used for optimizing blade geometry.
    TODO: add input arguments to init, so that user will have to pass geometry data and functions for initialisation.
    """

    def __init__(self, inp_args):
        # sections_radius, chord_lengths, chord_angles, dr = parse_sections()
        self.inp_args = inp_args

    def _power(
            self, wind_speed, rpm, power_only=True, _r=None, _delta=None, add_angle=None
    ):
        """
        Helper function that calculates power.

        If _r and _delta are provided, returns power with changing the blade section number [_r] for [_delta] degrees.

        If _r and _delta are not provided, just calculates power.

        :param _r: section number (int)
        :param _delta: change of twist - theta - for given section [deg]
        :return: power [W]
        """

        # creates a copy of the dict so that input dict is not affected
        _inp_args = {**self.inp_args}
        _inp_args["v"] = wind_speed
        _inp_args["rpm"] = rpm

        if add_angle != None:
            _inp_args["add_angle"] = add_angle

        if _r and _delta:
            chord_angles_copy = numpy.copy(self.inp_args["theta"])
            chord_angles_copy[_r] += _delta
            _inp_args["theta"] = chord_angles_copy
            res = calculate_power(_inp_args)
        else:
            res = calculate_power(_inp_args)

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

        :param _rpm:
        :param _wind_speed:
        :param _r: blade section number (int)
        :param _delta: change of twist - theta - for given section [deg]
        :return: True, False or None
        """

        power_orig = self._power(wind_speed=_wind_speed, rpm=_rpm)
        power_up = self._power(_r=_r, _delta=_delta,
                               wind_speed=_wind_speed, rpm=_rpm)
        power_down = self._power(
            _r=_r, _delta=-_delta, wind_speed=_wind_speed, rpm=_rpm
        )
        if power_up > power_down:
            return True
        if power_down > power_up:
            return False
        if power_up == power_down:
            if power_up == power_orig and power_down == power_orig:
                return None

    def optimize_angles(self, inp_args):
        """
        Optimizes angles of attack by changing twist angles.
        Prints results.
        """
        p = Printer(inp_args["return_print"])
        return_results = inp_args["return_results"]

        # creates copy of original chord angles
        chord_angles_orig = numpy.copy(inp_args["theta"])
        chord_angles = numpy.copy(inp_args["theta"])

        # init power
        power_preliminary = self._power(
            wind_speed=inp_args["target_speed"],
            rpm=inp_args["target_rpm"]
        )

        # printout
        p.print("---- Optimizing angles ----")
        p.print("Wind speed:", inp_args["target_speed"], "[m/s]")
        p.print("       Rpm:", inp_args["target_rpm"], "[RPM]")
        p.print("    Blades:", 5)
        p.print("Tip radius:", inp_args["R"], "[m]")
        p.print("Hub radius:", inp_args["Rhub"], "[m]")
        p.print("    Angles:", chord_angles_orig)
        p.print("---------------------------")
        p.print("Power at start:", power_preliminary, "[W]")

        # for every blade section
        for r in range(len(inp_args["r"])):

            p.print("Currently optimizing section", r)

            # initial power for comparison
            power = self._power(
                wind_speed=inp_args["target_speed"],
                rpm=inp_args["target_rpm"]
            )

            # grab delta
            delta = inp_args["delta_start"]

            while True:

                p.print("    Current delta", delta)

                # check if delta is smaller than limit
                if delta < inp_args["min_delta"]:
                    p.print("    Delta limit reached")
                    break

                # get direction
                d = self.direction(
                    r,
                    delta,
                    _wind_speed=inp_args["target_speed"],
                    _rpm=inp_args["target_rpm"],
                )

                # Decide whether to optimize (if any direction is better)
                run_optimize = False

                if d == True:
                    p.print("    Direction is up")
                    delta = +abs(delta)
                    run_optimize = True
                if d == False:
                    p.print("    Direction is down")
                    delta = -abs(delta)
                    run_optimize = True
                if d == None:
                    p.print("    No direction provides better results.")

                if run_optimize:

                    # save old value of power
                    power_old = power

                    while True:

                        p.print("        Changing section", r, "angle from",
                                inp_args["theta"][r],
                                "to",
                                inp_args["theta"][r] + delta
                                )

                        # increase section r angle by delta
                        inp_args["theta"][r] += delta

                        power = self._power(
                            wind_speed=inp_args["target_speed"],
                            rpm=inp_args["target_rpm"],
                        )

                        # check whether goal reached
                        if power <= power_old:
                            p.print("        New power", power,
                                    "<= old power", power_old)
                            break
                        else:
                            chord_angles[r] = inp_args["theta"][r]

                            p.print("        New power", power,
                                    "> old power", power_old)
                            power_old = power

                p.print("    Changing delta from", abs(delta), "to",
                        abs(delta) * inp_args["decrease_factor"])

                # decrease delta by given amount
                delta = abs(delta) * inp_args["decrease_factor"]

        power_new = self._power(
            wind_speed=inp_args["target_speed"], rpm=inp_args["target_rpm"]
        )

        p.print("power at start was", power_preliminary)
        p.print("power at end is", power_new)
        p.print(
            "Power increase", power_new - power_preliminary
        )
        percentage_increase = power_new / power_preliminary * 100
        p.print("Percentage",percentage_increase)
        old_angles = chord_angles_orig
        new_angles = chord_angles
        p.print("Old angles",old_angles)
        p.print("New angles",new_angles)
        p.print("----------------------")

        # reset chord angles so other functions work properly
        # self.chord_angles = numpy.copy(self.chord_angles_orig)
        inp_args["theta"] = chord_angles_orig
        p.print("!!!!EOF!!!!")
        return new_angles

    def optimize_pitch(
            self, inp_args
    ):
        """
        This function calculates the optimum pitch angle of the blade for the given wind speed and rotational velocity.

        :param inp_args:
        :return: best angle change [degrees]
        """
        p = Printer(inp_args["return_print"])
        return_results = inp_args["return_results"]
        chord_angles_orig = numpy.copy(inp_args["theta"])

        if inp_args["target_rpm"] != None:
            p.print(
                "Optimizing for wind speed of", inp_args[
                    "target_speed"], "m/s... and rpm", inp_args["target_rpm"], "RPM"
            )

            power_orig = self._power(
                wind_speed=inp_args["target_speed"], rpm=inp_args["target_rpm"]
            )

            p.print(
                "Without turning the blade, the power is: ", "%.2f" % round(power_orig), "[W] at wind speed", inp_args[
                    "target_speed"], "[m/s] and rotational velocity", inp_args["target_rpm"], "[RPM]"

            )

            current_add_angle = inp_args["min_add_angle"]
            results = []

            while current_add_angle <= inp_args["max_add_angle"]:
                p.print(
                    "Testing pitch: ", current_add_angle, "° at rpm", inp_args["target_rpm"], "[RPM]"
                )

                power_new = self._power(
                    add_angle=current_add_angle,
                    wind_speed=inp_args["target_speed"],
                    rpm=inp_args["target_rpm"],
                )

                results.append((current_add_angle, power_new))
                p.print("    power:", power_new)
                p.print("---")

                current_add_angle += inp_args["angle_step"]

            best_angle, best_power = sorted(
                results, key=lambda x: x[1], reverse=True)[0]

            p.print(
                "The best angle to turn would be:", "%.3f" % round(
                    best_angle, 2), "°, at this pitch angle the power is:", "%.2f" % round(best_power), "Watts"

            )
            p.print(
                "That is a power gain of", "%.2f" % round(
                    best_power - power_orig, 2), "Watts"

            )
            p.print("!!!!EOF!!!!")
            return best_angle
        else:
            p.print("Getting best power with original pitch.")
            power_orig = (None, None)  # power, rpm
            for rpm in numpy.linspace(start=100, stop=3000, num=30):
                res = self._power(
                    wind_speed=inp_args[
                        "target_speed"], rpm=rpm, power_only=False
                )
                if 0.0 < res["cp"] <= 0.6:
                    pwr = res["power"]
                    p.print("    rpm:", rpm, "p:", pwr)
                    if power_orig[0] == None and power_orig[1] == None:
                        power_orig = (pwr, rpm)
                    elif pwr > power_orig[0]:
                        power_orig = (pwr, rpm)
            p.print("Original power", power_orig)
            current_add_angle = inp_args["min_add_angle"]
            results = []
            p.print("Getting power for different pitches")

            while current_add_angle <= inp_args["max_add_angle"]:
                p.print("    Current add angle:",current_add_angle)
                best_result_angle = (None, None, None)
                for rpm in numpy.linspace(start=100, stop=3000, num=30):
                    res = self._power(
                        wind_speed=inp_args["target_speed"],
                        rpm=rpm,
                        add_angle=current_add_angle,
                        power_only=False,
                    )
                    if 0.0 < res["cp"] <= 0.6:
                        pwr = res["power"]
                        p.print("        rpm:",rpm,"p:",str(pwr))
                        if (
                                best_result_angle[0] == None
                                and best_result_angle[1] == None
                        ):
                            best_result_angle = (pwr, rpm, current_add_angle)
                        if pwr > best_result_angle[0]:
                            best_result_angle = (pwr, rpm, current_add_angle)
                    else:
                        pwr = res["power"]
                        p.print(
                            "        rpm:",rpm, "p:",pwr, "not possible / bad result"

                        )
                results.append(best_result_angle)
                current_add_angle += inp_args["angle_step"]
            best = sorted(results, key=lambda x: x[0], reverse=True)[0]
            p.print(
                "Best power is",best[0], "at rpm",best[1], "using pitch",best[2], "°")
            p.print("That is a power gain of",best[0] - power_orig[0], "Watts.")
        best_chord_angles = chord_angles_orig + best[2]
        p.print("!!!!EOF!!!!")
        return best_chord_angles
