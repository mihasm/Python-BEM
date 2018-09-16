__author__ = "Miha Smrekar"
__credits__ = ["Miha Smrekar"]
__license__ = "GPL"
__version__ = "0.2.3"
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
    return_print=None,
    print_all=None,
    print_out=None
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
        print_all=print_all,
        print_out=print_out,
        return_print=return_print
    )
    return results


class Optimizer:
    """
    Class used for optimizing blade geometry.
    TODO: add input arguments to init, so that user will have to pass geometry data and functions for initialisation.
    """

    def __init__(self, r, c, theta, dr, B, R, Rhub, f_c_L, f_c_D, return_print, print_out, print_all,*args,**kwargs):
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

        self.return_print = return_print
        self.print_out = print_out
        self.print_all = print_all

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
                return_print=self.return_print,
                print_all = self.print_all,
                print_out = self.print_out
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
                return_print=self.return_print,
                print_all = self.print_all,
                print_out = self.print_out
            )

        if power_only:
            return res["power"]
        else:
            return res

    def direction(self, _r, _delta, wind_speed, rpm, *args,**kwargs):
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

    def optimize_angles(self, target_speed, target_rpm, delta_start, decrease_factor, min_delta, return_print, return_results,*args,**kwargs):
        """
        Optimizes angles of attack by changing twist angles.
        Prints results.
        """
        power_preliminary = self._power(wind_speed=target_speed, rpm=target_rpm)

        return_print.append(
            "Optimizing angles for target speed "+\
            str(target_speed)+\
            " [m/s], target rpm "+\
            str(target_rpm)+\
            " [RPM], number of blades:"+\
            str(5)+\
            ", tip radius:"+\
            str(self.R)+" [m]\n"
            ", hub radius:"+\
            str(self.Rhub)+" [m]\n"
        )

        #min_delta = 0.01

        return_print.append("power at start is "+str(power_preliminary)+"\n")

        for r in range(len(self.sections_radius)):
            return_print.append("Currently optimizing section "+str(r)+"\n")
            power = self._power(wind_speed=target_speed, rpm=target_rpm)
            delta = delta_start
            current_delta = 0.0

            while True:
                return_print.append("    Current delta "+ str(delta)+"\n")
                if delta < min_delta:
                    return_print.append("    Delta limit reached\n")
                    break

                d = self.direction(r, delta, wind_speed=target_speed, rpm=target_rpm)
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
                            "        Changing blade section "+\
                            str(r)+\
                            " angle from "+\
                            str(self.chord_angles[r])+\
                            " to "+\
                            str(self.chord_angles[r]) + str(current_delta) + "\n"
                        )
                        self.chord_angles[r] += current_delta
                        power = self._power(wind_speed=target_speed, rpm=target_rpm)
                        if power <= power_old:
                            return_print.append(
                                "        New power "+\
                                str(power)+\
                                " is less than or equal to old power "+\
                                str(power_old)+"\n"
                            )
                            break
                        return_print.append(
                            "        New power "+\
                            str(power)+\
                            " is greater than old power "+\
                            str(power_old)+ "\n"
                        )
                        power_old = power
                    return_print.append("    Changing delta from "+str(delta)+" to "+str(delta * decrease_factor)+"\n")

                delta = delta * decrease_factor

        power_new = self._power(wind_speed=target_speed, rpm=target_rpm)

        return_print.append("power at start was "+str( power_preliminary)+"\n")
        return_print.append("power at end is "+str( power_new)+"\n")
        return_print.append("Power increase "+str( power_new - power_preliminary)+"\n")
        percentage_increase = power_new / power_preliminary * 100
        return_print.append("Percentage "+str( percentage_increase)+"\n")
        old_angles = self.chord_angles_orig
        new_angles = self.chord_angles
        return_print.append("Old angles "+str( old_angles)+"\n")
        return_print.append("New angles "+str( new_angles)+"\n")
        return_print.append("----------------------"+"\n")

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
        if target_rpm != None:
            return_print.append(
                "\nOptimizing for wind speed of "+\
                str(target_speed)+\
                "m/s... and rpm "+\
                str(target_rpm)+" RPM\n"
            )

            power_orig = self._power(wind_speed=target_speed, rpm=target_rpm)
            return_print.append(
                "Without turning the blade, the power is: "+\
                "%.2f" % round(power_orig)+\
                " [W] at wind speed "+\
                str(target_speed)+\
                "[m/s] and rotational velocity "+\
                str(target_rpm)+\
                "[RPM]"+"\n"
            )

            current_add_angle = min_add_angle
            results = []
            while current_add_angle <= max_add_angle:
                return_print.append("Testing pitch: "+str(current_add_angle)+"° at rpm "+ str(target_rpm)+" [RPM]\n")
                power = self._power(
                    add_angle=current_add_angle, wind_speed=target_speed, rpm=target_rpm
                )
                results.append((current_add_angle, power))
                current_add_angle += step
                return_print.append("    power:"+str(power)+"\n")
                return_print.append("---\n")

            best_angle, best_power = sorted(results, key=lambda x: x[1], reverse=True)[
                0
            ]

            return_print.append(
                "The best angle to turn would be:"+\
                "%.3f" % round(best_angle, 2)+\
                "°, at this pitch angle the power is:"+\
                "%.2f" % round(best_power)+\
                "Watts"+"\n"
            )
            return_print.append(
                "That is a power gain of "+\
                "%.2f" % round(best_power - power_orig, 2)+\
                "Watts"+"\n"
            )
            return_print.append("!!!!EOF!!!!")
            return best_angle
        else:
            return_print.append("Getting best power with original pitch.\n")
            power_orig = (None, None)  # power, rpm
            for rpm in numpy.linspace(start=100, stop=3000, num=30):
                res = self._power(wind_speed=target_speed, rpm=rpm, power_only=False)
                if 0.0 < res["cp"] <= 0.6:
                    p = res["power"]
                    return_print.append("    rpm:"+ str(rpm)+ " p:"+ str(p)+"\n")
                    if power_orig[0] == None and power_orig[1] == None:
                        power_orig = (p, rpm)
                    elif p > power_orig[0]:
                        power_orig = (p, rpm)
            return_print.append("Original power"+str(power_orig)+"\n")
            current_add_angle = min_add_angle
            results = []
            return_print.append("Getting power for different pitches"+"\n")
            while current_add_angle <= max_add_angle:
                return_print.append("    Current add angle:"+str(current_add_angle)+"\n")
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
                        return_print.append("        rpm:"+ str(rpm)+ "p:"+ str(p)+"\n")
                        if (
                            best_result_angle[0] == None
                            and best_result_angle[1] == None
                        ):
                            best_result_angle = (p, rpm, current_add_angle)
                        if p > best_result_angle[0]:
                            best_result_angle = (p, rpm, current_add_angle)
                    else:
                        p = res["power"]
                        return_print.append("        rpm:"+ str(rpm)+ " p:"+ str(p)+ " not possible / bad result"+"\n")
                results.append(best_result_angle)
                current_add_angle += step
            best = sorted(results, key=lambda x: x[0], reverse=True)[0]
            return_print.append("Best power is "+ str(best[0])+ " at rpm "+ str(best[1])+ " using pitch "+ str(best[2])+"°\n")
            return_print.append("That is a power gain of "+ str(best[0] - power_orig[0])+ " Watts."+"\n")
        best_chord_angles = self.chord_angles_orig + best[2]
        return_print.append("!!!!EOF!!!!")
        return best_chord_angles
