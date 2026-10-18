# Copyright (C) 2022 M. Smrekar

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import time
import traceback
from math import pi

import numpy as np
import scipy.interpolate
import scipy.optimize

from calculation import Calculator
from utils import Printer
from visualization import build_blade_lines


WHOLE_BLADE_MODE = 1
PENALTY_FACTOR = 1e5
ARRAY_SUM_VARIABLES = {"dQ", "dT", "dFn", "dFt", "dM"}


def optimize(inp_args, queue_pyqtgraph):
    """

    :param inp_args:
    :param queue_pyqtgraph:
    :return:
    """
    p = Printer(inp_args["return_print"])
    try:
        inp_args["v"] = inp_args["target_speed"]
        inp_args["rpm"] = inp_args["target_rpm"]
        inp_args["omega"] = 2 * pi * inp_args["rpm"] / 60

        input_variables = inp_args["optimization_inputs"]
        output_variables = inp_args["optimization_outputs"]
        target_variables = inp_args["optimization_targets"]

        p.print("Input variables:", input_variables)
        p.print("Output variables:", output_variables)
        p.print("Target_variables:", target_variables)

        if inp_args["turbine_type"] == 0:
            p.print("Turbine type: Wind turbine")
        elif inp_args["turbine_type"] == 1:
            p.print("Turbine type: Propeller")

        C = Calculator(inp_args)
        mode = int(inp_args.get("optimization_mode", 0) or 0)

        if mode == WHOLE_BLADE_MODE:
            optimize_whole_blade(inp_args, queue_pyqtgraph, p, C)
        else:
            optimize_section_by_section(inp_args, queue_pyqtgraph, p, C)

    except Exception as e:
        var = traceback.format_exc()
        p.print("Error in running optimizer: %s \n %s" % (str(e), var))
        inp_args["EOF"].value = True
        raise


def optimize_section_by_section(inp_args, queue_pyqtgraph, p, C):
    output_list = []
    input_variables = inp_args["optimization_inputs"]
    output_variables = inp_args["optimization_outputs"]
    target_variables = inp_args["optimization_targets"]
    optimized_theta = np.asarray(inp_args["theta"], dtype=float).copy()
    optimized_c = np.asarray(inp_args["c"], dtype=float).copy()
    preview_frame = 0

    if any(var[0] == "pitch" for var in input_variables):
        raise ValueError("Input variable 'pitch' is only supported in whole-blade optimization mode.")

    invalid_section_outputs = {"thrust", "power", "cp", "ct", "cq", "eff", "M", "Q", "blade_stall_percentage"}
    selected_section_outputs = {var for var, _coeff in output_variables}
    selected_section_targets = {var for var, _target, _coeff in target_variables}
    invalid_selected = (selected_section_outputs | selected_section_targets) & invalid_section_outputs
    if invalid_selected:
        raise ValueError("Whole-rotor variables %s are only supported in whole-blade optimization mode." %
                         sorted(invalid_selected))

    for n in range(len(inp_args["r"])):
        p.print("  Section_number is", n)

        list_queue_internal_x = []
        list_queue_internal_y = []

        section_inp_args = {
            "_r": inp_args["r"][n],
            "_c": inp_args["c"][n],
            "_theta": inp_args["theta"][n],
            "_dr": inp_args["dr"][n],
            "transition": C.transition_array[n],
            "_airfoil": C.airfoils_list[n],
            "_airfoil_prev": C.transition_foils[n][0],
            "_airfoil_next": C.transition_foils[n][1],
            "transition_coefficient": C.transition_foils[n][2],
            "max_thickness": C.max_thickness_array[n]
        }

        inputs_list = [i[0] for i in input_variables]
        bounds_list = [(b[1], b[2]) for b in input_variables]

        def fobj(input_numbers):
            nonlocal preview_frame
            for i in range(len(input_numbers)):
                section_inp_args[inputs_list[i]] = input_numbers[i]

            args = {**inp_args, **section_inp_args}

            d = C.calculate_section(**args, printer=p)
            if d is None or d is False:
                return 1e10

            fitness = 0

            for var, coeff in output_variables:
                value = d[var] * coeff
                fitness -= value
            for var, target_value, coeff in target_variables:
                comparison = abs(d[var] - target_value) * coeff
                fitness += comparison

            theta_preview = optimized_theta.copy()
            chord_preview = optimized_c.copy()
            theta_preview[n] = float(args["_theta"])
            chord_preview[n] = max(float(args["_c"]), 1e-4)
            preview_args = {**inp_args, "theta": list(theta_preview), "c": list(chord_preview)}
            preview_payload = {
                "frame": preview_frame,
                "data": build_blade_lines(preview_args, propeller_geom=inp_args["turbine_type"] == 1)["data"],
            }
            preview_frame += 1

            list_queue_internal_x.append(args["_theta"])
            list_queue_internal_y.append(fitness)
            queue_pyqtgraph[0] = [list_queue_internal_x, list_queue_internal_y, args["_theta"], fitness, preview_payload]

            return fitness

        decreasing = {"_theta"}

        if n > 0:
            for _vname in decreasing:
                if _vname in inputs_list:
                    index = inputs_list.index(_vname)
                    bounds_list[index] = (bounds_list[index][0], output_list[n - 1][index])

        p.print("Bounds:", bounds_list)

        initial_guess = [top for bottom, top in bounds_list]
        it = list(scipy.optimize.minimize(fobj, initial_guess, method="powell", bounds=bounds_list).x)

        p.print("best combination", it)
        output_list.append(it)

        for i, var_name in enumerate(inputs_list):
            if var_name == "_theta":
                optimized_theta[n] = it[i]
            elif var_name == "_c":
                optimized_c[n] = max(it[i], 1e-4)

    p.print("Final output:")
    p.print([v[0] for v in output_variables])
    for i in output_list:
        if len(i) > 1:
            p.print(i)
        else:
            p.print(i[0])

    final_args = {**inp_args, "theta": list(optimized_theta), "c": list(optimized_c)}
    final_result = C.run_array(**final_args)
    if final_result is not None and final_result is not False:
        attach_geometry_to_result(final_result, final_args)
        inp_args["return_results"].append(final_result)

    p.print("Done!")
    time.sleep(0.5)
    inp_args["EOF"].value = True


def optimize_whole_blade(inp_args, queue_pyqtgraph, p, C):
    baseline_theta = np.asarray(inp_args["theta"], dtype=float)
    baseline_c = np.asarray(inp_args["c"], dtype=float)
    radius = np.asarray(inp_args["r"], dtype=float)

    config = build_whole_blade_config(inp_args, baseline_theta, baseline_c, radius)
    if len(config["bounds"]) == 0:
        raise ValueError("Whole-blade optimization requires at least one active design variable.")

    p.print("Whole-blade optimization")
    p.print("Objective:", config["objective_label"])
    p.print("Solver:", config["solver"])
    p.print("Design variables:", [entry["label"] for entry in config["design_variables"]])

    history_x = []
    history_y = []
    best = {"fitness": np.inf, "result": None, "theta": baseline_theta.copy(), "c": baseline_c.copy(),
            "pitch": inp_args["pitch"], "iteration": 0}

    def fobj(x):
        theta, chord, pitch = decode_design_vector(x, config)
        args = {**inp_args, "theta": list(theta), "c": list(chord), "pitch": pitch}
        result = C.run_array(**args)

        if result is None or result is False:
            fitness = 1e12
        else:
            fitness = evaluate_whole_blade_fitness(result, theta, chord, pitch, config)

        history_x.append(len(history_x))
        history_y.append(fitness)

        if fitness < best["fitness"]:
            best["fitness"] = fitness
            best["result"] = result
            best["theta"] = np.asarray(theta, dtype=float)
            best["c"] = np.asarray(chord, dtype=float)
            best["pitch"] = pitch
            best["iteration"] = len(history_x) - 1
            preview_args = {**inp_args, "theta": list(best["theta"]), "c": list(best["c"]), "pitch": best["pitch"]}
            best["preview"] = {
                "frame": best["iteration"],
                "data": build_blade_lines(preview_args, propeller_geom=inp_args["turbine_type"] == 1)["data"],
            }

        queue_pyqtgraph[0] = [history_x, history_y, best["iteration"], best["fitness"], best.get("preview")]
        return fitness

    solver = config["solver"]
    initial_guess = np.asarray(config["initial_guess"], dtype=float)
    bounds = config["bounds"]

    if solver == "slsqp":
        result = scipy.optimize.minimize(
            fobj,
            initial_guess,
            method="SLSQP",
            bounds=bounds,
            options={"maxiter": config["num_iter"], "disp": False, "ftol": 1e-6},
        )
    else:
        mutation_upper = min(1.99, max(0.5, config["mut_coeff"]))
        result = scipy.optimize.differential_evolution(
            fobj,
            bounds=bounds,
            maxiter=config["num_iter"],
            popsize=config["population_size"],
            mutation=(0.5, mutation_upper),
            polish=True,
            disp=False,
        )

    final_theta, final_c, final_pitch = decode_design_vector(result.x, config)
    final_args = {**inp_args, "theta": list(final_theta), "c": list(final_c), "pitch": final_pitch}
    final_result = C.run_array(**final_args)
    final_fitness = evaluate_whole_blade_fitness(final_result, final_theta, final_c, final_pitch, config)

    if final_result is not None and final_result is not False and final_fitness <= best["fitness"]:
        best["fitness"] = final_fitness
        best["result"] = final_result
        best["theta"] = np.asarray(final_theta, dtype=float)
        best["c"] = np.asarray(final_c, dtype=float)
        best["pitch"] = final_pitch
        best["iteration"] = len(history_x)
        history_x.append(best["iteration"])
        history_y.append(best["fitness"])

    queue_pyqtgraph[0] = [history_x, history_y, best["iteration"], best["fitness"], best.get("preview")]

    p.print("Optimization success:", bool(getattr(result, "success", False)))
    p.print("Optimizer message:", getattr(result, "message", ""))
    p.print("Best fitness:", best["fitness"])
    p.print("Best pitch:", best["pitch"])
    p.print("Best chord:", [round(v, 5) for v in best["c"]])
    p.print("Best twist:", [round(v, 5) for v in best["theta"]])

    if best["result"] is not None:
        attach_geometry_to_result(best["result"], {**inp_args, "theta": list(best["theta"]), "c": list(best["c"]),
                                                    "pitch": best["pitch"]})
        p.print("Best thrust:", best["result"].get("thrust"))
        p.print("Best torque Q:", best["result"].get("Q"))
        p.print("Best power:", best["result"].get("power"))
        p.print("Best cp:", best["result"].get("cp"))
        p.print("Best eff:", best["result"].get("eff"))
        p.print("Best stall percentage:", best["result"].get("blade_stall_percentage"))
        inp_args["return_results"].append(best["result"])

    p.print("Done!")
    time.sleep(0.5)
    inp_args["EOF"].value = True


def build_whole_blade_config(inp_args, baseline_theta, baseline_c, radius):
    num_sections = len(radius)
    input_variables = inp_args["optimization_inputs"]
    objective_mode = int(inp_args.get("whole_blade_objective", 0) or 0)
    solver_mode = int(inp_args.get("whole_blade_solver", 0) or 0)

    twist_cp_count = max(2, int(inp_args.get("whole_blade_twist_cp", 5) or 5))
    chord_cp_count = max(2, int(inp_args.get("whole_blade_chord_cp", 4) or 4))
    root_locked = max(0, int(inp_args.get("whole_blade_root_cp_lock", 1) or 0))
    tip_locked = max(0, int(inp_args.get("whole_blade_tip_cp_lock", 1) or 0))

    config = {
        "objective_mode": objective_mode,
        "objective_label": get_objective_label(objective_mode),
        "solver": "slsqp" if solver_mode == 0 else "de",
        "num_iter": max(1, int(inp_args.get("num_iter", 30) or 30)),
        "population_size": max(1, int(inp_args.get("population", 10) or 10)),
        "mut_coeff": float(inp_args.get("mut_coeff", 0.8) or 0.8),
        "target_thrust": inp_args.get("whole_blade_target_thrust"),
        "max_power": inp_args.get("whole_blade_max_power"),
        "max_stall_percentage": inp_args.get("whole_blade_max_stall_percentage"),
        "monotonic_chord": bool(inp_args.get("whole_blade_monotonic_chord", False)),
        "output_variables": inp_args.get("optimization_outputs", []),
        "target_variables": inp_args.get("optimization_targets", []),
        "radius": np.asarray(radius, dtype=float),
        "baseline_theta": np.asarray(baseline_theta, dtype=float),
        "baseline_c": np.asarray(baseline_c, dtype=float),
        "baseline_pitch": float(inp_args.get("pitch", 0.0) or 0.0),
        "design_variables": [],
        "bounds": [],
        "initial_guess": [],
        "num_sections": num_sections,
    }

    normalized_radius = normalize_radius(radius)

    for name, lower, upper in input_variables:
        if name == "_theta":
            add_profile_design_variables(
                config=config,
                var_name=name,
                label_prefix="theta_cp",
                baseline_profile=baseline_theta,
                normalized_radius=normalized_radius,
                cp_count=twist_cp_count,
                lower=lower,
                upper=upper,
                root_locked=root_locked,
                tip_locked=tip_locked,
            )
        elif name == "_c":
            add_profile_design_variables(
                config=config,
                var_name=name,
                label_prefix="chord_cp",
                baseline_profile=baseline_c,
                normalized_radius=normalized_radius,
                cp_count=chord_cp_count,
                lower=max(lower, 1e-4),
                upper=max(upper, 1e-4),
                root_locked=root_locked,
                tip_locked=tip_locked,
            )
        elif name == "pitch":
            config["design_variables"].append({
                "type": "scalar",
                "var_name": "pitch",
                "index": 0,
                "label": "pitch",
            })
            config["bounds"].append((float(lower), float(upper)))
            config["initial_guess"].append(config["baseline_pitch"])

    return config


def add_profile_design_variables(config, var_name, label_prefix, baseline_profile, normalized_radius, cp_count,
                                 lower, upper, root_locked, tip_locked):
    cp_positions = np.linspace(0.0, 1.0, cp_count)
    baseline_values = np.interp(cp_positions, normalized_radius, baseline_profile)
    max_root_locked = min(root_locked, cp_count)
    max_tip_locked = min(tip_locked, cp_count - max_root_locked)

    config.setdefault("profile_definitions", {})[var_name] = {
        "cp_positions": cp_positions,
        "baseline_values": baseline_values,
        "normalized_radius": normalized_radius,
        "root_locked": max_root_locked,
        "tip_locked": max_tip_locked,
    }

    for i in range(cp_count):
        if i < max_root_locked or i >= cp_count - max_tip_locked:
            continue
        config["design_variables"].append({
            "type": "profile_cp",
            "var_name": var_name,
            "index": i,
            "label": "%s_%d" % (label_prefix, i),
        })
        config["bounds"].append((float(lower), float(upper)))
        value = float(np.clip(baseline_values[i], lower, upper))
        config["initial_guess"].append(value)


def normalize_radius(radius):
    radius = np.asarray(radius, dtype=float)
    if len(radius) == 1:
        return np.array([0.0])
    span = radius[-1] - radius[0]
    if span <= 0:
        return np.linspace(0.0, 1.0, len(radius))
    return (radius - radius[0]) / span


def decode_design_vector(x, config):
    profile_values = {}
    for name, definition in config.get("profile_definitions", {}).items():
        profile_values[name] = np.array(definition["baseline_values"], dtype=float)

    pitch = config["baseline_pitch"]

    for dv, value in zip(config["design_variables"], x):
        if dv["type"] == "profile_cp":
            profile_values[dv["var_name"]][dv["index"]] = value
        elif dv["type"] == "scalar" and dv["var_name"] == "pitch":
            pitch = float(value)

    theta = config["baseline_theta"].copy()
    chord = config["baseline_c"].copy()

    if "_theta" in profile_values:
        theta = evaluate_profile(profile_values["_theta"], config["profile_definitions"]["_theta"])
    if "_c" in profile_values:
        chord = evaluate_profile(profile_values["_c"], config["profile_definitions"]["_c"])

    chord = np.clip(chord, 1e-4, None)
    return theta, chord, pitch


def evaluate_profile(control_point_values, definition):
    cp_positions = definition["cp_positions"]
    normalized_radius = definition["normalized_radius"]

    if len(control_point_values) >= 3:
        interpolator = scipy.interpolate.PchipInterpolator(cp_positions, control_point_values)
        values = interpolator(normalized_radius)
    else:
        values = np.interp(normalized_radius, cp_positions, control_point_values)

    return np.asarray(values, dtype=float)


def evaluate_whole_blade_fitness(result, theta, chord, pitch, config):
    if result is None or result is False:
        return 1e12

    objective_mode = config["objective_mode"]
    fitness = 0.0

    if objective_mode == 0:
        target_thrust = require_value(config["target_thrust"], "Whole-blade target thrust")
        fitness += safe_abs(get_result_value(result, "Q"))
        fitness += normalized_penalty(get_result_value(result, "thrust") - target_thrust, target_thrust)
    elif objective_mode == 1:
        max_power = require_value(config["max_power"], "Whole-blade max power")
        fitness -= safe_abs(get_result_value(result, "thrust"))
        power_over = max(0.0, get_result_value(result, "power") - max_power)
        fitness += normalized_penalty(power_over, max_power)
    elif objective_mode == 2:
        fitness -= get_result_value(result, "cp")
    elif objective_mode == 3:
        target_thrust = require_value(config["target_thrust"], "Whole-blade target thrust")
        fitness -= get_result_value(result, "eff")
        fitness += normalized_penalty(get_result_value(result, "thrust") - target_thrust, target_thrust)
    else:
        for var, coeff in config["output_variables"]:
            fitness -= coeff * get_result_value(result, var)
        for var, target_value, coeff in config["target_variables"]:
            fitness += coeff * safe_abs(get_result_value(result, var) - target_value)

    max_stall_percentage = config.get("max_stall_percentage")
    if max_stall_percentage is not None:
        stall_over = max(0.0, get_result_value(result, "blade_stall_percentage") - max_stall_percentage)
        fitness += normalized_penalty(stall_over, max(max_stall_percentage, 1e-6))

    max_power = config.get("max_power")
    if max_power is not None and objective_mode != 1:
        power_over = max(0.0, get_result_value(result, "power") - max_power)
        fitness += normalized_penalty(power_over, max_power)

    if config.get("monotonic_chord"):
        idx = int(np.argmax(chord))
        if idx < len(chord) - 1:
            increases = np.maximum(np.diff(chord[idx:]), 0.0)
            if len(increases) > 0:
                fitness += PENALTY_FACTOR * float(np.sum(increases / max(np.max(chord), 1e-6)))

    return float(fitness)


def get_result_value(result, var):
    value = result[var]
    if isinstance(value, (list, tuple)):
        value = np.asarray(value)
    if isinstance(value, np.ndarray):
        if var in ARRAY_SUM_VARIABLES:
            return float(np.sum(value))
        return float(np.mean(value))
    return float(value)


def normalized_penalty(error_value, scale):
    return PENALTY_FACTOR * safe_abs(error_value) / max(safe_abs(scale), 1.0)


def safe_abs(value):
    return float(abs(value))


def require_value(value, label):
    if value is None:
        raise ValueError("%s is required for the selected whole-blade objective." % label)
    return float(value)


def get_objective_label(index):
    labels = {
        0: "Min torque at target thrust",
        1: "Max thrust at power limit",
        2: "Max Cp",
        3: "Max efficiency at target thrust",
        4: "Weighted custom objective",
    }
    return labels.get(index, "Weighted custom objective")


def attach_geometry_to_result(result, inp_args):
    result["optimized_geometry"] = {
        "r_in": list(inp_args.get("r", [])),
        "c_in": list(inp_args.get("c", [])),
        "theta_in": list(inp_args.get("theta", [])),
        "foils_in": list(inp_args.get("foils", inp_args.get("foils_in", []))),
        "pitch": inp_args.get("pitch"),
    }
