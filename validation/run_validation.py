# Python BEM - Blade Element Momentum Theory Software.

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

"""
Runs the BEM solver on every validation case and compares it with the measurements.

    python3 validation/run_validation.py                      # all cases, default solver variants
    python3 validation/run_validation.py --cases ntnu_bt1 mexico --variants buhl_prandtl buhl_prandtl_snel
    python3 validation/run_validation.py --list               # cases and variants

Writes validation/results/<case>/<experiment>.png and .csv, and validation/results/summary.md / summary.csv.
"""

import argparse
import math
import os
import sys
import time
import traceback
import multiprocessing
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import bem_validation as bv  # noqa: E402

RESULTS_DIR = HERE / "results"



def _arr(res, key):
    return np.asarray(res[key], float)


def _cn_chord(res):
    al = np.radians(_arr(res, "alpha"))
    return _arr(res, "Cl") * np.cos(al) + _arr(res, "Cd") * np.sin(al)


def _ct_chord(res):
    al = np.radians(_arr(res, "alpha"))
    return _arr(res, "Cl") * np.sin(al) - _arr(res, "Cd") * np.cos(al)


def _q_c(res, rho):
    return 0.5 * rho * _arr(res, "Vrel_norm") ** 2 * _arr(res, "c")


# spanwise column (fixed meaning, see validation/README.md) -> model value from the solver results
SPANWISE = {
    "fn_n_per_m": lambda res, rho: _arr(res, "dFn/n"),
    "ft_n_per_m": lambda res, rho: _arr(res, "dFt/n"),
    "cn_rotor": lambda res, rho: _arr(res, "Cn"),
    "ct_rotor": lambda res, rho: _arr(res, "Ct"),
    "cn_chord": lambda res, rho: _cn_chord(res),
    "ct_chord": lambda res, rho: _ct_chord(res),
    "fn_chord_n_per_m": lambda res, rho: _cn_chord(res) * _q_c(res, rho),
    "ft_chord_n_per_m": lambda res, rho: _ct_chord(res) * _q_c(res, rho),
    "cl": lambda res, rho: _arr(res, "Cl"),
    "cd": lambda res, rho: _arr(res, "Cd"),
    "alpha_deg": lambda res, rho: _arr(res, "alpha"),
    "phi_deg": lambda res, rho: _arr(res, "phi"),
    "a": lambda res, rho: _arr(res, "a"),
    "a_prime": lambda res, rho: _arr(res, "a'"),
    "re": lambda res, rho: _arr(res, "Re"),
}
SCALAR_COLUMNS = ["cp", "ct", "cq", "eta", "power_w", "thrust_n", "torque_nm",
                  "root_flap_moment_nm", "root_edge_moment_nm"]
LABELS = {
    "cp": "$C_P$", "ct": "$C_T$", "cq": "$C_Q$", "eta": r"$\eta$", "power_w": "P [W]", "thrust_n": "T [N]",
    "torque_nm": "Q [N m]", "root_flap_moment_nm": "flap moment [N m]", "root_edge_moment_nm": "edge moment [N m]",
    "tsr": r"$\lambda$", "J": "J", "v_ms": "v [m/s]", "r_m": "r [m]", "r_R": "r/R",
    "fn_n_per_m": "$F_n$ [N/m]", "ft_n_per_m": "$F_t$ [N/m]", "cn_rotor": "$C_n$ (rotor plane)",
    "ct_rotor": "$C_t$ (rotor plane)", "cn_chord": "$C_n$ (chord)", "ct_chord": "$C_t$ (chord)",
    "fn_chord_n_per_m": "$F_n$ chord [N/m]", "ft_chord_n_per_m": "$F_t$ chord [N/m]",
    "cl": "$C_l$", "alpha_deg": r"$\alpha$ [°]", "a": "a",
}
GROUP_KEYS = ["pitch_deg", "yaw_deg", "rpm", "v_ms"]  # conditions that make separate curves


def x_column(case, exp):
    cols = exp.data.columns
    if exp.kind == "spanwise":
        return "r_m" if "r_m" in cols else "r_R"
    for c in (exp.meta.get("x"), "tsr", "J", "v_ms"):
        if c and c in cols:
            return c
    raise ValueError("%s: no x column (tsr / J / v_ms)" % exp.file)


def y_columns(exp, x):
    wanted = exp.meta.get("y")
    if isinstance(wanted, str):
        wanted = [wanted]
    known = SCALAR_COLUMNS if exp.kind != "spanwise" else list(SPANWISE)
    cols = [c for c in (wanted or exp.data.columns) if c in exp.data.columns and c in known and c != x]
    return cols


def group_columns(exp, x):
    """Columns that define separate curves (e.g. pitch settings, rpm), excluding the x axis."""
    out = []
    for k in GROUP_KEYS:
        if k in exp.data.columns and k != x and exp.data[k].nunique() > 1:
            # v_ms / rpm also vary along a sweep (a tsr sweep at constant v, a J schedule with rpm = f(J)):
            # they make separate curves only when each value holds a curve of its own (3 points or more)
            if k in ("v_ms", "rpm") and exp.data.groupby(out + [k]).size().median() < 3:
                continue
            out.append(k)
    return out


def sweep_rows(group_df, x, n=None):
    """Dense rows along x for the model curve, holding the group's other conditions."""
    n = n or SWEEP_POINTS
    lo, hi = group_df[x].min(), group_df[x].max()
    span = hi - lo
    xs = np.linspace(max(lo - 0.05 * span, 1e-3 if x != "v_ms" else 0.5), hi + 0.05 * span, n)
    base = group_df.iloc[0]
    rows = []
    for xv in xs:
        row = base.copy()
        row[x] = xv
        # along a tsr/J sweep usually one of v and rpm is held (the one constant in the data) and the other
        # follows from x; when both vary (an rpm schedule over J), rpm is interpolated along x
        varying = [k for k in ("v_ms", "rpm") if k != x and k in group_df.columns and group_df[k].nunique() > 1]
        for k in varying:
            row[k] = np.nan
        if len(varying) == 2:
            d = group_df.sort_values(x)
            row["rpm"] = float(np.interp(xv, d[x], d["rpm"]))
        rows.append(row)
    return pd.DataFrame(rows)


def gauge_radius(case, exp, row):
    """Radius the measured root moments refer to: column or defaults `gauge_r_m`, else the rotor axis."""
    if "gauge_r_m" in row and pd.notna(row["gauge_r_m"]):
        return float(row["gauge_r_m"])
    return float({**case.meta.get("defaults", {}), **exp.defaults}.get("gauge_r_m", 0.0))


def run_scalar(runner, case, exp, df):
    out = []
    loads = any(c.startswith("root_") for c in df.columns)
    for _, row in df.iterrows():
        v, rpm, pitch, yaw, rho, nu = runner.conditions(exp, row)
        res = runner.run_point(v, rpm, pitch, yaw, rho, nu)
        s = bv.summarize(res, case, gauge_radius(case, exp, row) if loads else None)
        s.update({"v_ms_model": v, "rpm_model": rpm, "pitch_model": pitch, "yaw_model": yaw})
        out.append(s)
    return pd.DataFrame(out, index=df.index)


def run_experiment(runner, case, exp):
    """Returns (points, curves): model values at the measured rows, and dense model curves per group."""
    x = x_column(case, exp)
    df = exp.data
    if exp.kind == "spanwise":
        return run_spanwise(runner, case, exp, x)
    points = pd.concat([df.reset_index(drop=True),
                        run_scalar(runner, case, exp, df.reset_index(drop=True)).add_prefix("model_")],
                       axis=1)
    curves = []
    gcols = group_columns(exp, x)
    groups = df.groupby(gcols) if gcols else [((), df)]
    for key, g in groups:
        dense = sweep_rows(g, x).reset_index(drop=True)
        m = run_scalar(runner, case, exp, dense).add_prefix("model_")
        c = pd.concat([dense[[x] + gcols].reset_index(drop=True), m], axis=1)
        curves.append(c)
    return points, pd.concat(curves, ignore_index=True)


def run_spanwise(runner, case, exp, x):
    df = exp.data.reset_index(drop=True)
    R = float(case.meta["R"])
    keys = [k for k in ("v_ms", "rpm", "tsr", "J", "pitch_deg", "yaw_deg") if k in df.columns]
    groups = df.groupby(keys) if keys else [((), df)]
    pts, curves = [], []
    for key, g in groups:
        row = g.iloc[0]
        v, rpm, pitch, yaw, rho, nu = runner.conditions(exp, row)
        res = runner.run_point(v, rpm, pitch, yaw, rho, nu)
        ok = res is not None and "error" not in res
        r_model = np.asarray(res["r"]) if ok else np.array([])
        r_model_x = r_model if x == "r_m" else r_model / R
        p = g.copy()
        c = pd.DataFrame({x: r_model_x})
        for k in keys:
            c[k] = row[k]
        for col, model in SPANWISE.items():
            if col in g.columns:
                if ok:
                    vals = model(res, rho)
                    p["model_" + col] = np.interp(g[x].to_numpy(float), r_model_x, vals)
                    c["model_" + col] = vals
                else:
                    p["model_" + col] = np.nan
        p["model_converged"] = ok
        pts.append(p)
        curves.append(c)
    return pd.concat(pts, ignore_index=True), pd.concat(curves, ignore_index=True)


def metrics(points, ycols, x):
    out = {}
    for y in ycols:
        m = "model_" + y
        if m not in points:
            continue
        d = points[[x, y, m]].dropna()
        if d.empty:
            continue
        err = d[m] - d[y]
        scale = max(d[y].abs().max(), 1e-12)
        out[y] = {"n": len(d), "rmse": float(np.sqrt(np.mean(err ** 2))), "bias": float(err.mean()),
                  "nrmse_pct": float(100 * np.sqrt(np.mean(err ** 2)) / scale),
                  "max_abs": float(err.abs().max())}
        if y in ("cp", "eta"):
            i_e, i_m = d[y].idxmax(), d[m].idxmax()
            out[y]["peak_exp"] = float(d.loc[i_e, y])
            out[y]["peak_model"] = float(d.loc[i_m, m])
            out[y]["peak_err_pct"] = float(100 * (d.loc[i_m, m] - d.loc[i_e, y]) / max(abs(d.loc[i_e, y]), 1e-12))
            out[y]["x_at_peak_exp"] = float(d.loc[i_e, x])
            out[y]["x_at_peak_model"] = float(d.loc[i_m, x])
    return out


SWEEP_POINTS = 30


def count_points(case, exp):
    """How many solver points run_experiment will compute for one experiment file."""
    x = x_column(case, exp)
    if exp.kind == "spanwise":
        keys = [k for k in ("v_ms", "rpm", "tsr", "J", "pitch_deg", "yaw_deg") if k in exp.data.columns]
        return exp.data.groupby(keys).ngroups if keys else 1
    g = group_columns(exp, x)
    return len(exp.data) + (exp.data.groupby(g).ngroups if g else 1) * SWEEP_POINTS


def track_progress(runner, queue):
    """Wraps runner.run_point so every computed point is reported to the main process (True = converged)."""
    run_point = runner.run_point

    def tracked(v, rpm, pitch=0.0, yaw=0.0, rho=None, nu=None):
        res = run_point(v, rpm, pitch, yaw, rho, nu)
        if queue is not None:
            queue.put(res is not None and "error" not in res)
        return res

    runner.run_point = tracked


class ProgressBar:
    """One line for the whole run: points done of all jobs, failed points, elapsed time and ETA."""

    def __init__(self, total, jobs):
        self.total, self.jobs = max(total, 1), jobs
        self.done = self.failed = self.jobs_done = 0
        self.t0 = self.last = time.time()
        self.tty = sys.stdout.isatty()

    def update(self, force=False):
        now = time.time()
        # on a terminal redraw in place; in a log file one line every 10 s
        if not force and now - self.last < (0.2 if self.tty else 10.0):
            return
        self.last = now
        frac = min(self.done / self.total, 1.0)
        elapsed = now - self.t0
        eta = elapsed / self.done * (self.total - self.done) if self.done else 0
        bar = "#" * int(30 * frac) + "." * (30 - int(30 * frac))
        line = "[%s] %5.1f%%  %d/%d points  %d failed  jobs %d/%d  %s elapsed  ETA %s" % (
            bar, 100 * frac, self.done, self.total, self.failed, self.jobs_done, self.jobs,
            _hms(elapsed), _hms(eta))
        print(("\r" + line) if self.tty else line, end="" if self.tty else "\n", flush=True)


def _hms(seconds):
    seconds = int(seconds)
    return "%d:%02d:%02d" % (seconds // 3600, seconds // 60 % 60, seconds % 60)


def is_static(exp):
    """Hover / static thrust data: zero free-stream speed in every row."""
    d, defaults = exp.data, exp.defaults
    if "J" in d:
        return bool((d["J"].fillna(0) == 0).all())
    if "v_ms" in d:
        return bool((d["v_ms"].fillna(0) == 0).all())
    return defaults.get("v_ms", None) == 0


def count_job_points(case_path):
    case = bv.load_case(case_path)
    n = 0
    for exp in case.experiments:
        if exp.kind not in ("structural", "polar") and exp.data is not None and not is_static(exp):
            try:
                n += count_points(case, exp)
            except Exception:  # noqa: BLE001 - the experiment reports its own error when it runs
                pass
    return n


def run_case_variant(case_path, variant, n_sections, queue=None):
    t0 = time.time()
    case = bv.load_case(case_path)
    runner = bv.CaseRunner(case, variant, n_sections)
    out = {"case": case.id, "variant": variant, "experiments": {}, "errors": [], "notes": []}
    todo = [e for e in case.experiments if e.kind not in ("structural", "polar") and e.data is not None]
    for e in [e for e in todo if is_static(e)]:
        todo.remove(e)
        out["notes"].append("%s: static (v = 0) data — the solver's induction formulation needs v > 0, not run"
                            % e.file)
    track_progress(runner, queue)
    for exp in todo:
        try:
            points, curves = run_experiment(runner, case, exp)
            x = x_column(case, exp)
            out["experiments"][exp.file] = {"points": points, "curves": curves, "x": x,
                                            "y": y_columns(exp, x),
                                            "metrics": metrics(points, y_columns(exp, x), x)}
        except Exception:  # noqa: BLE001 - one broken file must not stop the rest
            out["errors"].append("%s: %s" % (exp.file, traceback.format_exc(limit=3)))
    out["seconds"] = time.time() - t0
    return out


def load_saved(case, variant, cdir):
    """The result of run_case_variant rebuilt from the CSVs a previous run wrote (for --replot)."""
    out = {"case": case.id, "variant": variant, "experiments": {}, "errors": []}
    for exp in case.experiments:
        stem = Path(exp.file).stem
        fp, fc = cdir / (stem + ".csv"), cdir / (stem + "_curves.csv")
        if not (fp.exists() and fc.exists()):
            continue
        pts, crv = pd.read_csv(fp), pd.read_csv(fc)
        pts = pts[pts["variant"] == variant].drop(columns="variant").reset_index(drop=True)
        crv = crv[crv["variant"] == variant].drop(columns="variant").reset_index(drop=True)
        if pts.empty:
            continue
        x = x_column(case, exp)
        out["experiments"][exp.file] = {"points": pts, "curves": crv, "x": x, "y": y_columns(exp, x),
                                        "metrics": metrics(pts, y_columns(exp, x), x)}
    return out


# ---------------------------------------------------------------- plotting

def plot_experiment(case, exp, by_variant, out_png):
    any_res = next(iter(by_variant.values()))
    x, ycols = any_res["x"], any_res["y"]
    if not ycols:
        return
    data = exp.data
    gcols = group_columns(exp, x) if exp.kind != "spanwise" else \
        [k for k in ("v_ms", "rpm", "tsr", "J", "pitch_deg", "yaw_deg") if k in data.columns and data[k].nunique() > 1]
    groups = list(data.groupby(gcols)) if gcols else [((), data)]
    n = len(ycols)
    if gcols:
        # one row of panels per operating group (pitch setting, rpm, wind speed ...)
        fig, axes = plt.subplots(len(groups), n, figsize=(4.6 * n, 3.4 * len(groups)), squeeze=False)
    else:
        fig, axes = plt.subplots(1, n, figsize=(4.8 * n, 4.0), squeeze=False)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    styles = ["-", "--", ":", "-.", (0, (5, 1)), (0, (3, 1, 1, 1))]
    panels = bool(gcols)
    by_group = False  # measurements black, every solver variant its own colour
    ranges = {}
    for gi, (key, g) in enumerate(groups):
        key = key if isinstance(key, tuple) else (key,)
        glabel = ", ".join("%s=%g" % (k, v) for k, v in zip(gcols, key))
        exp_color = colors[gi % len(colors)] if by_group else "k"
        row_axes = axes[gi] if panels else axes[0]
        for yi, y in enumerate(ycols):
            ax = row_axes[yi]
            vals = [g[y].to_numpy(float)]
            err = g[y + "_err"] if y + "_err" in g else None
            ax.errorbar(g[x], g[y], yerr=err, fmt="o", ms=4, color=exp_color, mfc="none", capsize=2,
                        label=("exp. " + glabel).strip() if by_group else "measured")
            for vi, (variant, res) in enumerate(by_variant.items()):
                c = res["curves"]
                if c is None or "model_" + y not in c:
                    continue
                sel = np.ones(len(c), bool)
                for k, v in zip(gcols, key):
                    sel &= np.isclose(c[k].astype(float), float(v))
                cc = c[sel].sort_values(x)
                vals.append(cc["model_" + y].to_numpy(float))
                ax.plot(cc[x], cc["model_" + y], linestyle=styles[vi % len(styles)],
                        color=exp_color if by_group else colors[vi % len(colors)], lw=1.4,
                        label=None if by_group else variant)
            ax.set_xlabel(LABELS.get(x, x))
            ax.set_ylabel(LABELS.get(y, y))
            ax.grid(alpha=0.3)
            if panels:
                ax.set_title(glabel, fontsize=9)
            # axis range from the values, not from the error bars (one huge error bar flattens everything)
            v = np.concatenate(vals)
            v = v[np.isfinite(v)]
            if v.size:
                lo, hi = v.min(), v.max()
                pad = 0.08 * (hi - lo or abs(hi) or 1.0)
                ax.set_ylim(lo - pad, hi + pad)
    handles, labels = axes[0][0].get_legend_handles_labels()
    if by_group:
        from matplotlib.lines import Line2D
        for vi, variant in enumerate(by_variant):
            handles.append(Line2D([], [], color="0.3", linestyle=styles[vi % len(styles)], lw=1.4))
            labels.append(variant)
    src = exp.meta.get("source", "")
    fig.suptitle("%s — %s" % (case.meta.get("name", case.id), Path(exp.file).stem))
    fig.tight_layout()
    rows_legend = math.ceil(len(labels) / 4)
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=min(4, len(labels)),
               fontsize=8, frameon=False)
    fig.text(0.01, -0.035 * (rows_legend + 1), ("source: " + src)[:200], fontsize=7, color="0.4")
    fig.savefig(out_png, dpi=130, bbox_inches="tight")
    plt.close(fig)


def plot_geometry(case, out_png, n_sections):
    r, dr, c, theta, foils = bv.blade_sections(case, n_sections)
    g = case.geometry
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
    ax[0].plot(g["r_m"], g["chord_m"], "o", mfc="none", label="table")
    ax[0].plot(r, c, ".-", lw=0.8, label="solver sections")
    ax[1].plot(g["r_m"], g["twist_deg"], "o", mfc="none")
    ax[1].plot(r, theta, ".-", lw=0.8)
    for i, f in enumerate(foils):
        if i == 0 or foils[i - 1] != f:
            ax[0].annotate(f, (r[i], c[i]), fontsize=7, rotation=60)
    ax[0].set_ylabel("chord [m]")
    ax[1].set_ylabel("twist [°]")
    for a in ax:
        a.set_xlabel("r [m]")
        a.grid(alpha=0.3)
    ax[0].legend(fontsize=8)
    fig.suptitle("%s — geometry" % case.meta.get("name", case.id))
    fig.tight_layout()
    fig.savefig(out_png, dpi=110)
    plt.close(fig)


# ---------------------------------------------------------------- main

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cases", nargs="*", help="case folder names (default: all)")
    ap.add_argument("--variants", nargs="*", default=bv.DEFAULT_VARIANTS)
    ap.add_argument("--sections", type=int, default=12)
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--out", default=str(RESULTS_DIR))
    ap.add_argument("--cases-dir", default=str(bv.CASES_DIR))
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--replot", action="store_true",
                    help="redraw charts and summary from the CSVs of the last run in --out, without the solver")
    a = ap.parse_args(argv)

    aero = ("performance", "power_curve", "spanwise", "loads")
    case_paths = [p for p in bv.list_cases(a.cases, a.cases_dir)
                  if any(e.kind in aero for e in bv.load_case(p).experiments)]
    if a.list:
        for p in case_paths:
            c = bv.load_case(p)
            print("%-28s %-10s %s" % (p.name, "propeller" if c.is_propeller else "turbine",
                                      ", ".join("%s(%s)" % (e.name, e.kind) for e in c.experiments)))
        print("\nvariants:", ", ".join(bv.VARIANTS))
        return 0
    for v in a.variants:
        if v not in bv.VARIANTS:
            ap.error("unknown variant %s (choose from %s)" % (v, ", ".join(bv.VARIANTS)))

    out_dir = Path(a.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(bv.ROOT)  # the solver writes airfoil .dat files under ./foils
    case_variants = {p: bv.variants_for(bv.load_case(p), a.variants) for p in case_paths}
    jobs = [(p, v) for p in case_paths for v in case_variants[p]]
    results = {}
    if a.replot:
        for p in case_paths:
            for v in case_variants[p]:
                results[(p, v)] = load_saved(bv.load_case(p), v, out_dir / p.name)
        jobs = []
    if jobs:
        points = {p: count_job_points(p) for p in case_paths}
        bar = ProgressBar(sum(points[p] for p, v in jobs), len(jobs))
        print("%d cases x variants = %d jobs, %d solver points, %d processes" % (
            len(case_paths), len(jobs), bar.total, a.jobs), flush=True)
        manager = multiprocessing.Manager()
        queue = manager.Queue()
        with ProcessPoolExecutor(max_workers=a.jobs) as pool:
            futs = {pool.submit(run_case_variant, p, v, a.sections, queue): (p, v) for p, v in jobs}
            pending = set(futs)
            while pending:
                finished, pending = wait(pending, timeout=0.2, return_when=FIRST_COMPLETED)
                while not queue.empty():
                    ok = queue.get()
                    bar.done += 1
                    bar.failed += 0 if ok else 1
                for f in finished:
                    p, v = futs[f]
                    try:
                        r = f.result()
                    except Exception:  # noqa: BLE001
                        r = {"case": p.name, "variant": v, "experiments": {},
                             "errors": [traceback.format_exc(limit=3)]}
                    results[(p, v)] = r
                    bar.jobs_done += 1
                    if r["errors"]:
                        print("\n%s / %s: %d experiment file(s) failed, see summary.md" % (
                            p.name, v, len(r["errors"])), flush=True)
                bar.update()
        bar.update(force=True)
        print("", flush=True)
        manager.shutdown()

    rows, md = [], ["# BEM validation results", "",
                    "Generated by `validation/run_validation.py` — variants: %s, %d sections." %
                    (", ".join(a.variants), a.sections), "Solver: `%s`." % bv.ROOT, ""]
    for p in case_paths:
        case = bv.load_case(p)
        cdir = out_dir / p.name
        cdir.mkdir(exist_ok=True)
        plot_geometry(case, cdir / "geometry.png", a.sections)
        md += ["## %s (`%s`)" % (case.meta.get("name", case.id), p.name), "", "![geometry](%s/geometry.png)" % p.name, ""]
        if not case_variants[p]:
            md += ["No requested variant applies to this rotor type.", ""]
        for v in case_variants[p]:
            for e in results[(p, v)]["errors"]:
                md += ["**error** (%s): `%s`" % (v, e.strip().splitlines()[-1]), ""]
        if case_variants[p]:
            for n in results[(p, case_variants[p][0])].get("notes", []):
                md += ["*%s*" % n, ""]
        for exp in case.experiments:
            by_variant = {v: results[(p, v)]["experiments"][exp.file] for v in case_variants[p]
                          if exp.file in results[(p, v)]["experiments"]}
            if not by_variant:
                continue
            stem = Path(exp.file).stem
            plot_experiment(case, exp, by_variant, cdir / (stem + ".png"))
            pts, crv = [], []
            for v, res in by_variant.items():
                q = res["points"].copy()
                q.insert(0, "variant", v)
                pts.append(q)
                q = res["curves"].copy()
                q.insert(0, "variant", v)
                crv.append(q)
                for y, m in res["metrics"].items():
                    rows.append({"case": p.name, "experiment": stem, "variant": v, "quantity": y, **m})
            pd.concat(pts).to_csv(cdir / (stem + ".csv"), index=False)
            pd.concat(crv).to_csv(cdir / (stem + "_curves.csv"), index=False)
            md += ["### %s — %s" % (stem, exp.kind), "", "Source: %s" % exp.meta.get("source", "?"), "",
                   "![%s](%s/%s.png)" % (stem, p.name, stem), ""]
            tab = [r for r in rows if r["case"] == p.name and r["experiment"] == stem]
            if tab:
                md += ["| variant | quantity | n | RMSE | bias | NRMSE % | peak exp | peak model | peak err % |",
                       "|---|---|---|---|---|---|---|---|---|"]
                for r in tab:
                    md.append("| %s | %s | %d | %.4g | %+.4g | %.1f | %s | %s | %s |" % (
                        r["variant"], r["quantity"], r["n"], r["rmse"], r["bias"], r["nrmse_pct"],
                        "%.3f" % r["peak_exp"] if "peak_exp" in r else "",
                        "%.3f" % r["peak_model"] if "peak_model" in r else "",
                        "%+.1f" % r["peak_err_pct"] if "peak_err_pct" in r else ""))
                md.append("")
    pd.DataFrame(rows).to_csv(out_dir / "summary.csv", index=False)
    if rows:
        s = pd.DataFrame(rows).groupby(["variant", "quantity"])["nrmse_pct"].agg(["mean", "median", "count"])
        md[5:5] = ["## Overall (NRMSE %, lower is better)", "", s.round(1).to_markdown(), ""]
    (out_dir / "summary.md").write_text("\n".join(md), encoding="utf-8")
    print("wrote", out_dir / "summary.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
