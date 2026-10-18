"""
Builds the case `structural_nrel_phase_vi` (structural data of the NREL UAE Phase VI blade) from
Hand et al. 2001, NREL/TP-500-29955 (sources/hand2001_NREL-TP-500-29955.pdf):
  Appendix A p. 75 (pdf p. 85): blade mass and centre of gravity (measured)
  Table A-9 p. 76 (pdf p. 86): estimated blade structural properties
  Table A-10 p. 79 (pdf p. 89): blade modal frequencies (measured)
  p. 34 (pdf p. 44): root strain-gauge station
The text of these pages as given by pdftotext is kept in extraction/hand2001_*_pdftotext.txt; values are parsed from it.
Run: python3 validation/cases/structural_nrel_phase_vi/extraction/build_case.py
"""
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.dirname(HERE)

rows = []
with open(os.path.join(HERE, "hand2001_tableA-9_pdftotext.txt")) as f:
    for line in f:
        tok = line.split()
        if len(tok) == 10:
            try:
                rows.append([float(t) for t in tok])
            except ValueError:
                pass
t = pd.DataFrame(rows, columns=["r_m", "mass_per_length_kg_m", "edge_inertia_kg_m", "flap_inertia_kg_m",
                                "cg_le_m", "ea_le_m", "gj_nm2", "ea_n", "ei_edge_nm2", "ei_flap_nm2"])
assert len(t) == 26, len(t)
t = t.sort_values("r_m").reset_index(drop=True)
# one convention for all structural cases: offsets from the pitch axis positive towards the TRAILING edge
t["cg_offset_m"] = -t.pop("cg_le_m") + 0.0
t["ea_offset_m"] = -t.pop("ea_le_m") + 0.0
hdr = """# NREL UAE Phase VI blade - ESTIMATED distributed structural properties (kind: structural; computed by NREL, not measured)
# Source: Hand et al. 2001, NREL/TP-500-29955, Appendix A, Table A-9 'Estimated blade structural properties' (p. 76,
#   pdf p. 86), tabulated (pdftotext, extraction/hand2001_tableA-9_pdftotext.txt), rows sorted root -> tip.
# The report says: estimates made 'with an experimental code that computes blade properties using the material properties',
#   'validated with hand calculations for the blade portion', but 'a satisfactory match between blade frequencies predicted
#   with the structural estimates and the measured frequencies has not been achieved' (pitch-shaft modelling).
# Blade construction: carbon-fibre D-spar (tapering thickness), fibreglass skin, aramid honeycomb, PVC tube (Fig. A-9) -
#   composite, so no single E/rho; r < 1.257 m is the steel pitch shaft / aluminium root insert / cylinder (axisymmetric
#   rows with zero offsets).
# Columns (structural column dictionary, validation/structural_sources.md): r_m distance from the rotor axis;
#   mass_per_length_kg_m; edge_inertia_kg_m = mass moment of inertia per length about the out-of-plane (y) axis
#   [source header says kgm4 - read as kg*m (kg m^2 per m), consistent with the magnitudes]; flap_inertia_kg_m about the
#   chordwise (x) axis; gj_nm2; ea_n [source header says Nm - read as N]; ei_edge_nm2; ei_flap_nm2;
#   cg_offset_m / ea_offset_m: chordwise distance of the section mass centre / elastic axis from the pitch-and-twist axis
#   (30 % chord), POSITIVE TOWARD THE TRAILING EDGE (the source gives them positive toward the leading edge - sign flipped).
"""
with open(os.path.join(CASE, "experiment", "distributed_properties.csv"), "w") as f:
    f.write(hdr)
    t.to_csv(f, index=False, float_format="%.6g")

# integrated values from the table (trapezoid; the table is not evenly spaced and stops at r = 5.0 m, tip 5.029 m)
m = np.trapezoid(t.mass_per_length_kg_m, t.r_m)
m1 = np.trapezoid(t.mass_per_length_kg_m * t.r_m, t.r_m)
blade = t[t.r_m >= 0.508]  # outboard of the blade root flange (Rhub 0.508 m in the aero case)
mb = np.trapezoid(blade.mass_per_length_kg_m, blade.r_m)
mb1 = np.trapezoid(blade.mass_per_length_kg_m * blade.r_m, blade.r_m)
integ = pd.DataFrame([
    dict(quantity="blade_mass_kg", value=60.2, table_a9_trapezoid=mb, measured=True,
         note="blade 3 with standard tip, probe plugs inserted (Hand 2001 p. 75); table integral from r = 0.508 m "
              "(blade root) to 5.0 m"),
    dict(quantity="cm_from_rotor_axis_m", value=2.266, table_a9_trapezoid=mb1 / mb, measured=True,
         note="blade 3 with standard tip (p. 75); table: same range"),
    dict(quantity="mass_r0.369_to_5.0_kg", value=np.nan, table_a9_trapezoid=m, measured=False,
         note="whole Table A-9 range incl. pitch shaft / root hardware"),
    dict(quantity="pitch_shaft_hardware_per_blade_kg", value=38.6, table_a9_trapezoid=np.nan, measured=True,
         note="pitch shaft, bull gear, instrumentation, bearings, nut, spacers, one blade (p. 75)"),
    dict(quantity="tip_mass_difference_tip_plate_kg", value=2.41, table_a9_trapezoid=np.nan, measured=True,
         note="tip plate minus standard tip (p. 75)"),
])
with open(os.path.join(CASE, "experiment", "integrated_properties.csv"), "w") as f:
    f.write("# NREL UAE Phase VI blade - measured mass and centre of gravity (kind: structural)\n"
            "# value: Hand et al. 2001, Appendix A p. 75 (pdf p. 85), tabulated text; measured=True for weighed values.\n"
            "# table_a9_trapezoid: trapezoidal integral of Table A-9 (estimated) by extraction/build_case.py, for\n"
            "#   comparison - the estimate and the weighing are independent.\n")
    integ.to_csv(f, index=False, float_format="%.6g")

modes = pd.DataFrame([
    dict(mode="edge_1", frequency_hz=9.062, blade="non-instrumented"), dict(mode="edge_1", frequency_hz=8.9, blade="instrumented"),
    dict(mode="flap_1", frequency_hz=7.313, blade="non-instrumented"), dict(mode="flap_1", frequency_hz=7.25, blade="instrumented"),
    dict(mode="flap_2", frequency_hz=30.062, blade="non-instrumented"), dict(mode="flap_2", frequency_hz=29.438, blade="instrumented"),
])
with open(os.path.join(CASE, "experiment", "modal_frequencies.csv"), "w") as f:
    f.write("# NREL UAE Phase VI blade - measured modal frequencies (kind: structural; for a future modal/deflection solver)\n"
            "# Source: Hand et al. 2001, Table A-10 (p. 79, pdf p. 89), tabulated. Blade on a dummy pitch shaft clamped at the\n"
            "#   bearing locations on a concrete block, trailing edge down, impact-hammer test (p. 77, Fig. A-10). Non-rotating.\n")
    modes.to_csv(f, index=False)
print(t.to_string())
print(integ.to_string())

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
ax[0].plot(t.r_m, t.mass_per_length_kg_m, "o-")
ax[0].set_xlabel("r [m]")
ax[0].set_ylabel("mass per length [kg/m]")
ax[0].set_title("Phase VI Table A-9 (estimated); blade mass 60.2 kg measured, table %.1f kg" % mb)
for c in ("ei_flap_nm2", "ei_edge_nm2", "gj_nm2", "ea_n"):
    ax[1].semilogy(t.r_m, t[c], "o-", ms=3, label=c)
ax[1].legend()
ax[1].set_xlabel("r [m]")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "check_structural_nrel_phase_vi.png"), dpi=100)
