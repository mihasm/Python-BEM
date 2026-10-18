# Structural (static) validation data — what exists and what it can check

Collected 2026-10-18 for the structural part of the program: `Calculator.statical_analysis`
(cantilever bending moments `Ms_n`/`Ms_t`, extreme-fibre bending stresses, centrifugal force/stress, equivalent stress),
`Calculator.get_crossection_data` (A, Ix, Iy, Ixy of a solid or constant-wall hollow airfoil section in the rotor frame)
and `bending_inertia.py` (polygon integrals, centroid, `generate_hollow_foil`). The program does not compute deflections
yet; several files below also carry deflection/slope/stiffness data for when it does.

## Column dictionary used by all structural files

One meaning per name in every `structural_*` case (SI units). Rows are sections, stations or quantities; which columns a
file has is listed in its header.

| column | meaning |
|---|---|
| `r_m` | radius from the rotor axis [m] (beam/centrifugal files: station on the beam, root at `r0_m`, tip at `R_m`) |
| `A_m2` | cross-section area [m²] |
| `xc_m`, `yc_m` | centroid; x along the chord from the leading edge, y normal to the chord towards the suction side [m] |
| `Ix_m4` | ∫(y−yc)² dA — about the centroidal axis parallel to the chord: **flapwise** (the program's `Ix` at blade angle 0) [m⁴] |
| `Iy_m4` | ∫(x−xc)² dA — **edgewise** (the program's `Iy`) [m⁴] |
| `Ixy_m4` | ∫(x−xc)(y−yc) dA [m⁴] |
| `y_up_m`, `y_low_m`, `x_le_m`, `x_te_m` | centroid to the outermost upper / lower fibre, to the leading / trailing edge [m] |
| `J_m4` | torsion constant [m⁴] |
| `wall_m` | constant wall thickness of a hollow section (0 = solid) [m] |
| `mass_per_length_kg_m` | mass per unit span [kg/m] |
| `ei_flap_nm2`, `ei_edge_nm2`, `gj_nm2`, `ea_n` | section stiffnesses [N·m², N] |
| `flap_inertia_kg_m`, `edge_inertia_kg_m` | mass moments of inertia per unit span [kg·m] |
| `cg_offset_m`, `ea_offset_m` | chordwise offset of mass centre / elastic axis from the pitch axis, + toward the trailing edge [m] |
| `shear_flap_n`, `moment_flap_nm` | shear force / bending moment out of the rotor plane at `r_m` from all load outboard (the program's `Ms_n`) |
| `stress_flap_pa` | extreme-fibre flapwise bending stress M·y/I [Pa] (the program's `stress_norm` is in MPa) |
| `deflection_flap_m`, `slope_flap_rad` | out-of-plane deflection and slope |
| `centrifugal_force_n`, `stress_centrifugal_pa` | axial force at `r_m` from the mass outboard (the program's `dFc`) and that force / A (`stress_cent`, MPa in the program) |
| `mass_outboard_kg` | blade mass between `r_m` and the tip |
| `frequency_hz` | natural frequency of the mode named in `mode` |

Measured blade root moments are NOT structural files: they live in the aerodynamic cases as `loads` files with the
rotor-level columns `root_flap_moment_nm` / `root_edge_moment_nm` and `gauge_r_m`.

## Cases created

### `cases/structural_analytical` — exact references (highest value for the current code)
`experiment/section_properties.csv` (73 rows), `beam_closed_form.csv`, `centrifugal_closed_form.csv`; airfoil files in
`airfoils/`; generator `extraction/make_reference.py` (independent of the solver; XFOIL BEND outputs in
`extraction/xfoil_bend/`), plot `extraction/check_structural_analytical.png`.

* Solid NACA 4-digit sections 0006…0024, 2412, 4412, 4415, 4418, 6412: A, centroid, Ix, Iy, Ixy, extreme fibres. Exact
  (dense polygon, converged to 1e-8), checked against closed-form integrals of the thickness polynomial (agree to 1e-10)
  and XFOIL BEND (agree to 2e-4, XFOIL's own spline/paneling).
  **Closed forms found:** for the NACA 4-digit thickness distribution (open TE, a4 = −0.1015) A = **0.68508·t·c**
  exactly (A = 2∫y_t dx = 10·t·c²·0.0685083) and Ix = **0.039407·c·t³** for the symmetric section; with camber
  K_I = Ix/(c t (t²+h²)) stays 0.0394 (0.03936 for 4412). Drela's rule of thumb (MIT 16.01 notes, "Area and Bending
  Inertia of Airfoil Sections", eqs. 6–11) is A ≈ 0.60·c·t, I ≈ 0.036·c·t(t²+h²) "for most common airfoils" — 12 % / 9 %
  below the NACA 4-digit values; it is in the file as `*_drela_rule` rows with a 15 % band. I did not find a primary
  source for the "0.685" and "0.036 c t³" pair quoted in the brief; 0.685 is the exact NACA 4-digit value (derived here),
  0.036 is Drela's generic coefficient.
* Hollow (constant wall) NACA 0012, 0018, 2412, 4412, 4415 with walls 2/5/10/20 mm on a 1 m chord, and at the program's
  scale (c = 0.1 m, walls 1–2 mm): exact inward offset of the outline (GEOS buffer), Bredt-Batho J, plus the thin-wall
  approximation from XFOIL's skin integrals for comparison.
* Rotated sections (NACA 4412, θ = 10/30/90°), with the rotation formulas in the header.
* Circle / tube (the program's cylindrical root) and ellipse: closed forms; hollow ellipse by exact offset.
* Chen, Yu & Capellaro 2010 (Wind Energy 13:497, doi 10.1002/we.372) **tabulated benchmarks**: aluminium tube (Table I,
  elasticity solution, 5 wall ratios) and the isotropic blade-like section (Fig. 6 geometry, Table IV: VABS, ANSYS,
  thin-wall analytical, PreComp). My exact geometry reproduces their A, Ix, Iy and centroid to < 0.02 %, which validates
  the offset method used for all hollow rows.
* Cantilever closed forms (uniform, linear and tip load; moment, shear, stress, deflection, slope; Roark Table 8.1,
  checked by numerical integration) and centrifugal force/stress of a uniform and a tapered solid NACA 0012 bar.

**Validates:** `get_crossection_data` (both `blade_design` values), `calculate_bending_inertia_2`, `calculate_centroid`,
`generate_hollow_foil`, extreme-fibre distances, the `Ms_n`/`Ms_t` sums, `stress_norm`, `dFc`/`stress_cent`; deflection
columns for a future solver.

**What the program gives today (run 2026-10-18, read-only, current working tree):** solid rows within 2·10⁻⁴ (= the
161-point .dat discretisation, column `*_datfile_rel_err`); hollow symmetric sections within 0.3 %; circle tube A +0.1…+0.8 %,
Iy +0.2…+1.2 %; NACA 4412/4415 with thin walls (2–5 mm on 1 m) A −1.3…−2.6 %, Iy −2.7…−5.2 % (`generate_hollow_foil`
drops/misplaces points where the upper/lower-surface split by the gradient sign is ambiguous on cambered sections); the
blade-like section was +182 % until its straight edges were given points (the offset is point-based) — with them A +0.16 %, Ix −0.03 %, Iy +0.35 %.
Two code observations (not changed): the moment sum `Σ dFn_j (r_j − r_i)` gives the section's own element zero arm (the
half-element outboard of r_i is lost, O(dr)); the centrifugal loop lumps element mass at its centre (midpoint rule,
O(dr²)). The closed-form files state this so a comparison uses fine sections.

### `cases/structural_nrel5mw` — reference blade properties (composite; reference model, no measurement)
Jonkman et al. 2009, NREL/TP-500-38060: Table 2-1 (49 stations: mass/length, EI flap/edge, GJ, EA, mass inertias,
CG offset), Table 2-2 (blade mass 17,740 kg, first/second mass moments, CG 20.475 m) and Table 3-1 geometry
(17 elements, chord, twist, airfoils), with the model's airfoil coordinates and polars from OpenFAST r-test.
Table 2-1 is identical (to every digit) to the r-test ElastoDyn and BeamDyn blade files. Lumping Table 2-1 × 1.04536 at
the 17 element centres — the program's rule — reproduces Table 2-2 to 1e-5 (trapezoid: −0.7 %).
`derived_equivalent_homogeneous.csv` (marked derived) gives the solid-section A/I of each element, and the equivalent
density/modulus/GFRP wall that a homogeneous model would need.
**Validates:** blade mass and centrifugal root force summation; order of magnitude of the section model on a real blade;
input for a future deflection/frequency solver. Quality: tabulated, exact, but design data of a composite blade, so a
homogeneous section only matches with fitted ρ/E/wall.

### `cases/structural_nrel_phase_vi` — Phase VI blade structural data (companion of `cases/nrel_phase_vi`)
Hand et al. 2001, NREL/TP-500-29955: measured blade mass 60.2 kg and CG 2.266 m (p. 75), estimated distributed
properties Table A-9 (26 stations; integrates to 57.3 kg for r ≥ 0.508 m), measured modal frequencies Table A-10.
Composite blade (carbon D-spar, glass skin, honeycomb).
**Validates:** blade mass / first mass moment (→ centrifugal force); future modal solver.

## Item 1 — measured root bending moments paired with aerodynamics (collected by the aero workers, not duplicated)

* **NREL Phase VI:** `cases/nrel_phase_vi/experiment/loads.csv` — mean root flap/edge moments vs wind speed, sequence S,
  gauges on the 80 mm steel pitch shaft at r = 0.432 m (Hand 2001 p. 34, channels B1RFB/B1REB/B3RFB/B3REB; gauge axes
  follow the tip chord — the file converts to rotor axes). The report also has an estimated aerodynamic root flap moment
  from the pressure taps (channel 817, Eq. 15) that includes a blade-mass correction.
* **MEXICO:** aluminium blades ("the chosen blade material is aluminium; both the external and the internal shape
  numerically milled", ECN-E--07-042 p. 15), 6 strain-gauge bridges for root flap/edge moments. The first campaign's
  gauges are reported as reliable on two blades only (p. 8–9) and "did not give correct results" according to the
  New MEXICO report (ECN-E--14-048; new root gauges were installed). The `cases/mexico` folder had no experiment files yet when
  this was written (sources only). No published blade mass / section drawing was found in the downloaded reports — the
  internal milling (pressure-tube channels) means the solid-section model would overestimate the mass.

## Item 2 — static bending tests of simple homogeneous blades: nothing usable found

Searched (Google Scholar/OpenAlex, web): small-turbine static tests (IEC 61400-2), timber, 3D-printed, aluminium blades,
Sandia CX-100/TX-100/NRT, propellers. Findings:

* **Peterson & Clausen 2004**, Wind Engineering 28(1):87 (doi 10.1260/0309524041210865): 1 m solid hoop-pine / radiata
  blades (E 12.1 / 10.5 GPa, ρ 545 / 500 kg/m³, Table 1) — static tests "with strains and deflections recorded" were
  done, but **no values or geometry are published** (only used to tune an FE model). Not usable.
* **Hartman & Biermann 1938, NACA Report 644** (https://ntrs.nasa.gov/citations/19930091717): solid aluminium-alloy
  Clark Y propellers (10 ft, 3 blades), blade-form curves (b/D, h/b, p/D vs r/R, Fig. 3) and **measured bending slope
  and twist at 0.70R and 0.85R vs V/nD and rpm in operation** (Figs. 7–13, ±0.1° accuracy). The only homogeneous blade
  with a measured operating deflection I found; it needs digitising of overlapping curves, the aero data of the same
  propellers (NACA TR 640/642), and a deflection solver with centrifugal stiffening — left as a future case (PDF was
  fetched to /tmp only).
* Sandia CX-100/TX-100 static tests (Paquette/Resor, OSTI 1111747) and NRT 13 m (NREL/TP-5000-79416, OSTI 1567986):
  composite; the reports were not reachable (osti.gov timed out, docs.nlr.gov serves a captcha). Recommended next try
  from a browser.
* Alber et al. 2022 (J. Phys. Conf. Ser. 2265 042070): 3D-printed PLA shell with a CFRP spar — not homogeneous, no
  deflection data. Shah et al. 2013 (flax vs glass small blade, Compos. B): composite, only a tip deflection number.
  Kim/Lee 2016 (MSSP, wooden blade specimen): dynamic, no static table.

## Item 3/4 references not reached
* PreComp user guide (Bir 2006, NREL/TP-500-38929) — captcha/timeout; Chen et al. 2010 contains the same tube/blade-like
  benchmarks with PreComp results, so nothing is lost.
* WindPACT/SNL blade property tables — not attempted after the NREL/OSTI servers failed.
