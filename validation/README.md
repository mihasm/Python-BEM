# Validation cases

Experimental rotors from the literature, prepared so the BEM solver can be run against them and compared.
Every case lives in `validation/cases/<case_id>/` and follows the layout below; `run_validation.py` reads it.

## Case layout

```
cases/<case_id>/
  case.json                 metadata, operating conditions, sources
  geometry.csv              blade sections
  airfoils/<name>.dat       airfoil coordinates, Selig format (name line, then x y from TE over the upper side to LE and back)
  airfoils/<name>_polars.csv  lift/drag polars
  experiment/<dataset>.csv  measured data (one file per figure/table)
  extraction/               audit trail: page renders, crops, overlays of digitised points on the original figure, scripts
  sources/                  downloaded papers/reports (git-ignored)
```

### Units and sign conventions (the solver's own)

- SI units everywhere: m, m/s, rpm, degrees, N, N·m, kg/m³, m²/s.
- `twist_deg` is the local angle between the chord line and the rotor plane, **without** the global pitch, so that
  - wind turbine: α = φ − (twist + pitch)
  - propeller:    α = (twist + pitch) − φ   (twist = blade angle β measured from the plane of rotation)
- `pitch_deg` (global blade pitch / setting) is a column in the experiment files, so a case can hold several pitch settings.
- Wind turbine coefficients: Cp = P / (½ρv³πR²), Ct = T / (½ρv²πR²), λ = ΩR/v.
- Propeller coefficients: Ct = T / (ρn²D⁴), Cp = P / (ρn³D⁵), Cq = Q / (ρn²D⁵), J = v/(nD), η = J·Ct/Cp.

### case.json

```json
{
  "id": "ntnu_bt1",
  "name": "NTNU Blind Test 1 model turbine",
  "turbine_type": 0,                    // 0 = wind/water turbine, 1 = propeller
  "B": 3,
  "R": 0.447, "Rhub": 0.045,            // m
  "fluid": {"rho": 1.2, "kin_viscosity": 1.5e-5, "medium": "air"},   // "water" → speed_of_sound 1480 m/s unless given
  "geometry_file": "geometry.csv",
  "airfoils": {"s826": {"coordinates": "airfoils/s826.dat", "polars": "airfoils/s826_polars.csv",
                         "polar_source": "experimental, DTU (Sarmast & Mikkelsen 2013), Table 3"}},
  "experiments": [
    {"file": "experiment/performance.csv", "kind": "performance",
     "x": "tsr", "y": ["cp", "ct"],
     "conditions": "v = 10 m/s constant, rpm varied; blockage-corrected? say which",
     "source": "Krogstad & Eriksen 2013, Fig. 5 (digitised)", "uncertainty": "±0.01 in Cp from digitising"}
  ],
  "sources": [{"key": "krogstad2013", "citation": "full reference", "doi": "...", "url": "..."}],
  "notes": "anything the solver needs to know: blockage, tunnel, Re range, tip shape, cylinder root, etc."
}
```

### geometry.csv

Columns: `r_m, chord_m, twist_deg, airfoil, thickness_ratio` (t/c, optional). One row per tabulated station, root to tip.
`airfoil` names a key of `case.json` → `airfoils`; `transition` is allowed between two named airfoils, like in the program.
Header lines starting with `#` say where the table came from.

### Polars (`<name>_polars.csv`)

Columns: `re, alpha_deg, cl, cd` (+ optional `cm`, `ncrit`, `source`). Long format, several Reynolds numbers in one file.
Measured polars are preferred over computed ones; a computed polar (XFOIL, Ncrit stated) is marked so in `source`.
Only the measured/attached range is stored; the runner extrapolates to ±180° (Montgomerie) itself.

### Experiment files — one fixed column dictionary

CSV, header comment lines (`#`) naming the source figure/table, how it was extracted (table / digitised) and the
estimated extraction uncertainty. **A column name always means the same quantity, in the same axes, units and sign,
in every case.** When a source uses other axes or normalisation, the case builder converts the values to these
columns (and may keep the source's own values in extra columns with other names, which the runner ignores). There
are no per-case switches.

Operating conditions (a column, or the experiment's `defaults` when constant for the whole file):

| column | meaning |
|---|---|
| `v_ms` | free-stream speed [m/s] |
| `rpm` | rotor speed [1/min] |
| `pitch_deg` | global blade pitch [°], solver convention (added to `twist_deg`) |
| `yaw_deg` | yaw angle [°] |
| `rho`, `kin_viscosity` | fluid density [kg/m³], kinematic viscosity [m²/s] |
| `tsr` | λ = ΩR/v (wind/water turbines) |
| `J` | v/(nD) (propellers) |
| `gauge_r_m` | radius [m] about which the root moments in the row are taken |

Rotor quantities:

| column | meaning |
|---|---|
| `cp`, `ct`, `cq` | power, thrust, torque coefficients (turbine or propeller definition above) |
| `eta` | propeller efficiency J·Ct/Cp |
| `power_w`, `thrust_n`, `torque_nm` | aerodynamic power [W], thrust [N], torque [N·m] |
| `root_flap_moment_nm` | bending moment of ONE blade about `gauge_r_m`, out of the rotor plane (from the thrust-direction force, downwind positive) [N·m] |
| `root_edge_moment_nm` | bending moment of ONE blade about `gauge_r_m`, in the rotor plane (from the force in the direction of rotation, positive) [N·m] |

Spanwise distributions (`r_m` or `r_R` as x, plus the operating point columns):

| column | meaning |
|---|---|
| `fn_n_per_m`, `ft_n_per_m` | force per unit span normal to the rotor plane (thrust direction) / in the rotor plane (direction of rotation) [N/m] |
| `cn_rotor`, `ct_rotor` | the same two forces as coefficients, divided by ½ρW²c |
| `cn_chord`, `ct_chord` | force coefficients normal to the chord / along the chord (positive towards the leading edge), divided by the local dynamic pressure and chord — what pressure taps give |
| `fn_chord_n_per_m`, `ft_chord_n_per_m` | the chord-axes forces per unit span [N/m] |
| `cl`, `cd`, `alpha_deg`, `phi_deg`, `a`, `a_prime`, `re` | local lift/drag coefficients, angle of attack, inflow angle, induction factors, Reynolds number |

Uncertainty: `<column>_err` (± absolute) next to any measured column — drawn as error bars. `<column>_std` may hold
the cyclic/turbulent standard deviation.

Structural (static tests): load case description, `r_m`, `deflection_m`, `strain`, `moment_nm`, … — document the
columns in the header.

### case.json → experiments

Each entry: `file`, `kind` (`performance` vs `tsr`/`J`, `power_curve` vs `v_ms`, `spanwise`, `loads`, `structural`,
`polar`), `x`, `y` (which columns to compare), `defaults` (constant conditions, same names as the columns), `source`,
`conditions` (prose), `uncertainty` (prose). For a performance file at least one of `v_ms` / `rpm` must be known per
row; the other follows from `tsr` or `J`.
