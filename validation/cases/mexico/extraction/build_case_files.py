"""Assembles the case files of validation/cases/mexico from the extraction products in this folder.

Run order: extract_mexnext3_figA4_A5.py, extract_risoe_A1_21_LER.py, digitise_raster_markers.py,
make_airfoil_coordinates.py, run_xfoil_polars.py, then this script.
"""
import os, math
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.abspath(os.path.join(HERE, '..'))
EXP = os.path.join(CASE, 'experiment')
AIR = os.path.join(CASE, 'airfoils')
os.makedirs(EXP, exist_ok=True)
R = 2.25
geo = pd.read_csv(os.path.join(CASE, 'geometry.csv'), comment='#')


def twist_at(r):
    return float(np.interp(r, geo.r_m, geo.twist_deg))


def chord_at(r):
    return float(np.interp(r, geo.r_m, geo.chord_m))


def write(path, header_lines, df, float_fmt='%.5g'):
    with open(path, 'w') as f:
        for h in header_lines:
            f.write('# ' + h + '\n')
        df.to_csv(f, index=False, float_format=float_fmt)
    print('wrote', os.path.relpath(path, CASE), len(df), 'rows')


def chord_to_rotor(fn_c, ft_c, theta_deg):
    """chord axes (normal to chord, downwind positive; along chord, TE->LE positive) -> rotor-plane axes
    (thrust direction; direction of rotation). theta = local blade angle = twist + pitch (chord vs rotor plane,
    leading edge turned upwind for theta > 0). Verified against the Mexnext-III weighting factors
    (ECN-E--18-003 App. C, Tables 1-2): Fax weights Ft/Fn = -tan(theta), torque weights Fn/Ft = tan(theta)."""
    t = np.radians(theta_deg)
    fn = fn_c * np.cos(t) - ft_c * np.sin(t)
    ft = fn_c * np.sin(t) + ft_c * np.cos(t)
    return fn, ft


# =====================================================================================================
# 1. Polars
# =====================================================================================================
def xfoil_csv(name, mode):
    d = pd.read_csv(os.path.join(HERE, 'xfoil', '%s_%s_all.csv' % (name, mode)))
    d['re'] = d['re'].astype(float)
    return d[['re', 'alpha_deg', 'cl', 'cd', 'cm', 'ncrit', 'source']]


hdr_x = ['XFOIL 6.97 (repository binary xfoil_executables/xfoil), 160 panels, ITER 200, Ncrit = 9, '
         'alpha -10..+20 deg in 0.5 deg steps from 0 deg both ways; non-converged points are absent.',
         'Coordinates: airfoils/%s.dat. Script: extraction/run_xfoil_polars.py, raw .pol/.log files in extraction/xfoil/.',
         'Re = 3e5, 5e5, 7e5, 1e6 cover the rotor: chord Re at 424.5 rpm is ~0.45e6 (r=0.56 m) ... ~0.75e6 (r=1.85 m).']
# DU 91-W2-250: forced transition at 5% (the blades carry zig-zag tape inboard in every MEXICO configuration)
du = xfoil_csv('du91_w2_250', 'trip')
write(os.path.join(AIR, 'du91_w2_250_polars.csv'),
      ['DU 91-W2-250 polars, COMPUTED (the measured Mexnext set, TU Delft LTT Re=0.5e6 with zig-zag tape, is not public).',
       'Forced transition x/c = 0.05 on both sides (XTR 0.05 0.05) = the zig-zag tape at 5% chord on the MEXICO blades.']
      + [h % 'du91_w2_250' if '%s' in h else h for h in hdr_x], du)
du_free = xfoil_csv('du91_w2_250', 'free')
write(os.path.join(AIR, 'du91_w2_250_polars_xfoil_free.csv'),
      ['ALTERNATIVE set (not referenced by case.json polars): DU 91-W2-250, XFOIL free transition Ncrit = 9.']
      + [h % 'du91_w2_250' if '%s' in h else h for h in hdr_x], du_free)
# measured clean DU at Re 1e6 (TU Delft, Timmer), tabulated in the Risø airfoil catalogue data set
m = pd.read_csv(os.path.join(HERE, 'raw_airfoil_data', 'du91w2250_cldm_tudelft.exp'), sep=r'\s+', comment='#',
                names=['alpha_deg', 'cl', 'cd', 'mcm'])
m['cm'] = -m['mcm']
m.insert(0, 're', 1.0e6)
m['source'] = 'measured, TU Delft LTT (Timmer), clean, Re=1e6; Risoe-R-1280 catalogue data DU91-W2-250/cldm.exp'
write(os.path.join(AIR, 'du91_w2_250_polars_measured_clean_re1e6.csv'),
      ['ALTERNATIVE set (not referenced by case.json polars): DU 91-W2-250 MEASURED, CLEAN (no tape), Re = 1.0e6,',
       'TU Delft low-speed low-turbulence tunnel (W.A. Timmer). Tabulated values from the data set of the Risoe Wind',
       'Turbine Airfoil Catalogue (Bertagnolio, Sorensen, Johansen, Fuglsang 2001, Risoe-R-1280(EN)), file',
       'DU91-W2-250/cldm.exp, mirrored at github.com/balagos/mjl-airfoil-data/RISO-data (copy: extraction/raw_airfoil_data/).',
       'The source column is "-cm .25c"; cm here = -(source), nose-up positive (cm ~ -0.12 at 0 deg).',
       'Clean and Re = 1e6 differ from the MEXICO blade (tripped, Re 0.4-0.6e6).'],
      m[['re', 'alpha_deg', 'cl', 'cd', 'cm', 'source']])

# NACA 64(3)-418: main = free transition (New MEXICO configuration 3 has the outboard 30% clean; Mexnext-III
# prescribes the CLEAN Re=0.7e6 data set for it); tripped set as alternative (2006 MEXICO and New MEXICO config 0)
na = xfoil_csv('naca64_418', 'free')
write(os.path.join(AIR, 'naca64_418_polars.csv'),
      ['NACA 64(3)-418 polars, COMPUTED (the measured Mexnext set, TU Delft Re=0.7e6 clean, is not public).',
       'FREE transition, Ncrit = 9: matches New MEXICO 2014 configuration 3 (zig-zag tape removed for r/R > 0.7), the',
       'configuration of the main data sets (Mexnext-III axial cases 1.1-1.3, clean TSR sweep). For the fully tripped',
       'configurations (MEXICO 2006, New MEXICO config 0) use naca64_418_polars_xfoil_trip.csv.']
      + [h % 'naca64_418' if '%s' in h else h for h in hdr_x], na)
na_t = xfoil_csv('naca64_418', 'trip')
write(os.path.join(AIR, 'naca64_418_polars_xfoil_trip.csv'),
      ['ALTERNATIVE set (not referenced by case.json polars): NACA 64(3)-418, XFOIL forced transition x/c = 0.05',
       'both sides (zig-zag tape as in MEXICO 2006 and New MEXICO configuration 0).']
      + [h % 'naca64_418' if '%s' in h else h for h in hdr_x], na_t)

# RISØ-A1-21: measured with zig-zag tape (LER) at Re 1.6e6, digitised exactly from the vector figure
ri = pd.read_csv(os.path.join(HERE, 'risoe_A1_21_LER_raw.csv'))
ri = ri[(ri.alpha_deg >= -10) & (ri.alpha_deg <= 30.5)].copy()
ri.insert(0, 're', 1.6e6)
ri['cm'] = ri['cm']
ri['source'] = 'measured, VELUX tunnel, zigzag tape (LER), Re=1.6e6; Risoe-R-1112 Fig. 6-8 (vector extraction)'
write(os.path.join(AIR, 'risoe_a1_21_polars.csv'),
      ['RISO-A1-21 polar, MEASURED: Fuglsang, Dahl & Antoniou 1999, "Wind tunnel tests of the Risoe-A1-18, Risoe-A1-21',
       'and Risoe-A1-24 airfoils", Risoe-R-1112(EN), section 6.2, run005 "LER" = 90-deg zigzag trip tape (3 mm wide,',
       '0.35 mm thick) at the leading edge; VELUX open-jet tunnel, Re = 1.6e6, chord 0.6 m, tunnel-corrected.',
       'This is the Re = 1.6e6 "rough" data set that Mexnext prescribes for RISO-A1-21 (ECN-E--18-003 App. C).',
       'Extraction: the figure (Fig. 6-8, PDF page 39) is vector graphics; the filled-circle markers were read from',
       'the PDF paths with PyMuPDF and the axes calibrated on the gridlines (residual < 0.001 in all three panels):',
       'extraction/extract_risoe_A1_21_LER.py, overlay extraction/overlay_risoe_fig6-8.png. Uncertainty from',
       'extraction: < 0.01 deg, < 0.002 in cl, < 0.0005 in cd. cm about 0.25c, nose-up positive (source convention).',
       'Reynolds number of the rotor at the RISO section is only ~0.6-0.7e6; this data set is at 1.6e6 (no measured',
       'low-Re data exist publicly; the coordinates are proprietary, so no XFOIL set was made either).'],
      ri[['re', 'alpha_deg', 'cl', 'cd', 'cm', 'source']])
smooth = pd.read_csv(os.path.join(HERE, 'raw_airfoil_data', 'risoA121_cldm_smooth_velux.exp'), sep=r'\s+',
                     comment='#', names=['alpha_deg', 'cl', 'cd', 'cm_source'])
smooth.insert(0, 're', 1.6e6)
smooth['source'] = 'measured, VELUX, smooth, Re=1.6e6; Risoe-R-1280 catalogue data RISO-A1-21/cldm.exp'
write(os.path.join(AIR, 'risoe_a1_21_polars_measured_smooth.csv'),
      ['ALTERNATIVE set (not referenced by case.json polars): RISO-A1-21 MEASURED SMOOTH (no tape), VELUX tunnel,',
       'Re = 1.6e6, run025 of Risoe-R-1112; tabulated in the Risoe Wind Turbine Airfoil Catalogue data set',
       '(Risoe-R-1280, RISO-A1-21/cldm.exp; copy in extraction/raw_airfoil_data/). Columns as in the source',
       '(aoa_c, cl, cd_c); cm_source = the catalogue cm_c column, sign convention not stated there (+0.07..0.15,',
       'while Risoe-R-1112 plots CM negative) - therefore not stored as cm.'], smooth)

# =====================================================================================================
# 2. New MEXICO 2014, Mexnext-III axial cases 1.1-1.3: spanwise loads + integrated rotor loads
# =====================================================================================================
A4 = pd.read_csv(os.path.join(HERE, 'mexnext3_figA4_raw.csv'), dtype={'case': str})
A4 = A4[A4.r_m < 2.25]
cond = {'1.1': (10.05, 1.197), '1.2': (15.06, 1.191), '1.3': (24.05, 1.195)}
RPM_C, PITCH = 425.1, -2.3
rows = []
for c, (v, rho) in cond.items():
    fn = A4[(A4.case == c) & (A4.quantity == 'fn')].sort_values('r_m')
    ft = A4[(A4.case == c) & (A4.quantity == 'ft')].sort_values('r_m')
    for (_, a), (_, b) in zip(fn.iterrows(), ft.iterrows()):
        r = round(a.r_m, 4)
        th = twist_at(r) + PITCH
        fnc, ftc = a.value, b.value
        fnc_e = 0.5 * (a.band_hi - a.band_lo)
        ftc_e = 0.5 * (b.band_hi - b.band_lo)
        fn_r, ft_r = chord_to_rotor(fnc, ftc, th)
        t = math.radians(th)
        fn_r_e = math.hypot(fnc_e * math.cos(t), ftc_e * math.sin(t))
        ft_r_e = math.hypot(fnc_e * math.sin(t), ftc_e * math.cos(t))
        rows.append(dict(case=c, r_m=r, r_R=round(r / R, 4), v_ms=v, rpm=RPM_C, pitch_deg=PITCH, yaw_deg=0.0, rho=rho,
                         tsr=round(RPM_C * math.pi / 30 * R / v, 3),
                         fn_n_per_m=fn_r, fn_n_per_m_err=fn_r_e, ft_n_per_m=ft_r, ft_n_per_m_err=ft_r_e,
                         fn_chord_n_per_m=fnc, fn_chord_n_per_m_err=fnc_e, ft_chord_n_per_m=ftc,
                         ft_chord_n_per_m_err=ftc_e, local_blade_angle_deg=round(th, 3)))
span14 = pd.DataFrame(rows)
write(os.path.join(EXP, 'newmexico2014_spanwise_axial.csv'), [
    'New MEXICO (DNW-LLF, 2014), Mexnext-III comparison cases 1.1 / 1.2 / 1.3 (axial flow, yaw 0, pitch -2.3 deg,',
    '425.1 rpm, V = 10.05 / 15.06 / 24.05 m/s, rho = 1.197 / 1.191 / 1.195 kg/m3; ECN-E--18-003 Table 4.1 and App. C.1).',
    'Model configuration 3: zig-zag tape (5% chord, both sides) on r/R < 0.7, outboard 30% clean.',
    'Measured: azimuth-averaged sectional forces from the pressure taps at r/R = 0.25, 0.35, 0.60, 0.82, 0.92, averaged',
    'over several datapoints (Table D.1); pressure-only (no skin friction), TE pressure = mean of last two sensors.',
    'Source: Boorsma, Schepers et al. 2018, "Final report of IEA Wind Task 29 Mexnext (Phase 3)", ECN-E--18-003,',
    'Figure A.4 (p. 149), black diamonds = "NewMexico"; *_err = half width of the grey "stdev(avg)" band',
    '(repeatability between datapoints), NOT the sensor accuracy (see ECN-E--14-048 Table 1: up to +-0.46 in Cn at',
    '25%R, +-0.03 at 92%R, worst case).',
    'Extraction: exact, from the PDF vector paths (PyMuPDF), axes calibrated on tick marks, residual < 0.006 N/m;',
    'extraction/extract_mexnext3_figA4_A5.py, overlay_mexnext3_figA4.png. Extraction uncertainty < 0.1 N/m.',
    'Check: integrating these values with the Mexnext weighting factors (App. C Tables 1-2) gives Fax = 974 / 1663 N',
    'for cases 1.1 / 1.2 = Table 12.3 "New MEXICO (2014)" exactly (974 / 1663 N).',
    'fn_chord_n_per_m / ft_chord_n_per_m = the source quantities (normal to chord, downwind positive / along chord,',
    'TE->LE positive). fn_n_per_m / ft_n_per_m = converted to rotor-plane axes (thrust direction / direction of',
    'rotation) with the local blade angle theta = twist(r) + pitch, local_blade_angle_deg column:',
    'fn = Fn_c cos(theta) - Ft_c sin(theta), ft = Fn_c sin(theta) + Ft_c cos(theta).'], span14)

# rotor integrals: Fig A.5 Cdax, Cp (and recomputed from the spanwise values with the weighting factors)
A5 = pd.read_csv(os.path.join(HERE, 'mexnext3_figA5_raw.csv'), dtype={'case': str})
wfn_ax = [0.847477, 1.170039, 1.584740, 1.079999, 0.607397]
wft_ax = [-0.179364, -0.162356, -0.069191, -0.001508, 0.011208]
wfn_q = [0.089458, 0.155253, 0.091073, 0.002579, -0.022949]
wft_q = [0.422679, 1.118850, 2.085914, 1.846798, 1.243645]
rows = []
A = math.pi * R ** 2
om = RPM_C * math.pi / 30
for c, (v, rho) in cond.items():
    s = span14[span14.case == c].sort_values('r_m')
    fax = np.dot(wfn_ax, s.fn_chord_n_per_m) + np.dot(wft_ax, s.ft_chord_n_per_m)
    q = np.dot(wfn_q, s.fn_chord_n_per_m) + np.dot(wft_q, s.ft_chord_n_per_m)
    cdax_fig = float(A5[(A5.case == c) & (A5.quantity == 'cdax')].value.iloc[0])
    cp_fig = float(A5[(A5.case == c) & (A5.quantity == 'cp')].value.iloc[0])
    rows.append(dict(case=c, v_ms=v, rpm=RPM_C, pitch_deg=PITCH, yaw_deg=0.0, rho=rho, tsr=round(om * R / v, 3),
                     thrust_n=fax, torque_nm=q, power_w=q * om,
                     ct=fax / (0.5 * rho * v ** 2 * A), cp=q * om / (0.5 * rho * v ** 3 * A),
                     ct_figA5=cdax_fig, cp_figA5=cp_fig))
integ = pd.DataFrame(rows)
write(os.path.join(EXP, 'newmexico2014_axial_rotor_loads.csv'), [
    'New MEXICO 2014, Mexnext-III axial cases 1.1-1.3 (conditions as in newmexico2014_spanwise_axial.csv).',
    'thrust_n / torque_nm = rotor axial force and aerodynamic torque from the PRESSURE sections, integrated with the',
    'Mexnext weighting factors (ECN-E--18-003 App. C, Tables 1 and 2; linear between the 5 sections, zero load at',
    'root and tip) applied to the Fn/Ft of newmexico2014_spanwise_axial.csv; ct = Cdax, cp as defined in the report.',
    'ct_figA5 / cp_figA5 = the values plotted in Fig. A.5 (p. 150, vector extraction): they agree to < 0.3%.',
    'Cross-check: Table 12.3 lists Fax 974 / 1663 N (same) and torque 68 / 301 N m (ours 68 / 317 N m; Fig. A.5 Cp',
    'agrees with ours, so Table 12.3 torque for case 1.2 comes from another reduction).',
    'CAUTION when comparing to a BEM rotor total: this integration (5 sections, linear, zero at root/tip) misses',
    'load: ECN-E--18-003 Table 3.1 estimates -5.8/-5.1/-2.5 % on Fax and -10.5/-5.6/+0.1 % on torque for cases',
    '1.1/1.2/1.3 from the radial spacing, and -24.5/-6.4/-1.8 % on torque from the missing skin friction (the',
    'pressure torque is too HIGH by that amount). The fair comparison is to apply the same weighting factors to the',
    'model Fn/Ft at the 5 radii (factors in case.json notes); the spanwise file is the primary data.'], integ)

# =====================================================================================================
# 3. New MEXICO 2014 TSR sweeps (424 rpm, pitch -2.3), balance / pressure Cdax, Cp
# =====================================================================================================
def pts(name, series):
    d = pd.read_csv(os.path.join(HERE, name + '_points.csv'))
    return d[d.series == series][['x_value', 'y_value']].rename(columns={'x_value': 'tsr', 'y_value': 'val'})


RPM_S = 424.5


def merge_series(parts):
    """parts: list of (column, DataFrame tsr/val) -> long frame with one value column filled per row."""
    out = []
    for col, d in parts:
        d = d.copy()
        d[col] = d.pop('val')
        out.append(d)
    m = pd.concat(out, ignore_index=True).sort_values('tsr')
    m['rpm'] = RPM_S
    m['pitch_deg'] = PITCH
    m['v_ms'] = RPM_S * math.pi / 30 * R / m['tsr']
    return m


# rough = configuration 0 (tape over the full span)
cp_rough_p = pts('m16063_fig4_0', 'rough (pressure)')
cp_rough_p = cp_rough_p[~((cp_rough_p.tsr > 7.0) & (cp_rough_p.tsr < 7.05))]          # line-crossing artefact
manual = pd.DataFrame({'tsr': [5.05, 5.27, 6.48, 6.73], 'val': [0.347, 0.364, 0.463, 0.467]})  # hidden under red
cp_rough_p = pd.concat([cp_rough_p, manual])
cp_rough_g = pts('m16063_fig4_0', 'rough (generator)')
cp_rough_g = cp_rough_g[~(((cp_rough_g.tsr < 4.12) & (cp_rough_g.val > 0.3)) | (cp_rough_g.val < 0.15))]
cp_rough_g = cp_rough_g.groupby((cp_rough_g.tsr.diff().abs() > 0.06).cumsum()).mean()       # merge split rings
rough = merge_series([('ct', pts('e14048_fig22_0', 'New MEXICO balance config 0')),
                      ('ct_pressure', pts('e14048_fig22_1', 'New MEXICO pressure config 0')),
                      ('cp', cp_rough_p), ('cp_generator', cp_rough_g)])
hdr_sweep = [
    'TSR sweep in axial flow, 424 rpm (stored as 424.5), pitch -2.3 deg, yaw 0; v_ms = Omega R / tsr (derived).',
    'ct = Cdax from the BALANCE: total axial force minus the estimated tower and nacelle drag (ECN-E--14-048 p. 45),',
    '     i.e. the full rotor thrust; ECN-E--14-048 (Boorsma & Schepers 2014, "New MEXICO experiment, preliminary',
    '     overview with initial validation") Fig. %s(a).',
    'ct_pressure = Cdax from the pressure sections (5 sections, linear, zero at root/tip, +6%% for the non-linear',
    '     root/tip load: ECN-E--14-048 p. 45-46), Fig. %s(b). Not compared by the runner (extra column).',
    'cp = Cp from the PRESSURE torque (pressure only, no skin friction -> too high), ECN-M--16-063 Fig. 4(b) filled',
    '     markers (Boorsma & Schepers 2016, "Rotor experiments in controlled conditions continued: New MEXICO",',
    '     J. Phys.: Conf. Ser. 753 022004). cp_generator = Cp from the generator torque, corrected for generator losses',
    '     by the authors (open markers; includes friction; probably still contains drivetrain losses) - extra column.',
    '     The true aerodynamic Cp lies between the two; the authors give no further reduction.',
    'Extraction: DIGITISED from raster (Excel) charts embedded in the PDFs: colour masks + marker cores, axes from',
    'gridlines; extraction/digitise_raster_markers.py, overlays *_overlay.png. Uncertainty +-0.05 in tsr,',
    '+-0.006 in Cdax, +-0.003 in Cp; points hidden under markers of the other series are missing (4 rough Cp',
    'pressure points were read by eye from the overlay: tsr 5.05/5.27/6.48/6.73, +-0.005).',
    'Tunnel speed: New MEXICO uses the re-calibrated DNW tunnel speed (correct).',
    'Rows hold one measured quantity each where the series were measured at different tsr; empty = not given.']
write(os.path.join(EXP, 'newmexico2014_tsr_sweep_rough.csv'),
      ['New MEXICO 2014, model configuration 0 = zig-zag tape (5% chord, both sides) on the FULL span.']
      + [h.replace('%s', '22').replace('%%', '%') for h in hdr_sweep],
      rough[['tsr', 'v_ms', 'rpm', 'pitch_deg', 'ct', 'ct_pressure', 'cp', 'cp_generator']])

cp_clean_p = pts('m16063_fig4_0', 'clean (pressure)')
gen_low = cp_clean_p[(cp_clean_p.tsr < 4.2) & (cp_clean_p.val < 0.24)]                     # open ring filled by overlap
cp_clean_p = cp_clean_p.drop(gen_low.index)
cp_clean_g = pd.concat([pts('m16063_fig4_0', 'clean (generator)'), gen_low])
cp_clean_g = cp_clean_g.sort_values('tsr')
cp_clean_g = cp_clean_g.groupby((cp_clean_g.tsr.diff().abs() > 0.06).cumsum()).mean()
clean = merge_series([('ct', pts('e14048_fig30_0', 'New MEXICO balance config 3 (clean outboard)')),
                      ('ct_pressure', pts('e14048_fig30_1', 'New MEXICO pressure config 3 (clean outboard)')),
                      ('cp', cp_clean_p), ('cp_generator', cp_clean_g)])
write(os.path.join(EXP, 'newmexico2014_tsr_sweep_clean.csv'),
      ['New MEXICO 2014, model configuration 3 = tape on r/R < 0.7 only, outboard 30% (NACA 64-418) CLEAN -',
       'the configuration of the Mexnext-III axial cases.']
      + [h.replace('%s', '30').replace('%%', '%') for h in hdr_sweep],
      clean[['tsr', 'v_ms', 'rpm', 'pitch_deg', 'ct', 'ct_pressure', 'cp', 'cp_generator']])

# =====================================================================================================
# 4. MEXICO 2006
# =====================================================================================================
m06 = merge_series([('ct', pts('e14048_fig22_0', 'MEXICO 2006 balance')),
                    ('ct_pressure', pts('e14048_fig22_1', 'MEXICO 2006 pressure'))])
write(os.path.join(EXP, 'mexico2006_tsr_sweep.csv'), [
    'MEXICO 2006 (first campaign, DNW-LLF, December 2006), axial flow, 424.5 rpm, pitch -2.3 deg, full-span zig-zag tape.',
    'Source: ECN-E--14-048 Fig. 22 (a) balance, (b) pressures - black squares "MEXICO". Digitised (raster), see',
    'newmexico2014_tsr_sweep_rough.csv for method and uncertainty (+-0.05 tsr, +-0.006 Cdax).',
    'CORRECTED values as plotted there: the 2006 tunnel speed is known to be over-indicated (DNW calibration error',
    'found in Mexnext / New MEXICO); in this figure 0.2 m/s was deducted from the 2006 speeds before computing tsr',
    'and Cdax (ECN-E--14-048 p. 42); the 2006 pressures here use the Nov. 2013 sensor re-calibration (+5..7 % loads,',
    'ECN-E--14-048 Table 2). ct = balance Cdax minus estimated tower/nacelle drag; ct_pressure = pressure Cdax',
    '(extra column). v_ms derived = Omega R / tsr. No Cp: the 2006 campaign had no torque measurement.'],
    m06[['tsr', 'v_ms', 'rpm', 'pitch_deg', 'ct', 'ct_pressure']])

# Bechmann et al. 2011 Table IV (2006 run 11 pnt 92/93/95) Cn, Ct at the 5 sections, raw 2006 processing
b4 = {9.97: [(0.25, 0.5593, -0.0136), (0.35, 0.6787, -0.0254), (0.60, 0.4972, -0.0092), (0.82, 0.5380, -0.0181),
             (0.92, 0.5475, -0.0189)],
      14.93: [(0.25, 1.1387, -0.1794), (0.35, 1.1700, -0.1065), (0.60, 0.9250, -0.1007), (0.82, 0.8433, -0.0783),
              (0.92, 0.8168, -0.0786)],
      24.02: [(0.25, 1.8972, -0.5085), (0.35, 1.5603, -0.2867), (0.60, 1.2430, -0.2794), (0.82, 1.0278, -0.1617),
              (0.92, 0.9315, -0.1475)]}
rho06 = {9.97: 1.19, 14.93: 1.19, 24.02: 1.20}
OM06 = 44.44
rows = []
for v, lst in b4.items():
    for rR, cn, ctb in lst:
        r = rR * R
        c = chord_at(r)
        q = 0.5 * rho06[v] * (v ** 2 + (OM06 * r) ** 2)
        fnc = cn * q * c
        ftc = -ctb * q * c
        th = twist_at(r) + PITCH
        fn_r, ft_r = chord_to_rotor(fnc, ftc, th)
        rows.append(dict(r_m=round(r, 4), r_R=rR, v_ms=round(v - 0.2, 2), v_tunnel_raw_ms=v, rpm=424.5,
                         pitch_deg=PITCH, yaw_deg=0.0, rho=rho06[v],
                         fn_n_per_m=fn_r, ft_n_per_m=ft_r, fn_chord_n_per_m=fnc, ft_chord_n_per_m=ftc,
                         cn_source=cn, ct_source=ctb, chord_m=c, local_blade_angle_deg=round(th, 3)))
write(os.path.join(EXP, 'mexico2006_spanwise_axial.csv'), [
    'MEXICO 2006, run 11 datapoints 92 / 93 / 95: axial flow, pitch -2.3 deg, 424.5 rpm (paper: 425.5 rpm and',
    'Omega = 44.44 rad/s, which is 424.4 rpm), full-span zig-zag tape.',
    'Source: Bechmann, Sorensen & Zahle 2011, "CFD simulations of the MEXICO rotor", Wind Energy 14:677-689,',
    'doi:10.1002/we.450, Table IV (TABULATED): cn_source, ct_source = normal / tangential force coefficients from the',
    'pressure distributions, normalised by 0.5 rho (W^2 + (Omega r)^2) c (paper eq. 2, no induction), ct_source',
    'positive towards the TRAILING edge (negative = driving). Tunnel speeds as given there (9.97/14.93/24.02 m/s) are',
    'the RAW 2006 values -> column v_tunnel_raw_ms; v_ms = raw - 0.2 m/s (correction used in ECN-E--14-048 p. 42;',
    'Mexnext-III found the 2006 tunnel speed wrongly calibrated). The 2006 pressures are NOT re-calibrated here',
    '(the 2013 calibration raises Fn by +4..6 %, Ft by +1..15 %, ECN-E--14-048 Table 2).',
    'Conversion to forces: fn_chord = cn_source*q*c, ft_chord = -ct_source*q*c (TE->LE positive), q with the raw',
    'speed, Omega = 44.44 rad/s and rho = 1.19/1.19/1.20 kg/m3 (MEXICO densities in ECN-E--14-048 Fig. 19-21);',
    'rotor-plane fn/ft as in newmexico2014_spanwise_axial.csv. Uncertainty: tabulated to 4 digits; the',
    'normalisation assumption (rho, Omega) adds ~1 %. Sensor accuracy: ECN-E--14-048 Table 1.'],
    pd.DataFrame(rows))

b3 = pd.DataFrame(dict(v_tunnel_raw_ms=[9.97, 14.93, 24.02], axial_force_balance_n=[1106, 2051, 3315],
                       tower_drag_est_n=[174, 390, 1010], thrust_n=[932, 1661, 2305], ct_source=[0.96, 0.77, 0.41]))
b3['v_ms'] = (b3.v_tunnel_raw_ms - 0.2).round(2)
b3['rpm'] = 424.5
b3['pitch_deg'] = PITCH
b3['rho'] = [1.19, 1.19, 1.20]
b3['tsr'] = (424.5 * math.pi / 30 * R / b3.v_ms).round(3)
write(os.path.join(EXP, 'mexico2006_balance_thrust.csv'), [
    'MEXICO 2006, the same three datapoints: rotor thrust estimated from the 6-component balance.',
    'Source: Bechmann, Sorensen & Zahle 2011, Table III (TABULATED): axial_force_balance_n = total axial force on the',
    'balance at the tower foot, tower_drag_est_n = their COARSE tower-drag estimate (single Cd from the parked',
    'run at 30 m/s, tower d ~0.5 m, h 5 m, Re near critical), thrust_n = difference; ct_source = their CT with',
    'rho = 1.225 and the raw speed. v_ms = raw - 0.2 m/s (tunnel-speed correction, see mexico2006_spanwise_axial.csv);',
    'rho column = MEXICO densities from ECN-E--14-048. Uncertainty: dominated by the tower drag estimate, of order',
    '+-100..300 N (the authors call it "a very rough estimate").',
    'Also known (not stored): ECN-E--07-042 Table 1 lists balance axial forces at pitch -1.2/-1.3 deg incl. tower',
    'drag (1039/1979/2487/3285/4091 N at 10/15/18/24/30 m/s).'],
    b3[['v_ms', 'v_tunnel_raw_ms', 'tsr', 'rpm', 'pitch_deg', 'rho', 'thrust_n', 'axial_force_balance_n',
        'tower_drag_est_n', 'ct_source']])
