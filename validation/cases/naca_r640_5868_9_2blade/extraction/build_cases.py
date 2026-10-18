"""Builds validation/cases/naca_r640_5868_9_{2,3,4}blade from the digitised NACA R-640 charts (Fig. 2 blade form,
Figs. 5/6/8 (2 blades), 9/11 (3 blades), 13/15 (4 blades)) and the XFOIL Clark Y family polars."""
import json, os, shutil, glob, numpy as np, pandas as pd
W = '/tmp/naca640/work'; XF = '/tmp/naca640/xf'
ROOT = '/Users/miha/Documents/dev/Python-BEM/validation/cases'
FT = 0.3048; D = 10 * FT; R = D / 2
NU = 1.4607e-5; RHO = 1.225; A_SOUND = 340.3
# ---------------- geometry (Fig. 2) ----------------
st = pd.read_csv(f'{W}/fig2_stations.csv')
uiuc = pd.read_csv(f'{W}/../nr640_geom.txt', sep=r'\s+')
FAM = [0.08, 0.09, 0.10, 0.12, 0.14, 0.18, 0.24]
rows = []
for _, s in st.iterrows():
    rR = s.r_R
    if rR > 0.95: continue
    beta15 = s.beta15
    rows.append(dict(r_R=rR, b_D=s['b/D'], h_b=s['h/b'], beta15=beta15))
# tip station (Fig. 2 curves end at r/R 0.95): UIUC NR640 design file, which reproduces Fig. 2 (see notes)
u1 = uiuc[uiuc['r/R'] == 1.0].iloc[0]
rows.append(dict(r_R=1.0, b_D=u1['c/R'] / 2, h_b=np.nan, beta15=u1['beta']))
G = pd.DataFrame(rows)
G['h_b'] = G['h_b'].where(G.r_R < 0.95, np.nan)
# thickness: Fig. 2 h/b from 0.30 to 0.947; 0.95 and 1.0 held at the 0.947 value (0.0825)
hb = G.h_b.copy(); hb[G.r_R >= 0.95] = 0.0812   # UIUC nr640_thick.txt value at r/R 0.95
def foil(rR, t):
    if rR < 0.26: return 'circle'
    if rR < 0.31: return 'transition'
    return 'clarky_t%02d' % round(100 * min(FAM, key=lambda f: abs(f - t)))
geo = pd.DataFrame({'r_m': (G.r_R * R).round(5), 'chord_m': (G.b_D * D).round(5),
                    'twist_deg': G.beta15.round(3), 'airfoil': [foil(r, t) for r, t in zip(G.r_R, hb)],
                    'thickness_ratio': [round(float(t), 4) if np.isfinite(t) else '' for t in hb]})
geo.loc[geo.airfoil == 'circle', 'thickness_ratio'] = ''
GEO_HDR = """# NACA R-640 (Hartman & Biermann 1938) propeller 5868-9 (Navy drawing), Clark Y sections, D = 10 ft (R = 1.524 m)
# Source: Fig. 2 'Blade-form curves for propellers 5868-9, 5868-R6 and 37-3647' (report p. 2), DIGITISED from the
#   NTRS scan (300 dpi) by ridge tracing (extraction/fig2_final.py, overlay extraction/fig2_overlay.png).
#   chord_m = (b/D) * D; twist_deg = beta = atan((P/D)/(pi r/R)) from the geometric-pitch curve 'P/D, set 15 deg at 0.75R'.
#   => twist_deg is the blade angle from the plane of rotation AT THE 15 deg SETTING (beta(0.75R) = 15.0 deg).
#   For a test at blade angle B at 0.75R the experiment files carry pitch_deg = B - 15.
#   The 25 and 35 deg P/D curves of Fig. 2 give the same beta - beta(0.75R) within 0.3 deg for r/R >= 0.3.
# Cross-check: UIUC Propeller Data Site vol. 2, NR640 design geometry (nr640_geom.txt, scaled replica of 5868-9 at 15 deg):
#   chord agrees within 1 %, beta within 0.1 deg, t/c within 0.004 at all common stations (extraction/check_geometry.csv).
# r/R = 1.00 row (chord, beta) is taken from that UIUC file (Fig. 2 curves stop at r/R ~0.95); tip shape rounded (Fig. 1).
# thickness_ratio = h/b (Fig. 2, r/R 0.30-0.947); 0.0812 at r/R 0.95-1.0 (UIUC nr640_thick.txt at 0.95). Root r/R 0.20-0.25: blade shank
#   (no h/b given, round/oval) -> 'circle'; r/R 0.30 (t/c 0.36) -> 'transition' circle -> clarky_t24.
# airfoil = nearest member of the generated Clark Y family (t/c 0.08 ... 0.24); see case.json airfoils.
# Digitising uncertainty: chord +-0.5 mm.. +-1.5 mm (+-0.0003 b/D), beta +-0.15 deg, h/b +-0.002.
"""
# ---------------- experiments ----------------
CHARTS = {2: ('f8_ct_2b', 'f6_cp_2b', 'Fig. 8 (C_T) and Fig. 6 (C_P), pp. 4-5', 0.30),
          3: ('f9_ct_3b', 'f11_cp_3b', 'Fig. 9 (C_T, p. 6) and Fig. 11 (C_P, p. 7)', 0.40),
          4: ('f13_ct_4b', 'f15_cp_4b', 'Fig. 13 (C_T, p. 8) and Fig. 15 (C_P, p. 9)', 0.45)}
RPM = {2: {15: 1200, 20: 1200, 25: 1000, 30: 1000, 35: 800, 40: 900, 45: 800},
       3: {15: 1000, 20: 1000, 25: 1000, 30: 900, 35: 800, 40: 800, 45: 700},
       4: {15: 1000, 20: 1000, 25: 900, 30: 800, 35: 700, 40: 700, 45: 600}}
def curve(name, lab):
    c = json.load(open(f'{W}/dig_{name}.json'))['curves'][str(lab)]
    x = np.array(c['x']); y = np.array(c['y']); o = np.argsort(x); x, y = x[o], y[o]
    x, i = np.unique(np.round(x, 4), return_index=True)
    return x, y[i]
def experiment(B):
    ctn, cpn, figs, jmin = CHARTS[B]
    out = []
    for ang in (15, 20, 25, 30, 35, 40, 45):
        xt, yt = curve(ctn, ang); xp, yp = curve(cpn, ang)
        jmax = min(xt.max(), xp.max())
        for J in np.round(np.arange(jmin, jmax + 1e-9, 0.05), 3):
            ct = float(np.interp(J, xt, yt)); cp = float(np.interp(J, xp, yp))
            if ct < 0.003 or cp < 0.003:   # near zero thrust eta is meaningless; the charts end at C_T = 0
                continue
            rpm = min(RPM[B][ang], 1000.0 / J)   # report p. 3: schedule below 115 mph, else rpm = 1000/(V/nD)
            n = rpm / 60; v = J * n * D
            out.append(dict(J=J, ct=round(ct, 5), cp=round(cp, 5), eta=round(J * ct / cp, 4), pitch_deg=ang - 15,
                            blade_angle_075R_deg=ang, rpm=round(rpm, 1), v_ms=round(v, 2),
                            tip_mach=round(np.hypot(np.pi * n * D, v) / A_SOUND, 3),
                            ct_err=0.0015, cp_err=0.002 if B == 2 else 0.003))
            e = out[-1]; e['eta_err'] = round(e['eta'] * np.hypot(e['ct_err'] / ct, e['cp_err'] / cp), 4)
    return pd.DataFrame(out)
# digitised eta (Fig. 5) for the 2-blade check
def eta_check(df):
    rows = []
    for ang in (15, 20, 25, 30, 35, 40, 45):
        xe, ye = curve('f5_eta_2b', ang)
        for _, r in df[df.blade_angle_075R_deg == ang].iterrows():
            if xe.min() <= r.J <= xe.max():
                rows.append(dict(angle=ang, J=r.J, eta_from_ct_cp=r.eta, eta_fig5=round(float(np.interp(r.J, xe, ye)), 4)))
    return pd.DataFrame(rows)
# ---------------- polars ----------------
def read_pol(f):
    L = open(f).read().splitlines()
    i = [k for k, l in enumerate(L) if l.strip().startswith('------')][-1]
    d = np.array([list(map(float, l.split()[:5])) for l in L[i + 1:] if len(l.split()) >= 5])
    return d
def polars(name):
    rows = []
    for f in sorted(glob.glob(f'{XF}/{name}_re*.pol')):
        re = float(f.split('_re')[-1][:-4]); d = read_pol(f)
        if len(d) == 0: continue
        d = d[np.argsort(d[:, 0])]; _, i = np.unique(d[:, 0], return_index=True); d = d[i]
        for a, cl, cd, cdp, cm in d:
            rows.append(dict(re=re, alpha_deg=a, cl=cl, cd=cd, cm=cm, ncrit=9, source='XFOIL 6.99 Ncrit=9'))
    return pd.DataFrame(rows)
# ---------------- write ----------------
for B in (2, 3, 4):
    cid = f'naca_r640_5868_9_{B}blade'; C = f'{ROOT}/{cid}'
    for sub in ('airfoils', 'experiment', 'extraction', 'sources'):
        os.makedirs(f'{C}/{sub}', exist_ok=True)
    with open(f'{C}/geometry.csv', 'w') as f:
        f.write(GEO_HDR); geo.to_csv(f, index=False)
    used = sorted(set(geo.airfoil) - {'circle', 'transition'} | {'clarky_t24'})
    airf = {}
    for n in used:
        shutil.copy(f'{XF}/{n}.dat', f'{C}/airfoils/{n}.dat')
        p = polars(n)
        with open(f'{C}/airfoils/{n}_polars.csv', 'w') as f:
            f.write(f'# {n}: Clark Y (UIUC coordinate database clarky.dat, t/c 11.7 %) with ALL y-ordinates scaled by t/c / 0.1171\n'
                    '#   (simple ordinate scaling, flat-bottom shape kept, camber scales with thickness) - the usual way propeller\n'
                    '#   Clark Y sections of other thickness were derived (cf. NACA TR 352 Clark Y thickness family).\n'
                    '# Polars COMPUTED with XFOIL 6.99 (repo binary), Ncrit = 9, free transition, Mach 0 (incompressible), 200 panels,\n'
                    '#   alpha sweep 0 -> -10 and 0 -> 20 deg in 0.5 deg steps; non-converged points are absent. Re = 0.5, 1, 2, 3 million\n'
                    '#   covers the full-scale blade (Re(0.75R) ~ 1.2e6 - 2.0e6, root ~0.6e6).\n')
            p.to_csv(f, index=False)
        airf[n] = {'coordinates': f'airfoils/{n}.dat', 'polars': f'airfoils/{n}_polars.csv',
                   'polar_source': 'XFOIL 6.99, Ncrit=9, Re 0.5e6-3e6, incompressible; coordinates = UIUC Clark Y y-scaled to t/c %.2f' % (int(n[-2:]) / 100)}
    airf['circle'] = {'coordinates': 'airfoils/circle.dat', 'polars': 'airfoils/circle_polars.csv',
                      'polar_source': 'blade shank r/R 0.20-0.25: cylinder, cl = 0, cd = 1.0 (subcritical circular cylinder, Re ~ 2e5-5e5)'}
    t = np.linspace(0, 2 * np.pi, 61)
    with open(f'{C}/airfoils/circle.dat', 'w') as f:
        f.write('circle (blade shank, generated, D = chord)\n')
        for a in t: f.write('%.6f %.6f\n' % (0.5 + 0.5 * np.cos(a), 0.5 * np.sin(a)))
    with open(f'{C}/airfoils/circle_polars.csv', 'w') as f:
        f.write('# circle: blade shank, cl = 0, constant cd = 1.0 (subcritical cylinder; the shank carries almost no load)\n')
        f.write('re,alpha_deg,cl,cd\n')
        for a in (-180, 0, 180): f.write(f'1000000,{a},0.0,1.0\n')
    df = experiment(B); ctn, cpn, figs, jmin = CHARTS[B]
    if B == 2:
        ec = eta_check(df); m = ec.set_index(['angle', 'J']).eta_fig5
        df['eta_fig5'] = [m.get((a, j), np.nan) for a, j in zip(df.blade_angle_075R_deg, df.J)]
    with open(f'{C}/experiment/performance.csv', 'w') as f:
        f.write(f"""# NACA R-640 (Hartman & Biermann 1938), propeller 5868-9, Clark Y, {B} blades, D = 10 ft, NACA 20-ft propeller-research tunnel
# Source: {figs}; curves 'traced from the original curves of faired test points' (report p. 3) - DIGITISED here
#   (300 dpi NTRS scan, ridge tracing calibrated on the printed grid, overlays extraction/ov_{ctn}.png, ov_{cpn}.png).
# Sampled every 0.05 in V/nD from J = {jmin} (the charts draw the low-J take-off part as DASHED, i.e. faired/extrapolated,
#   up to J ~0.25-0.41; it is left out) to the end of each curve (C_T -> 0). eta = J*ct/cp computed from the two charts.
# ct = T_e/(rho n^2 D^4) with T_e = shaft tension - change of nacelle drag due to slipstream (EFFECTIVE thrust, propeller
#   in front of a liquid-cooled engine nacelle, max. 43 x 38 in, no spinner); cp = P/(rho n^3 D^5) (shaft power).
# pitch_deg = blade angle at 0.75R - 15 deg (geometry.csv twist is the blade angle at the 15 deg setting).
# rpm: report p. 3 schedule (constant rpm while tunnel speed < 115 mph), at higher J rpm = 1000/(V/nD); v_ms = J n D.
#   tip_mach = helical tip Mach number (a = 340.3 m/s) - extra column, compressibility matters above ~0.6 (NACA TR 639).
# blade_angle_075R_deg, tip_mach: source-native/extra columns (ignored by the runner).
# eta_err = propagated from ct_err/cp_err (grows near zero thrust where eta = J ct/cp is ill-conditioned).
{'# eta_fig5 (2 blades only, extra column): eta read directly off Fig. 5; agrees with J*ct/cp within ~0.01 except near zero thrust.' if B == 2 else '#'}
# Uncertainty (digitising + scan): ct +-0.0015, cp +-{0.002 if B == 2 else 0.003}, J +-0.005; report: repeat tests within ~1 %.
""")
        df.to_csv(f, index=False)
    exps = [{'file': 'experiment/performance.csv', 'kind': 'performance', 'x': 'J', 'y': ['ct', 'cp', 'eta'],
             'defaults': {'rho': RHO, 'kin_viscosity': NU},
             'conditions': f'Full-scale 10-ft propeller ahead of a liquid-cooled engine nacelle (Curtiss Conqueror, geared 7:5), NACA 20-ft open-throat propeller-research tunnel, V up to ~115 mph (51 m/s); blade angles 15-45 deg at 0.75R (pitch_deg = angle - 15); rpm per row from the report schedule; thrust is EFFECTIVE thrust (shaft tension minus slipstream increment of nacelle drag); no tunnel-wall correction stated (open throat 20 ft, prop 10 ft).',
             'source': f'Hartman & Biermann 1938, NACA Report 640, {figs} (digitised)',
             'uncertainty': f'digitising ct +-0.0015, cp +-{0.002 if B == 2 else 0.003}, J +-0.005; repeat tests within ~1 % (report p. 3)'}]
    if B == 2:
        ec = eta_check(df)
        ec.to_csv(f'{C}/extraction/check_eta_fig5_vs_ct_cp.csv', index=False)
        m = ec.set_index(['angle', 'J']).eta_fig5
        df['eta_fig5'] = [m.get((a, j), np.nan) for a, j in zip(df.blade_angle_075R_deg, df.J)]
        print('eta check: mean |diff| %.4f, max %.4f (n=%d)' % ((ec.eta_from_ct_cp - ec.eta_fig5).abs().mean(),
              (ec.eta_from_ct_cp - ec.eta_fig5).abs().max(), len(ec)))
    case = {
        'id': cid, 'name': f'NACA 5868-9 propeller (Clark Y), {B} blades, 10 ft - NACA Report 640',
        'turbine_type': 1, 'B': B, 'R': round(R, 4), 'Rhub': round(0.2 * R, 4),
        'fluid': {'rho': RHO, 'kin_viscosity': NU, 'medium': 'air'},
        'geometry_file': 'geometry.csv', 'airfoils': airf, 'experiments': exps,
        'sources': [
            {'key': 'hartman1938', 'citation': 'Hartman, E.P., Biermann, D.: The aerodynamic characteristics of full-scale propellers having 2, 3, and 4 blades of Clark Y and R.A.F. 6 airfoil sections. NACA Report No. 640, 1938.',
             'url': 'https://ntrs.nasa.gov/citations/19930091715', 'ntrs': '19930091715'},
            {'key': 'uiuc_nr640', 'citation': 'Deters, R.W., Ananda, G.K., Selig, M.S.: Reynolds number effects on the performance of small-scale propellers. AIAA 2014-2151 (NR640 = scaled 5868-9), and UIUC Propeller Data Site vol. 2, NR640 design geometry files nr640_geom.txt / nr640_thick.txt',
             'doi': '10.2514/6.2014-2151', 'url': 'https://m-selig.ae.illinois.edu/props/volume-2/propDB-volume-2.html'},
            {'key': 'uiuc_clarky', 'citation': 'UIUC Airfoil Coordinates Database, clarky.dat', 'url': 'https://m-selig.ae.illinois.edu/ads/coord/clarky.dat'},
            {'key': 'biermann1938_tr639', 'citation': 'Biermann, D., Hartman, E.P.: The effect of compressibility on eight full-scale propellers operating in the take-off and climbing range. NACA Report 639, 1938 (tip-speed effects on these data).', 'url': 'https://ntrs.nasa.gov/citations/19930091714'}],
        'notes': ('Geometry digitised from R-640 Fig. 2 and cross-checked against the UIUC NR640 design table (agreement within 1 % chord, 0.1 deg beta). '
                  'twist_deg = blade angle beta from the plane of rotation at the 15 deg setting (beta(0.75R) = 15 deg exactly); pitch_deg in the experiment = (blade angle at 0.75R) - 15, so a 30 deg test has pitch_deg = 15. '
                  'Rhub = 0.2R = first tabulated station (blade shank, modelled as a cylinder with cd 1.0 to r/R 0.25, transition at 0.30); the propeller runs in front of a large nacelle (43 x 38 in, i.e. ~0.30-0.36 of D), no spinner, so the inner ~0.3R sits in the nacelle blockage - the measured thrust is EFFECTIVE thrust (minus slipstream drag increment on the nacelle) and eta is propulsive efficiency of the propeller-nacelle unit, both a little below isolated-propeller values. '
                  'Airfoils: Clark Y family generated by y-scaling the UIUC Clark Y; polars are XFOIL (Ncrit 9, incompressible), no measured section data used (measured Clark Y: NACA TR 502 Silverstein, 11.7 % only, Re up to 8.6e6). '
                  'Reynolds number at 0.75R ~1.2e6-2.0e6. COMPRESSIBILITY: tip speeds 600-1200 rpm -> helical tip Mach 0.29-0.6 (column tip_mach); the 2-blade 15/20 deg runs at 1200 rpm reach M_tip ~0.57-0.6 - per NACA TR 639 these coefficients are slightly lowered in the take-off range; the polars are incompressible, so add a Prandtl-Glauert correction if the solver has one. '
                  'Low-J take-off parts of the charts are drawn dashed (faired/extrapolated) and are excluded. The report states the tabulated test points were only available on request (not published), so the faired chart curves are the data.')}
    json.dump(case, open(f'{C}/case.json', 'w'), indent=2)
    # extraction audit trail + sources
    for f in glob.glob(f'{W}/*.py') + glob.glob(f'{W}/dig_*.json') + glob.glob(f'{W}/ov_f*.png') + [f'{W}/fig2_overlay.png', f'{W}/fig2_stations.csv', f'{W}/dash_ends.json']:
        if os.path.exists(f): shutil.copy(f, f'{C}/extraction/')
    shutil.copy('/tmp/naca640/build/build_cases.py', f'{C}/extraction/')
    shutil.copy(f'{XF}/make_foils_and_run.py', f'{C}/extraction/')
    for f in ['r640.pdf', 'deters2014.pdf', 'nr640_geom.txt', 'nr640_thick.txt']:
        shutil.copy(f'/tmp/naca640/{f}', f'{C}/sources/')
    shutil.copy(f'{XF}/clarky_uiuc.dat', f'{C}/sources/')
    print(cid, len(df), 'points', df.groupby('blade_angle_075R_deg').J.agg(['min', 'max', 'count']).to_dict('index'))
# geometry cross-check table
u = pd.read_csv(f'{W}/../nr640_geom.txt', sep=r'\s+'); th = pd.read_csv(f'{W}/../nr640_thick.txt', sep=r'\s+')
ck = pd.DataFrame({'r_R': G.r_R, 'c_R_fig2': 2 * G.b_D, 'beta15_fig2': G.beta15, 'h_b_fig2': G.h_b})
ck = ck.merge(u.rename(columns={'r/R': 'r_R', 'c/R': 'c_R_uiuc', 'beta': 'beta_uiuc'}), on='r_R', how='left')
ck = ck.merge(th.rename(columns={'r/R': 'r_R', 't/c': 'tc_uiuc'}), on='r_R', how='left')
for B in (2, 3, 4): ck.round(4).to_csv(f'{ROOT}/naca_r640_5868_9_{B}blade/extraction/check_geometry.csv', index=False)
print(ck.round(4).to_string())
