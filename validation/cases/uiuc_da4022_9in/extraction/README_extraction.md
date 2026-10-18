# Extraction trail — UIUC propeller cases (uiuc_*)

Built 2026-10-18. Nothing here is digitised: every number comes from a tabulated file.

## Raw downloads (copies of the files used are in `raw/`)
- UIUC Propeller Database, whole archive `https://m-selig.ae.illinois.edu/props/download/UIUC-propDB.zip`
  (md5 a41e484f1fd0fb6ff80b76e27410808b, matches the site). Vol. 1 = Brandt (APC props), Vol. 2 = Deters (DA props).
  Files `<prop>_geom.txt` (r/R, c/R, beta), `<prop>_<run>_<rpm>.txt` (J, CT, CP, eta), `<prop>_static_<run>.txt`
  (RPM, CT, CP). Geometry plots (`*_geom.png`) give the true diameter used for the coefficients.
- APC geometry files `PE0-FILES_WEB-202602.zipx` from https://www.apcprop.com/technical-information/file-downloads/
  (`*-PERF.PE0`, version v2025-1001): chord, twist, thickness ratio per station, airfoil sections.
- Airfoils: SDA1045/1075/1100 coordinates from Deters, Ananda & Selig 2014 (AIAA 2014-2151) Appendix A Tables 3–5
  (pdftotext of the paper, parsed by `build_uiuc_cases.py`); E63 and Clark-Y from the UIUC coordinate database
  (`e63.dat`, `clarky.dat`); NACA 4412 from the analytic 4-digit formulas. All normalised (LE at origin, LE–TE chord on
  the x axis, unit chord) before XFOIL.
- Clark-Y measured polar: UIUC LSAT Vol. 3 `volume03.zip`, `DRAG03.TXT`, block `CLARKYB.DRG` (`raw/lsat_vol3_DRAG03_clarkyb.txt`).
- OSU APC 10x7E data: `CT_OSU_UIUC.fig`, `CP_OSU_UIUC.fig` (MATLAB figures from the UIUC site); XData/YData of the
  series `10x7_{OSU}` read with scipy.io.loadmat — exact plotted values.

## Scripts
- `build_uiuc_cases.py` — builds all ten case folders (geometry, airfoils, XFOIL polars, experiments, case.json).
  `python3 build_uiuc_cases.py --raw /tmp/uiuc` (raw layout: `UIUC-propDB/`, `apc/PE0-FILES_WEB/`, `lsat/v03/`,
  `foils/`). XFOIL: repo binary `xfoil_executables/xfoil` (6.97), Ncrit 9 and 5, Re 20k–200k, results cached.
- `check_uiuc_cases.py` — sanity checks (geometry ends at R, eta = J·ct/cp, polars span 0°) and `check_<case>.png`.
- `clarky_xfoil_vs_lsat.py` (uiuc_apcsp_11x7 only) — XFOIL vs measured Clark-Y, `clarky_xfoil_vs_lsat.png`.

## Conversions
- r_m = r/R·R, chord_m = c/R·R (DA) or inch·0.0254 (APC); twist_deg = beta (UIUC) or APC TWIST (LE–TE datum), both
  the blade angle from the plane of rotation = solver convention; pitch_deg = hub offset of the DA "pitched" configurations.
- v_ms = J·n·D with the run's nominal rpm (from the file name) and the true diameter.
- DA tip: the PropellerScanner chord/twist collapse at r/R 0.95–1.0 is an image-edge artefact (the printed blades have
  a square or designed tip); replaced as stated in the geometry.csv header. Raw scan kept in `geometry_scan_raw.csv`.
