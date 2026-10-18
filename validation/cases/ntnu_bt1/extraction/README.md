# ntnu_bt1 — audit trail

Everything in the case is rebuilt by `python3 build_case.py` from `../sources/*.txt` (pdftotext -layout of the PDFs)
and the `*_raw.csv` files the extraction scripts below write. `python3 check_case.py` runs the sanity checks and
draws `check_ntnu_bt1.png`.

| script | source figure | method | overlay |
|---|---|---|---|
| `extract_krogstad_eriksen2013_fig3.py` | Krogstad & Eriksen 2013, Fig. 3 (BT1 measured CP, CT) | vector PDF, exact marker centres (PyMuPDF) | `overlay_krogstad_eriksen2013_fig3.png` |
| `extract_krogstad_adaramola2012_fig4.py` | Krogstad & Adaramola 2012, Fig. 4 (Re dependence, 7 speeds) | vector PDF, series by marker shape | `overlay_krogstad_adaramola2012_fig4.png` |
| `extract_krogstad_adaramola2012_fig5.py` | Krogstad & Adaramola 2012, Fig. 5 (yaw 0–50°) | vector PDF, series by marker style/colour | `overlay_krogstad_adaramola2012_fig5.png` |
| `digitise_bartl2018b_fig1.py` | Bartl et al. 2018 (WES 3:489), Fig. 1 (yaw sweep, TSR 6) | raster image, colour masks, grid-line calibration | `overlay_bartl2018b_fig1.png` |
| `digitise_bartl2019_fig4.py` | Bartl et al. 2019 (EJMB 75:180), Fig. 4 (NTNU S826 polars) | raster image, exact legend colours, marker-width centre | `overlay_bartl2019_fig4.png` |

Tabulated (copied by `build_case.py`, no digitising): geometry and S826 coordinates and XFOIL polar (BT5 invitation
Tables C.1, B.1, D.1), Krogstad & Lund 2012 Tables I (cross-check of the geometry) and II (CP, CT), Sarlak 2018
Tables A2/A3 (DTU polars), Bartl 2018a Table 2.

`pages/` holds the page renders / embedded figure images used.

Findings while extracting:
- pdftotext renders the minus sign of the Krogstad & Lund PDF (MathTime MT2SYS font) as `\x02`; the last twist
  value is −0.71674° and the last CP −0.013258.
- Krogstad & Adaramola 2012 Fig. 5: the CT panel's markers are at 1.024× the TSR of the CP/Cq panels for all six
  series — the CP-panel TSR is used (see the CSV header). The caption names panels (b)/(c) the wrong way round.
- Krogstad & Adaramola 2012 Fig. 4: its CT is ~0.13 above the rotor CT of all other sources → stored as
  `ct_incl_tower`, not `ct`.
