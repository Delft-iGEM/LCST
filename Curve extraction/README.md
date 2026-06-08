# LCST Extractor

A browser-based tool for standardised, reproducible extraction of Lower Critical Solution Temperature (LCST) from absorbance-versus-temperature graphs published as images. No installation, no server, no dependencies to manage — open `tm_extractor.html` in any modern browser and start working.

---

## The problem it solves

Absorbance LCST curves are published as figures in papers, patents, theses, and lab reports. Most times the raw data is unavailable and in some cases the LCST isn't specified. Extracting this transition temperature with a consistent and reproducible protocol is of upmost importance for the construction of a dataset that furthers the analytical and quantitative study of ELPs. Extracting LCST values by eye is subjective, hard to document, and not reproducible between analysts. This tool replaces eyeballing with a fully documented, algorithmically consistent workflow:

1. Calibrate the axes against known tick marks on the image
2. Click points along the curve
3. The tool fits a LOESS smoother, computes the derivative dA/dT numerically, and reports LCST as the temperature of the peak gradient (in accordance with the general ELP corpus) — the inflection point of the sigmoid

Every parameter, every clicked coordinate, and every computed value is saved and can be reloaded, verified, or handed to a collaborator. This tool isn't exclusive to LCST or absorbance as a function of temperature — it can be used to extract datapoints from any graphed curve.

---

## Features

### Image handling
- Load images via file picker, drag-and-drop, or paste (Ctrl+V)
- Rotate 90° clockwise for graphs in non-standard orientation
- Works with any image size or plot style

### Axis calibration
- **Multi-point regression** — click as many tick marks as you want per axis; a least-squares linear regression is fitted through all of them, averaging out individual click errors
- **Centroid snapping** — when hovering near a tick mark, the cursor snaps to the intensity-weighted centroid of the surrounding dark pixels, removing sub-pixel ambiguity from thick or blurry tick lines; snap radius is adjustable (2–20 px)
- **Log or linear scale** — each axis is independently configurable; log-scale axes are fitted in log₁₀ space and back-transformed, giving correct calibration for concentration, wavelength, or any multiplicative axis
- **Live R²** — the goodness-of-fit of the calibration regression is shown in real time; values close to 1.00000 confirm consistent clicking
- **Reconstructed axis overlay** — once ≥ 2 calibration points are set, a dashed line and tick marks are drawn over the image at the positions predicted by the regression, letting you visually confirm alignment before digitising

### Curve digitisation
- **Multiple curves per image** — add as many curves as needed; each has its own colour, name, and free-text sample/condition label
- Click points along a curve in any order; they are sorted automatically before fitting
- Undo last point or clear all points for the active curve independently

### LCST computation
- LOESS smoothing with adjustable bandwidth k (1–40); larger k smooths more, smaller k follows the raw data more closely
- First derivative dA/dT computed by central finite differences on the smoothed curve
- LCST defined as the temperature at which |dA/dT| is maximised
- Results update live as points are added or smoothing is adjusted

### Visualisation
- Live chart panel showing all smoothed curves and their derivatives side by side
- LCST markers on derivative plots
- Data table of smoothed values for all curves

### Session save and reload
- **Export** produces a `.zip` containing:
  - The original image
  - `session.json` — full machine-readable record (see below)
  - One `curve_<name>.csv` per curve with raw, smoothed, and derivative columns
  - `summary.csv` — one row per curve with LCST, max dA/dT, n points, smoothing k
  - `METHOD.txt` — a plain-English methods description ready to paste into a paper
- **Load session** restores a previously exported `.zip` or `session.json` exactly, including all calibration points, scale settings, curve points, labels, and smoothing parameter

---

## Workflow

```
1. Load image
2. Calibration panel → set X unit and scale (linear/log)
3. Mode: Calibrate X → type a tick value, click the tick mark → repeat for ≥ 3 ticks
4. Mode: Calibrate Y → same for the Y axis
5. Check the reconstructed axis overlay aligns with the image axes
6. Add a curve, set its name and sample label
7. Mode: Digitize → click ~20–40 points along the curve
8. Read LCST from the Results panel
9. Repeat steps 6–8 for additional curves
10. Export session
```

---

## Calibration design

### What is stored and what is purely visual

The X axis calibration maps **horizontal pixel position → temperature** using only `pixel_x` from each calibration click. The vertical position (`pixel_y`) of where you clicked on an X tick mark is stored in the export for traceability but plays no role in any calculation. The reconstructed axis overlay line is drawn at the mean `pixel_y` of your X calibration points — this is visual only and does not affect LCST.

Symmetrically, the Y axis calibration maps **vertical pixel position → absorbance** using only `pixel_y`.

This means the calibration is genuinely 1-dimensional per axis, as it should be for a standard rectangular plot.

### Regression model

For each axis the tool fits:

```
transformed_value = slope × pixel_coordinate + intercept
```

where `transformed_value` is the raw tick value for linear axes, or `log₁₀(tick value)` for log axes. After fitting, pixel coordinates are converted to data values by applying the inverse transform.

---

## Output format

### `session.json`

```json
{
  "format_version": "2.0",
  "software": "LCST Extractor v2",
  "exported": "2025-06-01T14:23:00.000Z",
  "method": {
    "smoothing_algorithm": "LOESS",
    "smoothing_bandwidth_k": 10,
    "tm_definition": "Temperature at which |dA/dT| is maximised — recorded as LCST (lower critical solution temperature)",
    "derivative_method": "Central finite differences on LOESS-smoothed data",
    "calibration": {
      "x_axis": {
        "unit": "°C",
        "scale": "linear",
        "n_points": 5,
        "regression": { "slope": -0.5432, "intercept": 312.1, "r2": 0.999987, "fit_space": "linear" },
        "pixel_points": [
          { "pixel_x": 102.3, "pixel_y": 487.1, "real_value": 20 },
          "..."
        ]
      },
      "y_axis": { "..." }
    }
  },
  "image": { "label": "Fig2A", "notes": "WT protein pH 7.4", "original_width_px": 1200, "original_height_px": 900 },
  "curves": [
    {
      "name": "Curve 1",
      "label": "WT 150 mM NaCl",
      "tm_result": {
        "LCST": 58.432,
        "LCST_unit": "°C",
        "max_dA_dT": 0.02341
      },
      "n_points": 34,
      "raw_points": [ { "pixel_x": 203.1, "pixel_y": 341.2, "°C": 42.1, "Absorbance": 0.2341 }, "..." ],
      "smoothed_curve": [ "..." ],
      "derivative": [ { "T": 42.1, "dA_dT": 0.00394 }, "..." ]
    }
  ]
}
```

### `curve_<name>.csv`

```
# LCST Extractor — Curve 1
# Label: WT 150 mM NaCl
# LCST = 58.4321 °C
# Smoothing bandwidth k = 10
# Method: LOESS smooth + central finite difference derivative
T_°C,Absorbance_raw,Absorbance_smooth,dAbsorbance_dT
20.0000,0.18210,0.18341,
22.5000,0.19103,0.19287,0.000394
...
# LCST,58.4321
# max_dAbsorbance_dT,0.023410
```

### `summary.csv`

```
# LCST Extractor — Summary
Curve,Label,n_points,LCST_°C,max_dA_dT,smoothing_k
Curve 1,WT 150 mM NaCl,34,58.4321,0.023410,10
```

---

## Methods text

The `METHOD.txt` file in each export is written to be pasted directly into a paper's methods section. It documents the digitisation procedure, smoothing algorithm, derivative computation, LCST definition, calibration coordinates, R² values, axis scale types, and all numerical results.

---

## Technical notes

### Coordinate system

All stored pixel coordinates are in **original image pixels**, independent of zoom level or display scale. Zoom and pan affect only the on-screen rendering; stored values are invariant to them.

### Centroid snapping

When snapping is enabled, the tool searches a circular region of radius r (in image pixels) around the clicked point. It finds the intensity-weighted centroid of pixels below a luminance threshold defined as `min_luminance + 0.3 × (max_luminance − min_luminance)` within the search window. This places the calibration point at the optical centre of the tick mark, not at the cursor tip.

### Derivative and LCST

```
dA/dT(i) = [A_smooth(i+1) − A_smooth(i-1)] / [T(i+1) − T(i-1)]
```

LCST is the temperature at index i where |dA/dT(i)| is maximum. The LOESS smoother is a locally weighted linear regression with a tricube kernel; bandwidth k controls the number of neighbouring points used in each local fit.

---

## Dependencies

The tool is entirely self-contained in a single HTML file. Two JavaScript libraries are loaded from CDN at runtime:

| Library | Version | Purpose |
|---|---|---|
| Chart.js | 4.4.1 | Result charts |
| JSZip | 3.10.1 | Session export/import |

No build step, no package manager, no backend.

---

## Browser compatibility

Any modern browser with Canvas 2D API and File API support. Tested in Chrome, Firefox, and Safari. The file input, drag-and-drop, and Ctrl+V paste are all independent image loading paths.

---

## Reproducibility guarantee

Two analysts working from the same `.zip` session file with the same smoothing parameter k will obtain an identical LCST, because the computation is fully deterministic from the stored pixel coordinates and calibration regression. The only source of inter-analyst variability is the initial choice of which points to click along the curve, which is documented in the exported pixel coordinates.

---

## Limitations

- Calibration assumes the graph axes are linear or log-scale and orthogonal (standard rectangular plot). Polar plots, ternary diagrams, or non-linear axis transformations other than log are not supported.
- LCST extraction is meaningful for sigmoidal absorbance transitions. Multi-state transitions with multiple inflection points will return only the single largest gradient peak; inspect the derivative chart to confirm the result is the peak you intended.
- The tool reads images; it does not read embedded data from PDF vector graphics or Excel charts. If the original data file is available, extract it directly rather than digitising.

---

## License

MIT
****
