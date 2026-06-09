# elp-lcst — LCST prediction for elastin-like polypeptides

Predicts the **LCST / inverse transition temperature (T_t)** of an elastin-like
polypeptide (ELP) directly from its amino-acid sequence, accounting for guest
composition, chain length, concentration and pH.

ELPs are built from the pentapeptide repeat **VPGXG** (Val-Pro-Gly-**X**aa-Gly).
On heating past T_t they undergo a sharp, reversible coacervation (the LCST
transition). T_t is governed by three things, all captured here:

1. **Guest-residue hydrophobicity** — hydrophobic X (Phe, Leu, Val) lowers T_t;
   polar/charged X (Ser, Glu⁻) raises it.
2. **Chain length** — longer ELPs transition at lower T_t.
3. **Concentration** — T_t falls logarithmically with concentration.

## The model

The predictor is **physics-informed** and needs no training data:

```
T_t = T_tc + (k / L) · ln(C_c / C)
```

* `T_tc` — the guest composition's mean **Urry characteristic transition
  temperature** (the Urry 1992 T_t-based hydrophobicity scale), with
  pH-dependent ionisation of D/E/H/K/C/Y side chains via Henderson–Hasselbalch.
* `k`, `C_c` — concentration-sensitivity and critical-concentration parameters
  from the **Chilkoti-group unified model** (McDaniel, Radford & Chilkoti,
  *Biomacromolecules* 2013, 14, 2866), fitted as smooth functions of guest
  hydrophobicity from their Table 1.
* `L` — chain length in pentapeptides; `C` — concentration (µM).

For Val/Ala hydrophobic ELPs the engine **is** the published unified model
(validated to r² > 0.99 against 120 measurements). For other guests the Urry
scale shifts `T_tc` correctly while the length/concentration physics is carried
over.

### Accuracy

| Regime | Expected error | Notes |
|---|---|---|
| Val/Ala ELPs | ~2–3 °C | engine reproduces the published unified model |
| poly(VPGVG) | ~3 °C | predicts 32 °C for a 120-mer at 25 µM (lit. ~25–32 °C) |
| Val/Ala/Gly (e.g. V5A2G3) | up to ~10 °C, over-predicts | Gly raises mean hydrophobicity but lowers T_t more than an equally-hydrophobic Ala mix; see below |
| charged / very hydrophobic guests | extrapolative | correct direction; magnitude uncertain — flagged in output |

Composition-trend interpolation of the Chilkoti Table 1 (leave-one-out): R² =
0.99 (T_tc), 0.95 (k), 0.95 (C_c). Independent check on **ELP[V5A2G3-90]**:
measured 50.0 °C, physics 59.5 °C (+9.5 °C).

**The systematic offset for Gly-rich / charged ELPs is removed by the
calibration layer.** Training the residual correction on the V5A2G3-90 point
(plus the unified V/A grid) brings that prediction to 50.0 °C while leaving the
V/A grid intact (cross-validated MAE 0.02 °C). For accurate predictions on a
specific ELP family, calibrate on a handful of your own measurements (below).

### Sources

* Urry, *Hydrophobicity scale for proteins based on inverse temperature
  transitions*, Biopolymers 1992, 32, 1243; J. Phys. Chem. B 1997, 101, 11007.
* McDaniel, Radford, Chilkoti, *A unified model for de novo design of ELPs with
  tunable inverse transition temperatures*, Biomacromolecules 2013, 14, 2866
  ([PMC3779073](https://pmc.ncbi.nlm.nih.gov/articles/PMC3779073/)).
* Meyer & Chilkoti, *Quantification of the effects of chain length and
  concentration…*, Biomacromolecules 2004, 5, 846.
* Christensen *et al.*, *Predicting transition temperatures of ELP fusion
  proteins*, Biomacromolecules 2013
  ([PMC3667497](https://pmc.ncbi.nlm.nih.gov/articles/PMC3667497/)).

## Install

Uses [uv](https://docs.astral.sh/uv/) for environment management.

```bash
uv sync --extra dev      # create venv + install
```

## Usage

### Command line

```bash
# From a raw sequence
uv run elp-lcst predict "$(python -c "print('VPGVG'*120)")" -c 25

# From the compact Chilkoti guest spec (V:A:G = 5:2:3, 9 repeats = 90 pentads)
uv run elp-lcst predict --spec V5A2G3 --repeats 9 -c 25 --pH 7.4

# JSON output
uv run elp-lcst predict --spec V1 --repeats 120 --json

# Accuracy report (fit quality + independent literature check)
uv run elp-lcst validate
```

### Python

```python
from elp_lcst import predict

p = predict("VPGVG" * 120, concentration_uM=25, pH=7.4)
print(p.tt_celsius)        # predicted LCST in °C
print(p.as_dict())         # full breakdown + warnings
```

## Calibrating on your own data

If you have measured T_t values (e.g. from the lab's LCST-extraction tool),
fit a residual correction for best accuracy on your constructs.

CSV with sequences:

```csv
sequence,concentration_uM,pH,tt_celsius
VPGVGVPGVG...,25,7.4,31.5
...
```

```bash
uv run elp-lcst calibrate my_data.csv -o calibrated_model.pkl --with-unified
```

```python
from elp_lcst.calibrate import CalibratedModel
cm = CalibratedModel.load("calibrated_model.pkl")
cm.predict_tt("VPGVG"*120, concentration_uM=25)
```

The calibration learns `T_t,measured − T_t,physics = f(features)`, so it removes
systematic bias where you have data while the physics model still governs
extrapolation. Cross-validated MAE/RMSE are reported when you train.

## Scope & caveats

* Designed for **canonical (VPGXG)ₙ ELPs**. Fusion proteins need the separate
  surface-index model (Christensen/Chilkoti 2013) — not implemented here.
* Assumes physiological-strength buffer (PBS). Ionic strength (Hofmeister
  salting-out) shifts T_t and is not modelled explicitly; calibrate for other
  buffers.
* Proline at the guest position is non-physical (it breaks the β-turn) and is
  flagged.

## Tests

```bash
uv run pytest -q
```
