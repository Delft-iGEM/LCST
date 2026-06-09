"""Chilkoti-group unified model for ELP transition temperature.

McDaniel, Radford & Chilkoti, "A unified model for de novo design of
elastin-like polypeptides with tunable inverse transition temperatures",
Biomacromolecules 2013, 14, 2866-2872.

For a fixed guest composition the transition temperature obeys

    T_t = T_tc + (k / L) * ln(C_c / C)                              (Eq. 1)

where
    T_tc = composition characteristic transition temperature (degC)
    k    = concentration-sensitivity constant (degC * pentapeptides)
    C_c  = composition critical concentration (uM)
    L    = chain length in number of pentapeptides
    C    = ELP concentration (uM)

The paper fits T_tc, k and C_c for a Val/Ala guest library, indexed by the
alanine fraction f_A. Those fitted values (Table 1, PBS) are tabulated below.
We fit each parameter to the exponential form A*exp(b*f_A) used in the paper so
the model is continuous in composition, and we map an arbitrary guest
composition onto an effective f_A through its mean Urry T_tc.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import PchipInterpolator

# Table 1 of McDaniel et al. 2013 (PBS). f_A = alanine fraction of the guest
# residues in a Val/Ala ELP library. C_c is non-monotone (the f_A=0.9 row sits
# above f_A=0.8), so we interpolate rather than force a monotone functional fit.
_F_A = np.array([0.0, 0.2, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
_TTC = np.array([20.2, 22.9, 27.7, 29.9, 33.8, 39.3, 44.7, 54.3])      # degC
_K = np.array([144.1, 201.2, 316.9, 387.0, 433.7, 484.7, 507.0, 623.0])  # degC*pent
_CC = np.array([22130.2, 13573.3, 6218.9, 3684.3, 2748.4, 2005.3, 2412.1, 843.9])  # uM

# Shape-preserving (monotone-where-data-are) interpolation through every
# calibrated composition; reproduces Table 1 exactly at the knots.
_TTC_I = PchipInterpolator(_F_A, _TTC, extrapolate=True)
_K_I = PchipInterpolator(_F_A, _K, extrapolate=True)
_CC_I = PchipInterpolator(_F_A, _CC, extrapolate=True)

# Endpoints of the Urry mean-T_tc axis that anchor the f_A mapping:
# pure Val guest -> f_A = 0, pure Ala guest -> f_A = 1.
TTC_VAL = 24.0
TTC_ALA = 45.0


def _clamp01(fa: float) -> float:
    return min(1.0, max(0.0, fa))


def ttc_of_fa(fa: float) -> float:
    return float(_TTC_I(_clamp01(fa)))


def k_of_fa(fa: float) -> float:
    return float(_K_I(_clamp01(fa)))


def cc_of_fa(fa: float) -> float:
    return float(_CC_I(_clamp01(fa)))


def effective_fa(mean_ttc: float) -> float:
    """Map a composition's mean Urry T_tc onto the Val/Ala f_A axis.

    f_A = 0 at the Val value (24 degC), f_A = 1 at the Ala value (45 degC).
    The result is clamped to [0, 1] only for the k / C_c sensitivity lookup;
    the composition's own mean T_tc is used directly as the base temperature,
    so highly hydrophobic or charged guests outside the Val/Ala range still
    shift the prediction correctly.
    """
    return (mean_ttc - TTC_VAL) / (TTC_ALA - TTC_VAL)


@dataclass
class ChilkotiParams:
    ttc: float   # base characteristic transition temperature (degC)
    k: float     # concentration sensitivity (degC * pentapeptides)
    cc: float    # critical concentration (uM)
    fa_eff: float


def params_for_composition(mean_ttc: float) -> ChilkotiParams:
    """Return (T_tc, k, C_c) for a guest composition of the given mean Urry T_tc.

    The base T_tc is the composition's own mean Urry value. The k and C_c
    concentration-sensitivity parameters are read from the Val/Ala fits at the
    effective f_A (clamped to the calibrated [0, 1] range).
    """
    fa = effective_fa(mean_ttc)
    fa_clamped = min(1.0, max(0.0, fa))
    return ChilkotiParams(
        ttc=mean_ttc,
        k=k_of_fa(fa_clamped),
        cc=cc_of_fa(fa_clamped),
        fa_eff=fa,
    )


def transition_temperature(
    mean_ttc: float, length_pentads: int, concentration_uM: float
) -> float:
    """Predict T_t (degC) via the unified equation T_t = T_tc + (k/L) ln(C_c/C)."""
    if length_pentads < 1:
        raise ValueError("Chain length must be >= 1 pentapeptide.")
    if concentration_uM <= 0:
        raise ValueError("Concentration must be positive.")
    p = params_for_composition(mean_ttc)
    return p.ttc + (p.k / length_pentads) * np.log(p.cc / concentration_uM)


def fit_quality() -> dict[str, float]:
    """Leave-one-out interpolation R^2 for the composition functions.

    The interpolators pass through Table 1 exactly, so a plain R^2 would be 1
    trivially. Instead we drop each interior knot, re-interpolate from the rest,
    and score the held-out predictions - an honest estimate of how well the
    composition trend is captured between measured points.
    """
    out = {}
    for name, y in (("Ttc", _TTC), ("k", _K), ("Cc", _CC)):
        preds, trues = [], []
        for i in range(1, len(_F_A) - 1):  # interior knots only
            mask = np.arange(len(_F_A)) != i
            interp = PchipInterpolator(_F_A[mask], y[mask], extrapolate=True)
            preds.append(float(interp(_F_A[i])))
            trues.append(float(y[i]))
        preds, trues = np.array(preds), np.array(trues)
        ss_res = float(np.sum((trues - preds) ** 2))
        ss_tot = float(np.sum((trues - trues.mean()) ** 2))
        out[name] = 1.0 - ss_res / ss_tot
    return out
