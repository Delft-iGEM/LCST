"""Physics-informed LCST predictor for ELP sequences.

Combines:
  * sequence parsing into pentads + guest residues (`sequence`),
  * the Urry T_tc hydrophobicity scale with pH-dependent ionisation (`urry`),
  * the Chilkoti unified length/concentration equation (`chilkoti`).

The prediction is deterministic and needs no training data. An optional
machine-learning calibration layer (`calibrate`) can refine it against
experimental measurements.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict

from . import chilkoti
from .sequence import ParsedELP, parse_elp
from .urry import GuestComposition, guest_composition


@dataclass
class Prediction:
    tt_celsius: float            # predicted transition temperature (LCST), degC
    n_pentads: int               # ELP chain length used
    mean_guest_ttc: float        # mean Urry hydrophobicity of guests (degC)
    effective_fraction_ala: float
    fraction_charged: float
    fraction_aromatic: float
    concentration_uM: float
    pH: float
    molecular_weight_da: float
    is_elp: bool
    motif_match_fraction: float
    warnings: list[str]

    def as_dict(self) -> dict:
        return asdict(self)


def _build_warnings(parsed: ParsedELP, comp: GuestComposition) -> list[str]:
    w: list[str] = []
    if not parsed.is_elp:
        w.append(
            "Sequence does not look like a canonical (VPGXG)n ELP "
            f"(only {parsed.motif_match_fraction:.0%} of pentads match the "
            "...PG.G consensus); prediction is extrapolative."
        )
    if parsed.n_pentads < 5:
        w.append(
            f"Very short ELP ({parsed.n_pentads} pentads); the unified model "
            "is calibrated for longer chains and over-predicts short-chain Tt."
        )
    if comp.fraction_charged > 0.05:
        w.append(
            "Guest residues include ionisable side chains; Tt is strongly "
            "pH- and salt-dependent. Verify pH and ionic strength."
        )
    fa = chilkoti.effective_fa(comp.mean_ttc)
    if fa < -0.2 or fa > 1.2:
        w.append(
            "Guest hydrophobicity is outside the Val/Ala calibration range; "
            "concentration/length sensitivity is extrapolated."
        )
    return w


def predict_from_parsed(
    parsed: ParsedELP,
    concentration_uM: float = 25.0,
    pH: float = 7.4,
) -> Prediction:
    comp = guest_composition(parsed.guests, pH=pH)
    tt = chilkoti.transition_temperature(
        mean_ttc=comp.mean_ttc,
        length_pentads=parsed.n_pentads,
        concentration_uM=concentration_uM,
    )
    return Prediction(
        tt_celsius=round(tt, 1),
        n_pentads=parsed.n_pentads,
        mean_guest_ttc=round(comp.mean_ttc, 2),
        effective_fraction_ala=round(chilkoti.effective_fa(comp.mean_ttc), 3),
        fraction_charged=round(comp.fraction_charged, 3),
        fraction_aromatic=round(comp.fraction_aromatic, 3),
        concentration_uM=concentration_uM,
        pH=pH,
        molecular_weight_da=round(parsed.molecular_weight, 1),
        is_elp=parsed.is_elp,
        motif_match_fraction=round(parsed.motif_match_fraction, 3),
        warnings=_build_warnings(parsed, comp),
    )


def predict(
    sequence: str,
    concentration_uM: float = 25.0,
    pH: float = 7.4,
) -> Prediction:
    """Predict the LCST (transition temperature) of an ELP sequence.

    Parameters
    ----------
    sequence : one-letter amino-acid sequence of the ELP.
    concentration_uM : ELP concentration in micromolar (default 25 uM).
    pH : solution pH (default 7.4; affects ionisable guest residues).
    """
    parsed = parse_elp(sequence)
    return predict_from_parsed(parsed, concentration_uM=concentration_uM, pH=pH)
