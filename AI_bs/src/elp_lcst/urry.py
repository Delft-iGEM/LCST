"""Urry T_t-based hydrophobicity scale for elastin-like polypeptides.

Urry and co-workers synthesised poly[f_v(VPGVG), f_x(VPGXG)] copolymers and
measured the inverse-transition temperature (T_t) as a function of the guest
mole fraction f_x. Extrapolating to f_x = 1 yields a characteristic transition
temperature (T_tc, in degC) for each guest residue X. A more hydrophobic guest
residue gives a *lower* T_tc; a polar / charged residue gives a *higher* one.

Reference values are the canonical Urry 1992 scale (Biopolymers 32, 1243;
J. Phys. Chem. B 1997, 101, 11007), as reproduced in the Chilkoti-group
reviews. Internal consistency check: the Val/Ala/Gly weighted average for the
classic ELP[V5A2G3] guest composition (0.5*V + 0.2*A + 0.3*G) is
0.5*24 + 0.2*45 + 0.3*55 = 37.5 degC, matching the 37.7 degC quoted by
Chilkoti et al. (Biomacromolecules 2013, 14, 2866) for that backbone.

For ionisable residues two values are stored - the protonated (neutral) and the
ionised (charged) form. The effective T_tc at a given pH is a
Henderson-Hasselbalch blend of the two, because partially ionised side chains
contribute an intermediate hydrophobicity.
"""

from __future__ import annotations

from dataclasses import dataclass

# Neutral-state characteristic transition temperature (degC). These are the
# values for the guest residue in its protonated / uncharged form.
URRY_TTC_NEUTRAL: dict[str, float] = {
    "W": -90.0,
    "Y": -55.0,
    "F": -30.0,
    "H": -10.0,   # neutral imidazole
    "P": -8.0,    # rarely used as a guest residue; included for completeness
    "L": 5.0,
    "I": 10.0,
    "M": 20.0,
    "E": 20.0,    # protonated -COOH
    "V": 24.0,
    "C": 30.0,    # neutral thiol
    "K": 35.0,    # neutral -NH2
    "D": 45.0,    # protonated -COOH
    "A": 45.0,
    "T": 50.0,
    "N": 50.0,
    "S": 50.0,
    "G": 55.0,
    "R": 60.0,
    "Q": 60.0,
}

# Ionised-state characteristic transition temperature (degC). Charging a side
# chain makes it far more hydrophilic, sharply raising T_tc.
URRY_TTC_CHARGED: dict[str, float] = {
    "E": 250.0,   # -COO(-)
    "D": 170.0,   # -COO(-)
    "K": 120.0,   # -NH3(+)
    "H": 90.0,    # imidazolium(+)
    "R": 60.0,    # guanidinium(+) (Arg is charged across the usable pH range)
    "C": 60.0,    # thiolate(-) (only above pH ~8.3)
    "Y": 40.0,    # phenolate(-) (only above pH ~10)
}

# Side-chain pKa values used to blend neutral/charged contributions.
# Sign convention: acids (D, E, C, Y) become charged ABOVE their pKa,
# bases (K, R, H) become charged BELOW their pKa.
SIDE_CHAIN_PKA: dict[str, float] = {
    "D": 3.65,
    "E": 4.25,
    "H": 6.00,
    "C": 8.30,
    "Y": 10.46,
    "K": 10.53,
    "R": 12.48,
}

ACIDIC = {"D", "E", "C", "Y"}
BASIC = {"K", "R", "H"}

# The 20 standard amino acids; X used as a wildcard guest placeholder.
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def fraction_charged(aa: str, pH: float) -> float:
    """Fraction of side chains that are ionised at the given pH (0..1)."""
    pKa = SIDE_CHAIN_PKA.get(aa)
    if pKa is None:
        return 0.0
    if aa in ACIDIC:
        # deprotonated (charged) fraction
        return 1.0 / (1.0 + 10.0 ** (pKa - pH))
    # basic: protonated (charged) fraction
    return 1.0 / (1.0 + 10.0 ** (pH - pKa))


def residue_ttc(aa: str, pH: float = 7.4) -> float:
    """Effective characteristic transition temperature (degC) of a guest residue.

    Blends the neutral and charged Urry values using the ionised fraction at
    the requested pH. Unknown / non-standard residues raise KeyError.
    """
    aa = aa.upper()
    neutral = URRY_TTC_NEUTRAL[aa]
    if aa not in SIDE_CHAIN_PKA:
        return neutral
    charged = URRY_TTC_CHARGED[aa]
    f = fraction_charged(aa, pH)
    return (1.0 - f) * neutral + f * charged


@dataclass(frozen=True)
class GuestComposition:
    """Mean Urry hydrophobicity descriptors of a set of guest residues."""

    mean_ttc: float          # composition-weighted characteristic T_t (degC)
    fraction_charged: float  # fraction of guest residues ionised at this pH
    fraction_aromatic: float # fraction that are F/Y/W
    n_guests: int


def guest_composition(guests: list[str], pH: float = 7.4) -> GuestComposition:
    """Summarise a list of guest residues into hydrophobicity descriptors."""
    if not guests:
        raise ValueError("No guest residues supplied.")
    guests = [g.upper() for g in guests]
    ttcs = [residue_ttc(g, pH) for g in guests]
    charged = sum(fraction_charged(g, pH) for g in guests) / len(guests)
    aromatic = sum(1 for g in guests if g in {"F", "Y", "W"}) / len(guests)
    return GuestComposition(
        mean_ttc=sum(ttcs) / len(ttcs),
        fraction_charged=charged,
        fraction_aromatic=aromatic,
        n_guests=len(guests),
    )
