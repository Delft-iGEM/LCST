"""Parse an ELP amino-acid sequence into pentapeptide repeats and guest residues.

Canonical ELP repeat: VPGXG (Val-Pro-Gly-Xaa-Gly), where Xaa ("the guest
residue") is any amino acid except proline. The transition temperature is set
mainly by the guest residues, the number of repeats, and the concentration.

The parser is tolerant of:
  * leading / trailing non-ELP residues (His-tags, linkers, Met start, etc.);
  * the repeat being written in any of the 5 possible reading frames;
  * minor motif deviations (e.g. IPGXG / APGXG / LPGXG variants), as long as
    the conserved ...PG.G... pattern is present.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

# Average monoisotopic-free (average) residue masses (Da), water already removed
# (i.e. residue masses in a polymer). Used only for reporting MW.
_RESIDUE_MASS = {
    "A": 71.08, "R": 156.19, "N": 114.10, "D": 115.09, "C": 103.14,
    "E": 129.12, "Q": 128.13, "G": 57.05, "H": 137.14, "I": 113.16,
    "L": 113.16, "K": 128.17, "M": 131.19, "F": 147.18, "P": 97.12,
    "S": 87.08, "T": 101.10, "V": 99.13, "W": 186.21, "Y": 163.18,
}
_WATER = 18.02


@dataclass
class ParsedELP:
    sequence: str            # cleaned input sequence
    guests: list[str]        # guest residue (X) of each detected pentad, in order
    n_pentads: int           # number of VPGXG repeats detected
    frame_offset: int        # offset (0-4) of the best reading frame
    motif_match_fraction: float  # share of pentads matching the ...PG.G consensus
    molecular_weight: float  # of the full cleaned sequence (Da)

    @property
    def is_elp(self) -> bool:
        return self.n_pentads >= 2 and self.motif_match_fraction >= 0.5


def _clean(sequence: str) -> str:
    seq = re.sub(r"\s+", "", sequence).upper()
    seq = seq.replace("-", "")
    if not seq:
        raise ValueError("Empty sequence.")
    bad = set(seq) - set("ACDEFGHIKLMNPQRSTVWYX")
    if bad:
        raise ValueError(f"Sequence contains non amino-acid characters: {sorted(bad)}")
    return seq


def _score_frame(seq: str, offset: int) -> tuple[int, int, list[str]]:
    """Score a reading frame by ...P.G.X.G consensus matches.

    Returns (n_matching_pentads, n_total_pentads, guest_residues).
    A pentad V-P-G-X-G matches the consensus when position 2 is P and
    positions 3 and 5 are G (1-indexed within the pentad).
    """
    guests: list[str] = []
    matches = 0
    total = 0
    for start in range(offset, len(seq) - 4, 5):
        pentad = seq[start : start + 5]
        if len(pentad) < 5:
            break
        total += 1
        # consensus: _ P G X G  -> pentad[1]=='P', pentad[2]=='G', pentad[4]=='G'
        if pentad[1] == "P" and pentad[2] == "G" and pentad[4] == "G":
            matches += 1
            guests.append(pentad[3])
    return matches, total, guests


def parse_elp(sequence: str) -> ParsedELP:
    """Parse a sequence into its ELP pentads and guest residues."""
    seq = _clean(sequence)

    best = None
    for offset in range(5):
        matches, total, guests = _score_frame(seq, offset)
        if total == 0:
            continue
        # Prefer the frame with the most consensus matches, then the most pentads.
        key = (matches, total)
        if best is None or key > best[0]:
            best = (key, offset, matches, total, guests)

    if best is None:
        raise ValueError(
            "Sequence too short to contain a single pentapeptide repeat."
        )

    _, offset, matches, total, all_guests = best
    # Keep only guests from consensus-matching pentads for the hydrophobicity calc.
    matching_guests = [
        seq[start + 3]
        for start in range(offset, len(seq) - 4, 5)
        if seq[start + 1] == "P" and seq[start + 2] == "G" and seq[start + 4] == "G"
    ]
    mw = sum(_RESIDUE_MASS[a] for a in seq if a in _RESIDUE_MASS) + _WATER

    return ParsedELP(
        sequence=seq,
        guests=matching_guests,
        n_pentads=len(matching_guests),
        frame_offset=offset,
        motif_match_fraction=(matches / total) if total else 0.0,
        molecular_weight=mw,
    )


def guests_from_repeat_spec(spec: str) -> list[str]:
    """Expand a compact guest spec like 'V5A2G3' into a guest-residue list.

    'V5A2G3' -> 5 Val + 2 Ala + 3 Gly guests (one repeat unit of 10 pentads).
    Useful for the classic Chilkoti ELP naming convention.
    """
    guests: list[str] = []
    for letter, count in re.findall(r"([A-Z])(\d+)", spec.upper()):
        guests.extend([letter] * int(count))
    if not guests:
        raise ValueError(f"Could not parse guest spec: {spec!r}")
    return guests
