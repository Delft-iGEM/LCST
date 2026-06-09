"""Validate the physics model against the literature dataset and Table 1."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from . import chilkoti
from .model import predict
from .sequence import guests_from_repeat_spec

_DATA = Path(__file__).resolve().parents[2] / "data" / "literature_validation.csv"


def _sequence_from_spec(guest_spec: str, n_repeats: int) -> str:
    guests = guests_from_repeat_spec(guest_spec) * n_repeats
    return "".join(f"VPG{g}G" for g in guests)


def validate_literature(path: str | Path = _DATA) -> pd.DataFrame:
    df = pd.read_csv(path)
    out = []
    for _, r in df.iterrows():
        seq = _sequence_from_spec(r["guest_spec"], int(r["n_repeats"]))
        pred = predict(seq, concentration_uM=float(r["concentration_uM"]),
                       pH=float(r["pH"]))
        out.append(dict(
            name=r["name"],
            n_pentads=pred.n_pentads,
            measured_tt=float(r["tt_celsius"]),
            predicted_tt=pred.tt_celsius,
            residual=round(pred.tt_celsius - float(r["tt_celsius"]), 1),
        ))
    return pd.DataFrame(out)


def report() -> str:
    lines = []
    fq = chilkoti.fit_quality()
    lines.append("Composition fit quality (leave-one-out interpolation of Chilkoti Table 1):")
    for k, v in fq.items():
        lines.append(f"  {k:>4s}(f_A): R^2 = {v:.3f}")
    lines.append("")
    df = validate_literature()
    lines.append("Independent literature check:")
    for _, r in df.iterrows():
        lines.append(
            f"  {r['name']:<18s} n={r['n_pentads']:<4d} "
            f"measured={r['measured_tt']:5.1f}  "
            f"predicted={r['predicted_tt']:5.1f}  "
            f"residual={r['residual']:+5.1f} degC"
        )
    mae = df["residual"].abs().mean()
    lines.append("")
    lines.append(f"Mean absolute residual (independent points): {mae:.1f} degC")
    return "\n".join(lines)
