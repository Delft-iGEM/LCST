"""Optional machine-learning calibration layer.

The physics model (`model.predict`) is unbiased for Val/Ala hydrophobic ELPs
but can carry a systematic offset for compositions outside that range
(e.g. Gly-rich or charged guests). This module learns a *residual correction*

    Tt_measured - Tt_physics = f(features)

from experimental data, so the physics model still governs extrapolation while
the ML layer removes systematic bias where data exist. Training on data that
already match the physics model yields a ~zero correction, i.e. the calibrated
model gracefully reduces to the physics model.
"""

from __future__ import annotations

import pickle
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import KFold, cross_val_predict

from . import chilkoti
from .model import predict as physics_predict
from .sequence import parse_elp
from .urry import guest_composition

FEATURES = [
    "mean_guest_ttc",
    "log10_n",
    "log10_conc",
    "fraction_charged",
    "fraction_aromatic",
    "pH",
]


def featurize(mean_ttc, n_pentads, conc_uM, pH, frac_charged, frac_aromatic):
    return {
        "mean_guest_ttc": mean_ttc,
        "log10_n": np.log10(n_pentads),
        "log10_conc": np.log10(conc_uM),
        "fraction_charged": frac_charged,
        "fraction_aromatic": frac_aromatic,
        "pH": pH,
    }


def features_from_sequence(sequence: str, conc_uM: float, pH: float) -> dict:
    parsed = parse_elp(sequence)
    comp = guest_composition(parsed.guests, pH=pH)
    return featurize(
        comp.mean_ttc, parsed.n_pentads, conc_uM, pH,
        comp.fraction_charged, comp.fraction_aromatic,
    )


def generate_unified_dataset() -> pd.DataFrame:
    """High-fidelity Val/Ala training grid from the Chilkoti unified equation.

    These rows reproduce the published model (r^2 > 0.99 vs experiment across
    120 measurements), giving dense, trustworthy coverage of the
    composition / length / concentration space for hydrophobic ELPs.
    """
    rows = []
    fas = [0.0, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    lengths = [20, 30, 60, 90, 120, 150, 180]
    concs = [1, 5, 10, 25, 50, 100, 250, 500]
    for fa in fas:
        mean_ttc = chilkoti.ttc_of_fa(fa)
        for n in lengths:
            for c in concs:
                tt = chilkoti.transition_temperature(mean_ttc, n, c)
                rows.append(
                    dict(
                        source="unified_model(VA)",
                        mean_guest_ttc=mean_ttc,
                        n_pentads=n,
                        concentration_uM=c,
                        pH=7.4,
                        fraction_charged=0.0,
                        fraction_aromatic=0.0,
                        tt_celsius=tt,
                    )
                )
    return pd.DataFrame(rows)


@dataclass
class CalibratedModel:
    regressor: object
    cv_mae: float
    cv_rmse: float
    n_train: int

    def predict_tt(self, sequence: str, concentration_uM=25.0, pH=7.4) -> float:
        base = physics_predict(sequence, concentration_uM, pH).tt_celsius
        feats = features_from_sequence(sequence, concentration_uM, pH)
        x = np.array([[feats[f] for f in FEATURES]])
        return float(base + self.regressor.predict(x)[0])

    def save(self, path: str | Path) -> None:
        Path(path).write_bytes(pickle.dumps(self))

    @staticmethod
    def load(path: str | Path) -> "CalibratedModel":
        return pickle.loads(Path(path).read_bytes())


def _physics_tt_for_row(row) -> float:
    """Physics Tt for a dataset row that may specify features directly."""
    p = chilkoti.transition_temperature(
        row["mean_guest_ttc"], int(row["n_pentads"]), row["concentration_uM"]
    )
    return p


def train(df: pd.DataFrame, random_state: int = 0) -> CalibratedModel:
    """Train the residual-correction regressor on a dataset.

    Required columns: mean_guest_ttc, n_pentads, concentration_uM, pH,
    fraction_charged, fraction_aromatic, tt_celsius.
    (Use `dataset_from_sequences` to build these from raw sequences.)
    """
    df = df.copy()
    physics = np.array([_physics_tt_for_row(r) for _, r in df.iterrows()])
    residual = df["tt_celsius"].to_numpy() - physics

    X = np.column_stack([
        df["mean_guest_ttc"].to_numpy(),
        np.log10(df["n_pentads"].to_numpy()),
        np.log10(df["concentration_uM"].to_numpy()),
        df["fraction_charged"].to_numpy(),
        df["fraction_aromatic"].to_numpy(),
        df["pH"].to_numpy(),
    ])

    reg = GradientBoostingRegressor(
        n_estimators=300, max_depth=2, learning_rate=0.05,
        subsample=0.9, random_state=random_state,
    )

    n_splits = min(5, len(df))
    if n_splits >= 2:
        kf = KFold(n_splits=n_splits, shuffle=True, random_state=random_state)
        resid_cv = cross_val_predict(reg, X, residual, cv=kf)
        tt_cv = physics + resid_cv
        err = tt_cv - df["tt_celsius"].to_numpy()
        cv_mae = float(np.mean(np.abs(err)))
        cv_rmse = float(np.sqrt(np.mean(err**2)))
    else:
        cv_mae = cv_rmse = float("nan")

    reg.fit(X, residual)
    return CalibratedModel(reg, cv_mae, cv_rmse, len(df))


def dataset_from_sequences(df: pd.DataFrame, pH_default: float = 7.4) -> pd.DataFrame:
    """Build a training frame from rows of (sequence, concentration_uM, tt_celsius).

    Optional columns: pH. Computes the guest descriptors from each sequence.
    """
    out = []
    for _, r in df.iterrows():
        pH = float(r.get("pH", pH_default))
        feats = features_from_sequence(r["sequence"], float(r["concentration_uM"]), pH)
        parsed = parse_elp(r["sequence"])
        comp = guest_composition(parsed.guests, pH=pH)
        out.append(dict(
            source=r.get("source", "user"),
            sequence=r["sequence"],
            mean_guest_ttc=comp.mean_ttc,
            n_pentads=parsed.n_pentads,
            concentration_uM=float(r["concentration_uM"]),
            pH=pH,
            fraction_charged=comp.fraction_charged,
            fraction_aromatic=comp.fraction_aromatic,
            tt_celsius=float(r["tt_celsius"]),
        ))
    return pd.DataFrame(out)
