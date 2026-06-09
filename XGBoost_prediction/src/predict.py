"""
predict.py

This file loads the trained LCST prediction model and predicts LCST
for a new ELP sequence.
"""

import joblib
import pandas as pd

from features import make_feature_table


def predict_lcst(sequence, concentration_uM=25, salt_mM=150, pH=7.4):
    """
    Predicts the LCST temperature of a new ELP sequence.

    Parameters
    ----------
    sequence : str
        Full amino acid sequence.
    concentration_uM : float
        Polymer concentration in micromolar.
    salt_mM : float
        Salt concentration in millimolar.
    pH : float
        Solution pH.

    Returns
    -------
    float
        Predicted LCST in degrees Celsius.
    """

    model_path = "models/lcst_xgboost_model.joblib"

    # Load trained model
    model = joblib.load(model_path)

    # Put the input sequence and conditions into a dataframe
    new_data = pd.DataFrame([
        {
            "sequence": sequence,
            "concentration_uM": concentration_uM,
            "salt_mM": salt_mM,
            "pH": pH,
        }
    ])

    # Convert sequence to numerical features
    X_new = make_feature_table(new_data)

    # Predict LCST
    prediction = model.predict(X_new)

    return prediction[0]


if __name__ == "__main__":
    # Example ELP sequence.
    # This is VPGVG repeated 4 times.
    example_sequence = "VPGVGVPGVGVPGVGVPGVG"

    predicted_lcst = predict_lcst(
        sequence=example_sequence,
        concentration_uM=25,
        salt_mM=150,
        pH=7.4
    )

    print(f"Predicted LCST: {predicted_lcst:.2f} °C")