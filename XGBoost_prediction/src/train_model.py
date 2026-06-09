"""
train_model.py

This file trains an XGBoost regression model to predict LCST temperature.
"""

import os
import joblib
import pandas as pd

from xgboost import XGBRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, mean_squared_error

from features import make_feature_table


def train_model():
    """
    Loads LCST data, converts sequences to features, trains XGBoost,
    evaluates the model, and saves it.
    """

    # File paths
    data_path = "data/elp_lcst_data.csv"
    model_path = "models/lcst_xgboost_model.joblib"

    # Load data
    data = pd.read_csv(data_path)

    # Check that required columns exist
    required_columns = ["sequence", "LCST_C"]

    for column in required_columns:
        if column not in data.columns:
            raise ValueError(f"Missing required column: {column}")

    # Convert sequences into numerical features
    X = make_feature_table(data)

    # Target/output value
    y = data["LCST_C"]

    # Split into training and test data
    # 80% training, 20% testing
    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=0.2,
        random_state=1
    )

    # Create the XGBoost regression model
    model = XGBRegressor(
        n_estimators=300,
        learning_rate=0.05,
        max_depth=3,
        random_state=1,
        objective="reg:squarederror"
    )

    # Train the model
    model.fit(X_train, y_train)

    # Make predictions on the test set
    predictions = model.predict(X_test)

    # Evaluate the model
    mae = mean_absolute_error(y_test, predictions)
    mse = mean_squared_error(y_test, predictions)
    rmse = mse ** 0.5

    print("Model evaluation:")
    print(f"Mean Absolute Error: {mae:.2f} °C")
    print(f"Root Mean Squared Error: {rmse:.2f} °C")

    # Save the trained model
    os.makedirs("models", exist_ok=True)
    joblib.dump(model, model_path)

    print(f"Model saved to: {model_path}")


if __name__ == "__main__":
    train_model()