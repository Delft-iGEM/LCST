"""
features.py

This file contains helper functions for converting ELP amino acid sequences
into numerical features for machine learning.
"""

import pandas as pd


# Basic amino acid hydrophobicity values.
# These are Kyte-Doolittle hydrophobicity values.
HYDROPHOBICITY = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
}


# Approximate amino acid charges at neutral pH.
CHARGE = {
    "D": -1,
    "E": -1,
    "K": 1,
    "R": 1,
    "H": 0.1,
}


AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")


def clean_sequence(sequence):
    """
    Cleans an amino acid sequence.

    Example:
    '(VPGVG)40' should NOT be given to this function directly.
    This function expects the fully written sequence, like:
    'VPGVGVPGVGVPGVG'

    Parameters
    ----------
    sequence : str
        Amino acid sequence.

    Returns
    -------
    str
        Cleaned uppercase sequence.
    """
    return sequence.replace(" ", "").replace("\n", "").upper()


def sequence_to_features(sequence):
    """
    Converts one amino acid sequence into a dictionary of numerical features.

    Parameters
    ----------
    sequence : str
        Full amino acid sequence.

    Returns
    -------
    dict
        Numerical features describing the sequence.
    """

    sequence = clean_sequence(sequence)
    length = len(sequence)

    if length == 0:
        raise ValueError("Sequence is empty.")

    features = {}

    # Basic length feature
    features["sequence_length"] = length

    # Amino acid fractions
    for aa in AMINO_ACIDS:
        features[f"fraction_{aa}"] = sequence.count(aa) / length

    # Average hydrophobicity
    hydrophobicity_sum = 0
    for aa in sequence:
        hydrophobicity_sum += HYDROPHOBICITY.get(aa, 0)

    features["avg_hydrophobicity"] = hydrophobicity_sum / length

    # Net charge
    net_charge = 0
    for aa in sequence:
        net_charge += CHARGE.get(aa, 0)

    features["net_charge"] = net_charge

    # Fraction of charged residues
    charged_count = 0
    for aa in sequence:
        if aa in ["D", "E", "K", "R", "H"]:
            charged_count += 1

    features["fraction_charged"] = charged_count / length

    # ELP-specific features
    features["num_VPGXG_motifs"] = count_vpgxg_motifs(sequence)

    return features


def count_vpgxg_motifs(sequence):
    """
    Counts how many VPGXG motifs appear in the sequence.

    VPGXG means:
    V - P - G - any amino acid - G

    Example:
    VPGVG counts as one motif.
    VPGAG counts as one motif.
    VPGKG counts as one motif.
    """

    sequence = clean_sequence(sequence)
    count = 0

    for i in range(len(sequence) - 4):
        five_mer = sequence[i:i + 5]

        if (
            five_mer[0] == "V"
            and five_mer[1] == "P"
            and five_mer[2] == "G"
            and five_mer[4] == "G"
        ):
            count += 1

    return count


def make_feature_table(dataframe):
    """
    Converts a dataframe with ELP sequences into a machine learning feature table.

    The input dataframe must have a column called 'sequence'.

    Optional experimental columns can include:
    - concentration_uM
    - salt_mM
    - pH

    Parameters
    ----------
    dataframe : pandas.DataFrame
        Data table containing at least a 'sequence' column.

    Returns
    -------
    pandas.DataFrame
        Feature table for machine learning.
    """

    all_features = []

    for sequence in dataframe["sequence"]:
        features = sequence_to_features(sequence)
        all_features.append(features)

    feature_table = pd.DataFrame(all_features)

    # Add experimental condition columns if they exist
    optional_columns = ["concentration_uM", "salt_mM", "pH"]

    for column in optional_columns:
        if column in dataframe.columns:
            feature_table[column] = dataframe[column]

    return feature_table