import math

import pytest

from elp_lcst import chilkoti
from elp_lcst.model import predict


def test_monotonic_in_hydrophobicity():
    # More hydrophobic guests (Phe) -> lower Tt than polar guests (Ser).
    phe = predict("VPGFG" * 30).tt_celsius
    val = predict("VPGVG" * 30).tt_celsius
    ser = predict("VPGSG" * 30).tt_celsius
    assert phe < val < ser


def test_tt_decreases_with_length():
    short = predict("VPGVG" * 30).tt_celsius
    long = predict("VPGVG" * 150).tt_celsius
    assert long < short


def test_tt_decreases_with_concentration():
    lo = predict("VPGVG" * 90, concentration_uM=1).tt_celsius
    hi = predict("VPGVG" * 90, concentration_uM=500).tt_celsius
    assert hi < lo


def test_reproduces_unified_model_table1():
    # The interpolators pass through Chilkoti Table 1 exactly at the knots.
    assert chilkoti.ttc_of_fa(0.0) == pytest.approx(20.2, abs=0.05)
    assert chilkoti.ttc_of_fa(1.0) == pytest.approx(54.3, abs=0.05)
    assert chilkoti.k_of_fa(0.5) == pytest.approx(316.9, abs=0.05)
    assert chilkoti.cc_of_fa(0.7) == pytest.approx(2748.4, abs=0.05)


def test_charged_guest_raises_tt_and_is_ph_sensitive():
    # Glutamate guest: deprotonated (charged) at pH 7 -> very high Tt;
    # protonated at low pH -> much lower Tt.
    high_ph = predict("VPGEG" * 40, pH=7.4).tt_celsius
    low_ph = predict("VPGEG" * 40, pH=3.0).tt_celsius
    assert high_ph > low_ph


def test_fit_quality_reasonable():
    fq = chilkoti.fit_quality()
    assert fq["Ttc"] > 0.9
    assert fq["k"] > 0.95
    assert fq["Cc"] > 0.9
