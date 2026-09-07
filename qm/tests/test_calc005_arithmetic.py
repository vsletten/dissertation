"""Analytical sign/unit/detailed-balance gates for CALC-005."""

import math

import pytest

from quarry.calc005 import (
    detailed_balance_rates,
    standard_state_1m_from_1bar_kj_mol,
)


def test_one_molar_standard_state_correction_has_declared_sign_and_value():
    correction = standard_state_1m_from_1bar_kj_mol(temperature_k=298.15)

    assert correction == pytest.approx(7.958500693927389, abs=1e-12)


def test_positive_detachment_stabilization_slows_reverse_rate():
    rates = detailed_balance_rates(
        attachment_rate_s=1.7e5,
        detachment_free_energy_kj_mol=42.0,
        temperature_k=298.15,
        activity=0.037,
    )

    assert rates["detachment_rate_s"] < rates["attachment_rate_s"]
    assert rates["kinetic_ratio"] == pytest.approx(
        rates["thermodynamic_ratio"], rel=1e-14
    )
    assert rates["detachment_rate_s"] == pytest.approx(
        1.7e5 * math.exp(-42.0 / (rates["r_kj_mol_k"] * 298.15))
    )


@pytest.mark.parametrize("field", ["attachment_rate_s", "temperature_k", "activity"])
def test_detailed_balance_rejects_nonpositive_inputs(field):
    values = {
        "attachment_rate_s": 1.0,
        "detachment_free_energy_kj_mol": 0.0,
        "temperature_k": 298.15,
        "activity": 1.0,
    }
    values[field] = 0.0
    with pytest.raises(ValueError, match=field):
        detailed_balance_rates(**values)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("attachment_rate_s", float("nan")),
        ("detachment_free_energy_kj_mol", float("nan")),
        ("temperature_k", float("inf")),
        ("activity", True),
    ],
)
def test_detailed_balance_rejects_nonfinite_and_boolean_inputs(field, value):
    values = {
        "attachment_rate_s": 1.0,
        "detachment_free_energy_kj_mol": 0.0,
        "temperature_k": 298.15,
        "activity": 1.0,
    }
    values[field] = value
    with pytest.raises((TypeError, ValueError), match=field):
        detailed_balance_rates(**values)


@pytest.mark.parametrize("value", [float("nan"), float("inf"), True])
def test_standard_state_rejects_nonfinite_or_boolean_temperature(value):
    with pytest.raises((TypeError, ValueError), match="temperature_k"):
        standard_state_1m_from_1bar_kj_mol(temperature_k=value)


def test_standard_state_rejects_finite_inputs_with_overflowing_result():
    with pytest.raises(ValueError, match="nonfinite"):
        standard_state_1m_from_1bar_kj_mol(
            temperature_k=1e308,
            concentration_mol_m3=1e308,
            pressure_pa=1.0,
        )
