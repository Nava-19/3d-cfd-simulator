"""
Tests unitarios para el módulo de Atmósfera Estándar Internacional (ISA 1976).
Valida las propiedades calculadas frente a las tablas oficiales de la OACI / NASA.
"""

import math
import pytest
from src.physics.atmosphere import AtmosphereISA


@pytest.fixture
def isa() -> AtmosphereISA:
    return AtmosphereISA()


def test_sea_level_conditions(isa: AtmosphereISA):
    """Valida los valores estándar al nivel del mar (h = 0 m)."""
    t, p, rho = isa.get_state_at_altitude(0.0)
    assert pytest.approx(t, rel=1e-4) == 288.15
    assert pytest.approx(p, rel=1e-4) == 101325.0
    assert pytest.approx(rho, rel=1e-3) == 1.2250


def test_troposphere_5000m(isa: AtmosphereISA):
    """Valida condiciones a 5000 m (Troposfera)."""
    t, p, rho = isa.get_state_at_altitude(5000.0)
    assert pytest.approx(t, rel=1e-4) == 255.65
    assert pytest.approx(p, rel=5e-3) == 54048.0
    assert pytest.approx(rho, rel=5e-3) == 0.7364


def test_tropopause_11000m(isa: AtmosphereISA):
    """Valida condiciones en la tropopausa a 11000 m."""
    t, p, rho = isa.get_state_at_altitude(11000.0)
    assert pytest.approx(t, rel=1e-4) == 216.65
    assert pytest.approx(p, rel=5e-3) == 22632.0
    assert pytest.approx(rho, rel=5e-3) == 0.3639


def test_lower_stratosphere_20000m(isa: AtmosphereISA):
    """Valida condiciones en la estratosfera baja isotérmica a 20000 m."""
    t, p, rho = isa.get_state_at_altitude(20000.0)
    assert pytest.approx(t, rel=1e-4) == 216.65
    assert pytest.approx(p, rel=5e-3) == 5474.9
    assert pytest.approx(rho, rel=5e-3) == 0.0880


def test_speed_of_sound_sea_level(isa: AtmosphereISA):
    """Valida la velocidad del sonido a nivel del mar (~340.3 m/s)."""
    cond = isa.compute_flight_condition(mach=1.0, altitude_km=0.0)
    assert pytest.approx(cond.speed_of_sound_m_s, rel=1e-3) == 340.294
    assert pytest.approx(cond.velocity_magnitude_m_s, rel=1e-3) == 340.294


def test_flight_condition_supersonic_with_aoa(isa: AtmosphereISA):
    """Valida el cálculo de vector velocidad a Mach 2.5 con alfa = 4 grados."""
    cond = isa.compute_flight_condition(mach=2.5, altitude_km=10.0, alpha_deg=4.0)

    # Comprobar módulo de velocidad
    expected_u_mag = 2.5 * cond.speed_of_sound_m_s
    assert pytest.approx(cond.velocity_magnitude_m_s, rel=1e-4) == expected_u_mag

    # Proyección cartesiana
    ux, uy, uz = cond.velocity_vector_m_s
    assert pytest.approx(ux, rel=1e-4) == expected_u_mag * math.cos(math.radians(4.0))
    assert pytest.approx(uy, rel=1e-4) == expected_u_mag * math.sin(math.radians(4.0))
    assert uz == 0.0

    # Temperatura de estancamiento Mach 2.5: T0 = T * (1 + 0.2 * 2.5^2) = T * 2.25
    assert pytest.approx(cond.stagnation_temperature_k, rel=1e-4) == cond.temperature_k * 2.25


def test_altitude_out_of_bounds(isa: AtmosphereISA):
    """Valida que lance error para altitudes negativas o mayores a 32 km."""
    with pytest.raises(ValueError):
        isa.get_state_at_altitude(-100.0)

    with pytest.raises(ValueError):
        isa.get_state_at_altitude(35000.0)
