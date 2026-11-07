"""
Tests unitarios para el solver analítico de Taylor-Maccoll.
Compara la solución analítica con tablas clásicas de flujo compresible supersónico (NACA 1135).
"""

import pytest
from src.physics.taylor_maccoll import TaylorMaccollSolver


def test_taylor_maccoll_cone_15deg_mach_2():
    """Valida un cono de 15 grados a Mach 2.0 frente a NACA 1135."""
    solver = TaylorMaccollSolver(gamma=1.4)
    res = solver.solve(mach_inf=2.0, cone_angle_deg=15.0)

    # Ángulo de choque teórico débil: ~33.9 grados
    assert 33.0 < res.shock_angle_deg < 35.0

    # Relación de presión en superficie NACA 1135: ~1.57
    assert 1.50 < res.pressure_ratio_surface < 1.65

    # Mach en superficie: ~1.71
    assert 1.65 < res.mach_surface < 1.85


def test_taylor_maccoll_cone_20deg_mach_3():
    """Valida un cono de 20 grados a Mach 3.0 frente a NACA 1135."""
    solver = TaylorMaccollSolver(gamma=1.4)
    res = solver.solve(mach_inf=3.0, cone_angle_deg=20.0)

    # Ángulo de choque teórico débil: ~29.6 grados
    assert 28.5 < res.shock_angle_deg < 30.5

    # Relación de presión en superficie para cono 20 deg a Mach 3 NACA 1135: ~3.06
    assert 2.95 < res.pressure_ratio_surface < 3.20

    # Mach en superficie: ~2.29
    assert 2.20 < res.mach_surface < 2.40
