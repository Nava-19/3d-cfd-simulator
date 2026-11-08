"""
Tests unitarios y de integración para la Interfaz Gráfica de Usuario (Fase 8 y 8.5):
    - FlightControlPanel (reactividad en vivo con AtmosphereISA)
    - TelemetryPlotsPanel (streaming de métricas y lectura histórica de casos)
    - Viewport3DPanel y CFDMainWindow (cono 360°, esfera con Bow Shock, cono 3D de choque, 16 planos de partículas)
"""

import os
import sys
import math

os.environ["QT_API"] = "pyqt6"

from PyQt6.QtWidgets import QApplication
import pytest

from src.backend.log_parser import IterationMetrics
from src.gui.control_panel import FlightControlPanel
from src.gui.main_window import CFDMainWindow
from src.gui.plots_panel import TelemetryPlotsPanel


@pytest.fixture(scope="session")
def qapp():
    """Garantiza una única instancia de QApplication para toda la sesión de pytest."""
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    return app


def test_control_panel_updates_isa_in_realtime(qapp):
    panel = FlightControlPanel()
    panel.spin_mach.setValue(3.0)
    panel.spin_altitude.setValue(11.0)

    cond = panel.current_condition
    assert cond is not None
    assert math.isclose(cond.mach, 3.0)
    assert math.isclose(cond.altitude_m, 11000.0)
    assert math.isclose(cond.temperature_k, 216.65, rel_tol=1e-3)
    assert math.isclose(cond.stagnation_temperature_k, 606.62, rel_tol=1e-3)


def test_plots_panel_streaming_and_history_loading(qapp):
    plots = TelemetryPlotsPanel()
    m1 = IterationMetrics(
        time=1e-4,
        courant_max=0.48,
        residuals={"Ux": 1e-3, "Uy": 2e-3, "e": 5e-4},
        cd=0.145,
        cl=0.012,
        cm_pitch=-0.004
    )
    plots.append_iteration(m1)
    assert len(plots.times_res) == 1
    assert math.isclose(plots.cd_vals[0], 0.145)


def test_main_window_loads_cone_and_sphere_with_3d_shock_and_particles(qapp):
    window = CFDMainWindow(off_screen=True)
    vp = window.viewport_panel

    # 1. Verificar carga inicial del Cono Supersónico M=2.0
    assert vp.surface_mesh is not None
    assert vp.schlieren_slice is not None
    assert vp.streamlines_mesh is not None
    assert vp.control_volume_mesh is not None
    assert vp.theo_shock_surface is not None
    assert vp.particle_system is not None
    assert vp.particle_system.cloud_polydata.n_points >= 500

    # 2. Cargar el caso de la Esfera Supersónica (Bow Shock + Estela) mediante el botón del panel
    window.control_panel._on_load_sphere_benchmark()
    assert vp.is_blunt_sphere is True
    assert vp.theo_shock_surface is not None
    # El vértice del Bow Shock desprendido en la esfera debe estar aguas arriba de x=0 (x < -0.02 m)
    xmin_shock = vp.theo_shock_surface.bounds[0]
    assert xmin_shock < -0.02
    assert "Bow Shock" in window.plots_panel.console.toPlainText()
