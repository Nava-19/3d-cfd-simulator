"""
Tests unitarios y de integración para la Interfaz Gráfica de Usuario (Fase 8 y 8.5):
    - FlightControlPanel (reactividad en vivo con AtmosphereISA)
    - TelemetryPlotsPanel (streaming de métricas y lectura histórica de casos)
    - Viewport3DPanel y CFDMainWindow (cono 360°, volumen de control, partículas en vivo)
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


def test_main_window_loads_360_cone_particles_and_control_volume(qapp):
    window = CFDMainWindow(off_screen=True)
    vp = window.viewport_panel

    assert vp.surface_mesh is not None
    assert vp.schlieren_slice is not None
    assert vp.streamlines_mesh is not None
    assert vp.control_volume_mesh is not None
    assert vp.particle_system is not None

    # Verificar que el cono 360° cubre simétricamente +Y, -Y, +Z y -Z
    _, _, ymin, ymax, zmin, zmax = vp.surface_mesh.bounds
    assert ymin < -0.10 and ymax > 0.10
    assert zmin < -0.10 and zmax > 0.10

    # Probar conmutación a Corte Seccional 180° y verificar que sigue cubriendo +Y y -Y
    vp.combo_body_mode.setCurrentIndex(1)
    _, _, ymin_180, ymax_180, _, _ = vp.surface_mesh.bounds
    assert ymin_180 < -0.10 and ymax_180 > 0.10

    # Probar un tick de animación de partículas
    vp._on_animation_tick()

    assert len(window.plots_panel.times_force) > 100
    assert "Taylor-Maccoll" in window.plots_panel.console.toPlainText()
