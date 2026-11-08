"""
Tests unitarios y de integración para la Interfaz Gráfica de Usuario (Fase 8):
    - FlightControlPanel (reactividad en vivo con AtmosphereISA)
    - TelemetryPlotsPanel (streaming de métricas y lectura histórica de casos)
    - Viewport3DPanel y CFDMainWindow (integración completa 3D off-screen con PyVista/VTK)
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
    # En la tropopausa (11 km), T_inf = 216.65 K y T_0(M=3) = 216.65 * (1 + 0.2 * 9) = 606.62 K
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


def test_main_window_loads_and_renders_benchmark_3d_offscreen(qapp):
    window = CFDMainWindow(off_screen=True)

    # Al inicializarse, CFDMainWindow detecta y carga automáticamente cases/benchmark_cone_m2
    assert window.viewport_panel.surface_mesh is not None
    assert window.viewport_panel.schlieren_slice is not None
    assert window.viewport_panel.streamlines_mesh is not None
    assert window.viewport_panel.aero_metrics is not None

    # Verificar que las curvas de fuerzas del caso se cargaron en el panel inferior
    assert len(window.plots_panel.times_force) > 100
    assert len(window.plots_panel.cd_vals) > 100

    # Verificar que el reporte de validación Taylor-Maccoll se imprimió en la consola
    console_text = window.plots_panel.console.toPlainText()
    assert "Taylor-Maccoll" in console_text
