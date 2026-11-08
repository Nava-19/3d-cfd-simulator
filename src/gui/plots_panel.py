"""
Panel Inferior de Telemetría y Gráficas de Convergencia en Tiempo Real (Fase 8.3).
Integra pyqtgraph acelerado para trazar en directo:
    1. Estabilidad numérica (Courant máximo Co_max y residuos logarítmicos).
    2. Coeficientes aerodinámicos (CD, CL, Cm_pitch) estabilizándose hacia estado estacionario.
    3. Consola de telemetría WSL2 y tabla de validación analítica Taylor-Maccoll.
"""

import math
from pathlib import Path
import numpy as np
import pyqtgraph as pg
from PyQt6.QtWidgets import (
    QPlainTextEdit,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from src.backend.log_parser import IterationMetrics, OpenFOAMLogParser


# Configuración global de estilo aeroespacial oscuro para PyQtGraph
pg.setConfigOptions(antialias=True, background="#111827", foreground="#E5E7EB")


class TelemetryPlotsPanel(QWidget):
    """Panel inferior con pestañas para gráficas en vivo y consola de OpenFOAM."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)

        # Buffers de series temporales
        self.times_res: list[float] = []
        self.courant_vals: list[float] = []
        self.residuals_log: dict[str, list[float]] = {
            "Ux": [],
            "Uy": [],
            "e": [],
            "k": [],
            "omega": [],
        }

        self.times_force: list[float] = []
        self.cd_vals: list[float] = []
        self.cl_vals: list[float] = []
        self.cm_vals: list[float] = []

        self._build_ui()

    def _build_ui(self) -> None:
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self.tabs = QTabWidget()
        layout.addWidget(self.tabs)

        # ------------------------------------------------------------------
        # Pestaña 1: Coeficientes Aerodinámicos (CD, CL, Cm)
        # ------------------------------------------------------------------
        self.plot_forces = pg.PlotWidget(title="Evolución Temporal de Coeficientes Aerodinámicos (CD, CL, Cm)")
        self.plot_forces.setLabel("bottom", "Tiempo Físico (t)", units="s")
        self.plot_forces.setLabel("left", "Coeficiente Adimensional [-]")
        self.plot_forces.showGrid(x=True, y=True, alpha=0.25)
        self.plot_forces.addLegend(offset=(10, 10))

        self.curve_cd = self.plot_forces.plot(
            pen=pg.mkPen(color="#EF4444", width=2.5), name="CD (Resistencia Total)"
        )
        self.curve_cl = self.plot_forces.plot(
            pen=pg.mkPen(color="#3B82F6", width=2.0), name="CL (Sustentación)"
        )
        self.curve_cm = self.plot_forces.plot(
            pen=pg.mkPen(color="#10B981", width=2.0, style=pg.QtCore.Qt.PenStyle.DashLine),
            name="Cm,pitch (Momento Cabeceo)"
        )
        self.tabs.addTab(self.plot_forces, "📈 Coeficientes Aerodinámicos (CD / CL / Cm)")

        # ------------------------------------------------------------------
        # Pestaña 2: Estabilidad Numérica (CFL / Courant) y Residuos
        # ------------------------------------------------------------------
        self.plot_stability = pg.PlotWidget(
            title="Estabilidad Numérica (Courant CFL) y Residuos Logarítmicos log10(R)"
        )
        self.plot_stability.setLabel("bottom", "Tiempo Físico (t)", units="s")
        self.plot_stability.setLabel("left", "Courant Max / log10(Residuo)")
        self.plot_stability.showGrid(x=True, y=True, alpha=0.25)
        self.plot_stability.addLegend(offset=(10, 10))

        # Línea horizontal de límite CFL = 0.75 (Kurganov-Tadmor)
        self.cfl_limit_line = pg.InfiniteLine(
            pos=0.75,
            angle=0,
            pen=pg.mkPen(color="#F87171", width=1.5, style=pg.QtCore.Qt.PenStyle.DotLine),
            label="Límite CFL = 0.75",
            labelOpts={"position": 0.9, "color": "#F87171"}
        )
        self.plot_stability.addItem(self.cfl_limit_line)

        self.curve_courant = self.plot_stability.plot(
            pen=pg.mkPen(color="#FBBF24", width=2.2), name="Courant Máximo (Co_max)"
        )
        self.res_curves: dict[str, pg.PlotDataItem] = {
            "Ux": self.plot_stability.plot(pen=pg.mkPen("#60A5FA", width=1.5), name="log10(R_Ux)"),
            "Uy": self.plot_stability.plot(pen=pg.mkPen("#34D399", width=1.5), name="log10(R_Uy)"),
            "e": self.plot_stability.plot(pen=pg.mkPen("#F472B6", width=1.5), name="log10(R_e)"),
            "k": self.plot_stability.plot(pen=pg.mkPen("#A78BFA", width=1.5), name="log10(R_k)"),
            "omega": self.plot_stability.plot(pen=pg.mkPen("#FB923C", width=1.5), name="log10(R_omega)"),
        }
        self.tabs.addTab(self.plot_stability, "🔬 Convergencia y Número de Courant (CFL)")

        # ------------------------------------------------------------------
        # Pestaña 3: Consola de Streaming WSL2 y Validación Analítica
        # ------------------------------------------------------------------
        self.console = QPlainTextEdit()
        self.console.setReadOnly(True)
        self.console.setMaximumBlockCount(4000)
        self.console.setStyleSheet(
            "background-color: #0B0F19; color: #D1D5DB; font-family: Consolas, monospace; font-size: 12px;"
        )
        self.tabs.addTab(self.console, "🖥️ Consola OpenFOAM & Reporte Taylor-Maccoll")

    def clear_telemetry(self) -> None:
        """Limpia todos los buffers históricos y reinicia los lienzos."""
        self.times_res.clear()
        self.courant_vals.clear()
        for k in self.residuals_log:
            self.residuals_log[k].clear()

        self.times_force.clear()
        self.cd_vals.clear()
        self.cl_vals.clear()
        self.cm_vals.clear()

        self.curve_cd.setData([], [])
        self.curve_cl.setData([], [])
        self.curve_cm.setData([], [])
        self.curve_courant.setData([], [])
        for curve in self.res_curves.values():
            curve.setData([], [])

        self.console.clear()

    def append_log_line(self, line: str) -> None:
        """Añade una línea de texto a la consola en vivo."""
        self.console.appendPlainText(line)

    def append_iteration(self, metrics: IterationMetrics) -> None:
        """Actualiza las curvas en tiempo real al recibir un paquete IterationMetrics."""
        self.times_res.append(metrics.time)
        self.courant_vals.append(metrics.courant_max)

        for field_name, res_list in self.residuals_log.items():
            val = metrics.residuals.get(field_name)
            if val is not None and val > 0.0:
                res_list.append(math.log10(max(val, 1e-16)))
            elif res_list:
                res_list.append(res_list[-1])
            else:
                res_list.append(float("nan"))

        if metrics.cd is not None:
            self.times_force.append(metrics.time)
            self.cd_vals.append(metrics.cd)
            self.cl_vals.append(metrics.cl if metrics.cl is not None else 0.0)
            self.cm_vals.append(metrics.cm_pitch if metrics.cm_pitch is not None else 0.0)

        # Refrescar gráficas cada pocas iteraciones para mantener 60 FPS en la GUI
        if len(self.times_res) % 5 == 0 or len(self.times_res) <= 5:
            self._refresh_plots()

    def _refresh_plots(self) -> None:
        t_res = np.asarray(self.times_res, dtype=float)
        self.curve_courant.setData(t_res, np.asarray(self.courant_vals, dtype=float))

        for field_name, curve in self.res_curves.items():
            arr = np.asarray(self.residuals_log[field_name], dtype=float)
            valid = ~np.isnan(arr)
            if np.any(valid):
                curve.setData(t_res[valid], arr[valid])

        if self.times_force:
            t_f = np.asarray(self.times_force, dtype=float)
            self.curve_cd.setData(t_f, np.asarray(self.cd_vals, dtype=float))
            self.curve_cl.setData(t_f, np.asarray(self.cl_vals, dtype=float))
            self.curve_cm.setData(t_f, np.asarray(self.cm_vals, dtype=float))

    def load_history_from_case(self, case_dir: Path | str) -> None:
        """
        Lee los archivos históricos (coefficient.dat y log.rhoCentralFoam) de un caso
        ya simulado en disco para poblar las gráficas instantáneamente.
        """
        self.clear_telemetry()
        cpath = Path(case_dir)

        # 1. Leer log.rhoCentralFoam si existe
        log_file = cpath / "log.rhoCentralFoam"
        if log_file.exists():
            parser = OpenFOAMLogParser()
            with open(log_file, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    m = parser.parse_line(line)
                    if m is not None:
                        self.times_res.append(m.time)
                        self.courant_vals.append(m.courant_max)
                        for field_name, res_list in self.residuals_log.items():
                            val = m.residuals.get(field_name)
                            if val is not None and val > 0.0:
                                res_list.append(math.log10(max(val, 1e-16)))
                            else:
                                res_list.append(float("nan"))
                        if m.cd is not None:
                            self.times_force.append(m.time)
                            self.cd_vals.append(m.cd)
                            self.cl_vals.append(m.cl if m.cl is not None else 0.0)
                            self.cm_vals.append(m.cm_pitch if m.cm_pitch is not None else 0.0)

        # 2. Si no se extrajeron fuerzas del log, buscar en postProcessing/forces_breakdown/0/coefficient.dat
        if not self.times_force:
            coeff_files = list(cpath.glob("postProcessing/**/coefficient.dat"))
            if coeff_files:
                with open(coeff_files[0], "r", encoding="utf-8", errors="ignore") as f:
                    for line in f:
                        stripped = line.strip()
                        if not stripped or stripped.startswith("#"):
                            continue
                        parts = stripped.split()
                        if len(parts) >= 8:
                            try:
                                self.times_force.append(float(parts[0]))
                                self.cd_vals.append(float(parts[1]))
                                self.cl_vals.append(float(parts[4]))
                                self.cm_vals.append(float(parts[7]))
                            except ValueError:
                                continue

        self._refresh_plots()
