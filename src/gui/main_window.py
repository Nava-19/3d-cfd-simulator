"""
Ventana Principal de la Estación de Simulación CFD 3D (Fase 8.1).
Orquesta la comunicación reactiva (Señales y Slots Qt) entre:
    - FlightControlPanel (Atmósfera ISA 1976 + Controles de Vuelo)
    - OpenFOAMRunner (Ejecución asíncrona en WSL2 en hilo QThread independiente)
    - TelemetryPlotsPanel (Gráficas PyQtGraph en tiempo real)
    - Viewport3DPanel (Motor de renderizado 3D GPU PyVista/VTK)
    - SupersonicConeBenchmark & AeroDBExporter (Validación analítica y exportación GNC)
"""

from pathlib import Path
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QLabel,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QSplitter,
    QStatusBar,
    QVBoxLayout,
    QWidget,
)

from src.backend.case_generator import OpenFOAMCaseGenerator
from src.backend.runner import OpenFOAMRunner
from src.gui.control_panel import FlightControlPanel
from src.gui.plots_panel import TelemetryPlotsPanel
from src.gui.viewport_panel import Viewport3DPanel
from src.physics.taylor_maccoll import TaylorMaccollSolver
from src.postprocess.aerodb_exporter import AeroDBExporter, AeroDBRecord
from src.validation.cone_benchmark import SupersonicConeBenchmark


DARK_AEROSPACE_STYLESHEET = """
QMainWindow, QWidget {
    background-color: #111827;
    color: #F3F4F6;
    font-family: "Segoe UI", Arial, sans-serif;
    font-size: 12px;
}
QGroupBox {
    border: 1px solid #374151;
    border-radius: 6px;
    margin-top: 8px;
    padding-top: 10px;
    font-weight: bold;
    color: #93C5FD;
}
QGroupBox::title {
    subcontrol-origin: margin;
    left: 10px;
    padding: 0 4px;
}
QDoubleSpinBox, QComboBox {
    background-color: #1F2937;
    border: 1px solid #4B5563;
    border-radius: 4px;
    padding: 4px 6px;
    color: #F9FAFB;
}
QPushButton {
    background-color: #1F2937;
    border: 1px solid #4B5563;
    border-radius: 4px;
    padding: 6px 10px;
    color: #E5E7EB;
}
QPushButton:hover {
    background-color: #374151;
    border-color: #60A5FA;
}
QPushButton:disabled {
    background-color: #111827;
    color: #6B7280;
    border-color: #374151;
}
QTabWidget::pane {
    border: 1px solid #374151;
}
QTabBar::tab {
    background-color: #1F2937;
    color: #9CA3AF;
    padding: 6px 14px;
    border: 1px solid #374151;
    border-bottom: none;
    border-top-left-radius: 4px;
    border-top-right-radius: 4px;
}
QTabBar::tab:selected {
    background-color: #111827;
    color: #60A5FA;
    font-weight: bold;
}
QProgressBar {
    border: 1px solid #374151;
    border-radius: 4px;
    text-align: center;
    background-color: #1F2937;
    color: white;
}
QProgressBar::chunk {
    background-color: #2563EB;
}
"""


class CFDMainWindow(QMainWindow):
    """Ventana principal de la estación de simulación aerodinámica y aerotérmica."""

    def __init__(self, off_screen: bool = False) -> None:
        super().__init__()
        self.setWindowTitle(
            "Simulador CFD 3D Compresible — Aerodinámica & Aeroheating de Cohetes (OpenFOAM v2312 + PyVista)"
        )
        self.resize(1540, 920)
        self.setStyleSheet(DARK_AEROSPACE_STYLESHEET)

        self.project_root = Path(__file__).resolve().parent.parent.parent
        self.active_case_dir: Path | None = None
        self.runner: OpenFOAMRunner | None = None
        self.benchmark = SupersonicConeBenchmark()
        self.tm_solver = TaylorMaccollSolver()
        self.case_generator = OpenFOAMCaseGenerator()

        self._build_layout(off_screen=off_screen)
        self._connect_signals()

        # Si el caso validado benchmark_cone_m2 ya existe en disco, cargarlo al iniciar
        default_bench = self.project_root / "cases" / "benchmark_cone_m2"
        if (default_bench / "0.003").exists() or len(list(default_bench.glob("0.*"))) > 0:
            self.load_case_directory(default_bench)

    def _build_layout(self, off_screen: bool = False) -> None:
        central = QWidget(self)
        self.setCentralWidget(central)
        root_layout = QVBoxLayout(central)
        root_layout.setContentsMargins(6, 6, 6, 6)

        # Splitter Horizontal Principal: [Panel Izquierdo | Área Derecha (Visor 3D + Gráficas)]
        self.main_splitter = QSplitter(Qt.Orientation.Horizontal)

        self.control_panel = FlightControlPanel(self)
        self.control_panel.setMinimumWidth(320)
        self.control_panel.setMaximumWidth(390)

        self.right_splitter = QSplitter(Qt.Orientation.Vertical)
        self.viewport_panel = Viewport3DPanel(self, off_screen=off_screen)
        self.plots_panel = TelemetryPlotsPanel(self)

        self.right_splitter.addWidget(self.viewport_panel)
        self.right_splitter.addWidget(self.plots_panel)
        self.right_splitter.setStretchFactor(0, 3)
        self.right_splitter.setStretchFactor(1, 2)

        self.main_splitter.addWidget(self.control_panel)
        self.main_splitter.addWidget(self.right_splitter)
        self.main_splitter.setStretchFactor(0, 0)
        self.main_splitter.setStretchFactor(1, 1)

        root_layout.addWidget(self.main_splitter)

        # Barra de Estado Inferior
        status_bar = QStatusBar(self)
        self.setStatusBar(status_bar)

        self.lbl_status = QLabel("Listo. Seleccione condiciones de vuelo o inspeccione el caso en 3D.")
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        self.progress_bar.setFixedWidth(240)

        status_bar.addWidget(self.lbl_status, stretch=1)
        status_bar.addPermanentWidget(self.progress_bar)

    def _connect_signals(self) -> None:
        self.control_panel.start_simulation_requested.connect(self.start_simulation)
        self.control_panel.abort_simulation_requested.connect(self.abort_simulation)
        self.control_panel.load_case_requested.connect(self.load_case_directory)
        self.control_panel.export_aerodb_requested.connect(self.export_current_to_aerodb)

    def start_simulation(self, config: dict) -> None:
        """Prepara el caso en disco y lanza OpenFOAM de forma asíncrona en WSL2."""
        cond = config["flight_condition"]
        mode = config["mode"]
        cone_angle = config["cone_angle_deg"]

        self.plots_panel.clear_telemetry()
        self.progress_bar.setValue(0)

        if mode == "validation":
            case_name = f"cone_m{cond.mach:.2f}_theta{cone_angle:.0f}".replace(".", "p")
            case_dir = self.project_root / "cases" / case_name
            out_dir, _, tm_res = self.benchmark.setup_case(
                case_dir=case_dir,
                mach_inf=cond.mach,
                cone_angle_deg=cone_angle
            )
            self.active_case_dir = out_dir
            self.plots_panel.append_log_line(
                f"=== MODO VALIDACIÓN TAYLOR-MACCOLL (M={cond.mach:.2f}, θc={cone_angle:.1f}°) ===\n"
                f" -> Ángulo de Onda de Choque Analítico (β_TM): {tm_res.shock_angle_deg:.3f}°\n"
                f" -> Salto de Presión Exacto (pc/p∞):           {tm_res.pressure_ratio_surface:.4f}\n"
                f" -> Mach en Superficie Exacto (Mc):            {tm_res.mach_surface:.4f}\n"
            )
            flow_transit = 0.5 / cond.velocity_magnitude_m_s
            end_time = round(4.0 * flow_transit, 5)
            run_blockmesh = True
        else:
            case_name = f"rocket_m{cond.mach:.2f}_h{cond.altitude_m:.0f}_a{cond.alpha_deg:.1f}".replace(".", "p")
            case_dir = self.project_root / "cases" / case_name
            end_time = 0.02
            out_dir = self.case_generator.generate_case(
                output_case_dir=case_dir,
                flight_cond=cond,
                end_time=end_time
            )
            self.active_case_dir = out_dir
            run_blockmesh = (out_dir / "system" / "blockMeshDict").exists()

        self.control_panel.set_running_state(True)
        self.runner = OpenFOAMRunner(
            case_dir_windows=self.active_case_dir,
            run_blockmesh=run_blockmesh,
            end_time=end_time
        )
        self.runner.status_changed.connect(self.lbl_status.setText)
        self.runner.log_line_received.connect(self.plots_panel.append_log_line)
        self.runner.iteration_updated.connect(self.plots_panel.append_iteration)
        self.runner.progress_changed.connect(lambda val: self.progress_bar.setValue(int(val)))
        self.runner.finished_success.connect(self._on_simulation_finished)
        self.runner.error_occurred.connect(self._on_simulation_error)
        self.runner.start()

    def abort_simulation(self) -> None:
        if self.runner and self.runner.isRunning():
            self.runner.cancel()
            self.control_panel.set_running_state(False)
            self.lbl_status.setText("Simulación cancelada por el usuario.")

    def _on_simulation_finished(self, elapsed_s: float) -> None:
        self.control_panel.set_running_state(False)
        self.progress_bar.setValue(100)
        if self.active_case_dir:
            self.load_case_directory(self.active_case_dir)

    def _on_simulation_error(self, err_msg: str) -> None:
        self.control_panel.set_running_state(False)
        self.plots_panel.append_log_line(f"\n[ERROR CRÍTICO]: {err_msg}")
        QMessageBox.critical(self, "Error en Simulación OpenFOAM", err_msg)

    def load_case_directory(self, case_dir: Path | str) -> None:
        """Carga las gráficas históricas y renderiza el campo 3D de un caso en disco."""
        cpath = Path(case_dir)
        if not cpath.exists():
            QMessageBox.warning(self, "Carpeta no encontrada", f"No existe el directorio:\n{cpath}")
            return

        self.active_case_dir = cpath
        self.lbl_status.setText(f"Cargando resultados 3D desde {cpath.name}...")

        # 1. Poblar gráficas de convergencia y fuerzas desde los logs
        self.plots_panel.load_history_from_case(cpath)

        # 2. Calcular ángulo teórico de Taylor-Maccoll si estamos en modo validación supersónica
        cond = self.control_panel.current_condition
        assert cond is not None
        theo_beta: float | None = None
        cone_angle = self.control_panel.spin_cone_angle.value()

        if cond.mach > 1.15:
            try:
                tm_res = self.tm_solver.solve(mach_inf=cond.mach, cone_angle_deg=cone_angle)
                theo_beta = tm_res.shock_angle_deg
            except ValueError:
                theo_beta = None

        # 3. Renderizar escena 3D en el visor PyVista
        try:
            aero_metrics = self.viewport_panel.load_openfoam_case(
                case_dir=cpath,
                t_inf_k=cond.temperature_k,
                p_inf_pa=cond.pressure_pa,
                q_inf_pa=cond.dynamic_pressure_pa,
                mach_inf=cond.mach,
                theo_shock_angle_deg=theo_beta
            )
            self.plots_panel.append_log_line(
                f"\n=== POST-PROCESADO 3D COMPLETADO ({cpath.name}) ===\n"
                f" -> Temperatura Máxima de Pared (CFD):      {aero_metrics.max_wall_temperature_k:.2f} K\n"
                f" -> Temperatura de Remanso Teórica (T0):    {aero_metrics.theoretical_stagnation_temp_k:.2f} K\n"
                f" -> Temp. Recuperación Turbulenta (T_rec):  {aero_metrics.theoretical_recovery_temp_turb_k:.2f} K\n"
                f" -> Coordenada Punto Caliente (x, y, z):    {aero_metrics.hotspot_coordinates_m}\n"
            )

            # Si es un caso de cono, ejecutar evaluación cuantitativa frente a Taylor-Maccoll
            if "cone" in cpath.name.lower() and cond.mach > 1.2:
                report = self.benchmark.evaluate_cfd_results(
                    case_dir=cpath,
                    mach_inf=cond.mach,
                    cone_angle_deg=cone_angle
                )
                self.plots_panel.append_log_line("\n" + report.to_markdown_table())

            self.progress_bar.setValue(100)
            self.lbl_status.setText(
                f"Caso '{cpath.name}' renderizado en GPU | T_max = {aero_metrics.max_wall_temperature_k:.1f} K"
            )

        except Exception as e:
            self.lbl_status.setText(f"Error al cargar malla 3D: {e}")
            self.plots_panel.append_log_line(f"[Error VTK/PyVista]: {e}")

    def export_current_to_aerodb(self) -> None:
        """Exporta el estado estacionario actual a CSV y JSON mediante AeroDBExporter."""
        cond = self.control_panel.current_condition
        aero = self.viewport_panel.aero_metrics
        if cond is None or aero is None:
            QMessageBox.information(
                self,
                "Sin datos 3D cargados",
                "Primero carga o simula un caso para exportar sus métricas a AeroDB."
            )
            return

        cd_last = self.plots_panel.cd_vals[-1] if self.plots_panel.cd_vals else 0.0
        cl_last = self.plots_panel.cl_vals[-1] if self.plots_panel.cl_vals else 0.0
        cm_last = self.plots_panel.cm_vals[-1] if self.plots_panel.cm_vals else 0.0

        x_cp = AeroDBExporter.compute_center_of_pressure(
            cd=cd_last,
            cl=cl_last,
            cm_pitch=cm_last,
            alpha_deg=cond.alpha_deg
        )

        case_label = self.active_case_dir.name if self.active_case_dir else "manual_case"
        record = AeroDBRecord(
            case_name=case_label,
            mach=cond.mach,
            altitude_km=cond.altitude_m / 1000.0,
            alpha_deg=cond.alpha_deg,
            p_inf_pa=cond.pressure_pa,
            t_inf_k=cond.temperature_k,
            u_inf_m_s=cond.velocity_magnitude_m_s,
            cd=cd_last,
            cl=cl_last,
            cm_pitch=cm_last,
            cd_pressure=cd_last,
            cd_friction=0.0,
            x_cp_m=x_cp,
            max_wall_temp_k=aero.max_wall_temperature_k,
            stagnation_temp_theo_k=aero.theoretical_stagnation_temp_k
        )

        results_dir = self.project_root / "results"
        csv_path = AeroDBExporter.append_csv(record, results_dir / "aerodb.csv")
        json_path = AeroDBExporter.export_json(record, results_dir / f"{case_label}_summary.json")

        self.lbl_status.setText(f"Exportado con éxito a {csv_path.name} y {json_path.name}")
        QMessageBox.information(
            self,
            "AeroDB Exportado",
            f"Registro guardado con éxito en:\n- {csv_path}\n- {json_path}"
        )
