"""
Panel de Control de Vuelo y Telemetría Atmosférica ISA 1976 (Fase 8.2).
Permite configurar (Mach, Altitud, Ángulo de Ataque), alternar entre el Modo de Validación
Analítica (Cono Taylor-Maccoll) y Simulación de Cohete 3D, e inspeccionar las variables
termodinámicas de entrada en tiempo real antes de lanzar OpenFOAM.
"""

from pathlib import Path
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from src.physics.atmosphere import AtmosphereISA, FlightCondition
from src.postprocess.aeroheating import compute_theoretical_thermal_limits


class FlightControlPanel(QWidget):
    """Barra lateral izquierda para el control de parámetros físicos y ejecución."""

    # Señales emitidas hacia la ventana principal
    flight_condition_changed = pyqtSignal(object)   # Emite FlightCondition actualizado
    start_simulation_requested = pyqtSignal(dict)   # Emite diccionario de configuración
    abort_simulation_requested = pyqtSignal()       # Solicita detener el runner WSL2
    load_case_requested = pyqtSignal(str)           # Solicita abrir un directorio de caso OpenFOAM
    export_aerodb_requested = pyqtSignal()          # Solicita exportar resultados a CSV/JSON

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.isa = AtmosphereISA()
        self.stl_files: list[Path] = []
        self.current_condition: FlightCondition | None = None

        self._build_ui()
        self._connect_signals()
        self._update_isa_preview()

    def _build_ui(self) -> None:
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(10, 10, 10, 10)
        main_layout.setSpacing(10)

        # 1. Selector de Modo de Simulación
        mode_group = QGroupBox("1. Modo de Operación")
        mode_layout = QVBoxLayout(mode_group)

        self.combo_mode = QComboBox()
        self.combo_mode.addItems([
            "Modo Validación (Cono Taylor-Maccoll)",
            "Simulación Completa Cohete 3D"
        ])
        mode_layout.addWidget(self.combo_mode)

        self.btn_import_stl = QPushButton("📐 Importar STLs del Cohete...")
        self.btn_import_stl.setEnabled(False)
        self.lbl_stl_status = QLabel("Geometría: Cuña Axisimétrica (Cono 15°)")
        self.lbl_stl_status.setStyleSheet("color: #9CA3AF; font-size: 11px;")
        mode_layout.addWidget(self.btn_import_stl)
        mode_layout.addWidget(self.lbl_stl_status)
        main_layout.addWidget(mode_group)

        # 2. Parámetros Cinemáticos de Vuelo
        flight_group = QGroupBox("2. Condiciones de Vuelo")
        flight_form = QFormLayout(flight_group)

        self.spin_mach = QDoubleSpinBox()
        self.spin_mach.setRange(0.20, 4.00)
        self.spin_mach.setSingleStep(0.10)
        self.spin_mach.setDecimals(2)
        self.spin_mach.setValue(2.00)

        self.spin_altitude = QDoubleSpinBox()
        self.spin_altitude.setRange(0.0, 30.0)
        self.spin_altitude.setSingleStep(0.5)
        self.spin_altitude.setDecimals(1)
        self.spin_altitude.setSuffix(" km")
        self.spin_altitude.setValue(0.0)

        self.spin_alpha = QDoubleSpinBox()
        self.spin_alpha.setRange(0.0, 15.0)
        self.spin_alpha.setSingleStep(0.5)
        self.spin_alpha.setDecimals(1)
        self.spin_alpha.setSuffix(" °")
        self.spin_alpha.setValue(0.0)

        self.spin_cone_angle = QDoubleSpinBox()
        self.spin_cone_angle.setRange(5.0, 35.0)
        self.spin_cone_angle.setSingleStep(1.0)
        self.spin_cone_angle.setDecimals(1)
        self.spin_cone_angle.setSuffix(" °")
        self.spin_cone_angle.setValue(15.0)

        flight_form.addRow("Número de Mach (M∞):", self.spin_mach)
        flight_form.addRow("Altitud ISA (h):", self.spin_altitude)
        flight_form.addRow("Ángulo de Ataque (α):", self.spin_alpha)
        flight_form.addRow("Semiángulo Cono (θc):", self.spin_cone_angle)
        main_layout.addWidget(flight_group)

        # 3. Monitor Termodinámico ISA 1976 en Tiempo Real
        isa_group = QGroupBox("3. Atmósfera ISA 1976 & Aeroheating")
        isa_form = QFormLayout(isa_group)

        self.lbl_temp = QLabel("-")
        self.lbl_press = QLabel("-")
        self.lbl_rho = QLabel("-")
        self.lbl_sound = QLabel("-")
        self.lbl_vel = QLabel("-")
        self.lbl_qinf = QLabel("-")
        self.lbl_tstag = QLabel("-")
        self.lbl_trec = QLabel("-")

        for lbl in [
            self.lbl_temp, self.lbl_press, self.lbl_rho, self.lbl_sound,
            self.lbl_vel, self.lbl_qinf, self.lbl_tstag, self.lbl_trec
        ]:
            lbl.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            lbl.setStyleSheet("font-family: Consolas, monospace; font-weight: bold; color: #60A5FA;")

        self.lbl_tstag.setStyleSheet("font-family: Consolas, monospace; font-weight: bold; color: #F87171;")
        self.lbl_trec.setStyleSheet("font-family: Consolas, monospace; font-weight: bold; color: #FBBF24;")

        isa_form.addRow("Temp. Estática (T∞):", self.lbl_temp)
        isa_form.addRow("Presión Estática (p∞):", self.lbl_press)
        isa_form.addRow("Densidad (ρ∞):", self.lbl_rho)
        isa_form.addRow("Vel. Sonido (a∞):", self.lbl_sound)
        isa_form.addRow("Velocidad Flujo (U∞):", self.lbl_vel)
        isa_form.addRow("Presión Dinámica (q∞):", self.lbl_qinf)
        isa_form.addRow("Temp. Remanso (T0):", self.lbl_tstag)
        isa_form.addRow("Temp. Recup. Turb (Trec):", self.lbl_trec)
        main_layout.addWidget(isa_group)

        # 4. Controles de Ejecución y Carga
        exec_group = QGroupBox("4. Control del Solver (WSL2)")
        exec_layout = QVBoxLayout(exec_group)

        self.btn_run = QPushButton("🚀 Iniciar Simulación CFD")
        self.btn_run.setStyleSheet(
            "background-color: #2563EB; color: white; font-weight: bold; padding: 8px; border-radius: 4px;"
        )

        self.btn_abort = QPushButton("🛑 Abortar Simulación")
        self.btn_abort.setEnabled(False)
        self.btn_abort.setStyleSheet(
            "background-color: #DC2626; color: white; font-weight: bold; padding: 6px; border-radius: 4px;"
        )

        self.btn_load_benchmark = QPushButton("📂 Cargar Caso Cono M=2.0 Validado")
        self.btn_load_custom = QPushButton("📁 Abrir Otra Carpeta de Caso...")
        self.btn_export_db = QPushButton("💾 Exportar Punto a AeroDB (CSV/JSON)")

        exec_layout.addWidget(self.btn_run)
        exec_layout.addWidget(self.btn_abort)
        exec_layout.addWidget(self.btn_load_benchmark)
        exec_layout.addWidget(self.btn_load_custom)
        exec_layout.addWidget(self.btn_export_db)
        main_layout.addWidget(exec_group)

        main_layout.addStretch(1)

    def _connect_signals(self) -> None:
        self.combo_mode.currentIndexChanged.connect(self._on_mode_changed)
        self.spin_mach.valueChanged.connect(self._update_isa_preview)
        self.spin_altitude.valueChanged.connect(self._update_isa_preview)
        self.spin_alpha.valueChanged.connect(self._update_isa_preview)
        self.spin_cone_angle.valueChanged.connect(self._update_isa_preview)

        self.btn_import_stl.clicked.connect(self._on_import_stl)
        self.btn_run.clicked.connect(self._on_start_clicked)
        self.btn_abort.clicked.connect(self.abort_simulation_requested.emit)
        self.btn_load_benchmark.clicked.connect(self._on_load_default_benchmark)
        self.btn_load_custom.clicked.connect(self._on_load_custom_case)
        self.btn_export_db.clicked.connect(self.export_aerodb_requested.emit)

    def _on_mode_changed(self, index: int) -> None:
        is_validation = index == 0
        self.btn_import_stl.setEnabled(not is_validation)
        self.spin_cone_angle.setEnabled(is_validation)
        if is_validation:
            self.spin_alpha.setValue(0.0)
            self.spin_alpha.setEnabled(False)
            self.lbl_stl_status.setText(
                f"Geometría: Cuña Axisimétrica (Cono {self.spin_cone_angle.value():.1f}°)"
            )
        else:
            self.spin_alpha.setEnabled(True)
            if self.stl_files:
                self.lbl_stl_status.setText(f"{len(self.stl_files)} archivos STL cargados")
            else:
                self.lbl_stl_status.setText("Pendiente importar STLs (nose, fuselage, fins)")
        self._update_isa_preview()

    def _update_isa_preview(self) -> None:
        mach = self.spin_mach.value()
        alt_km = self.spin_altitude.value()
        alpha = self.spin_alpha.value()

        cond = self.isa.compute_flight_condition(mach=mach, altitude_km=alt_km, alpha_deg=alpha)
        self.current_condition = cond

        t_stag, t_rec_turb, _ = compute_theoretical_thermal_limits(
            t_inf_k=cond.temperature_k,
            mach_inf=cond.mach
        )

        self.lbl_temp.setText(f"{cond.temperature_k:.2f} K ({cond.temperature_k - 273.15:.1f} °C)")
        self.lbl_press.setText(f"{cond.pressure_pa:,.1f} Pa")
        self.lbl_rho.setText(f"{cond.density_kg_m3:.4f} kg/m³")
        self.lbl_sound.setText(f"{cond.speed_of_sound_m_s:.1f} m/s")
        self.lbl_vel.setText(f"{cond.velocity_magnitude_m_s:.1f} m/s")
        self.lbl_qinf.setText(f"{cond.dynamic_pressure_pa / 1000.0:.2f} kPa")
        self.lbl_tstag.setText(f"{t_stag:.1f} K ({t_stag - 273.15:.0f} °C)")
        self.lbl_trec.setText(f"{t_rec_turb:.1f} K ({t_rec_turb - 273.15:.0f} °C)")

        if self.combo_mode.currentIndex() == 0:
            self.lbl_stl_status.setText(
                f"Geometría: Cuña Axisimétrica (Cono {self.spin_cone_angle.value():.1f}°)"
            )

        self.flight_condition_changed.emit(cond)

    def _on_import_stl(self) -> None:
        files, _ = QFileDialog.getOpenFileNames(
            self,
            "Seleccionar parches STL del Cohete (nose, fuselage, fins, base)",
            str(Path.cwd()),
            "Archivos STL (*.stl)"
        )
        if files:
            self.stl_files = [Path(f) for f in files]
            names = ", ".join(p.stem for p in self.stl_files)
            self.lbl_stl_status.setText(f"STLs ({len(self.stl_files)}): {names}")

    def _on_start_clicked(self) -> None:
        if self.current_condition is None:
            self._update_isa_preview()

        config = {
            "mode": "validation" if self.combo_mode.currentIndex() == 0 else "rocket",
            "flight_condition": self.current_condition,
            "cone_angle_deg": self.spin_cone_angle.value(),
            "stl_files": self.stl_files
        }
        self.start_simulation_requested.emit(config)

    def _on_load_default_benchmark(self) -> None:
        default_dir = Path(__file__).resolve().parent.parent.parent / "cases" / "benchmark_cone_m2"
        self.load_case_requested.emit(str(default_dir))

    def _on_load_custom_case(self) -> None:
        folder = QFileDialog.getExistingDirectory(
            self,
            "Seleccionar Directorio de Caso OpenFOAM",
            str(Path(__file__).resolve().parent.parent.parent / "cases")
        )
        if folder:
            self.load_case_requested.emit(folder)

    def set_running_state(self, is_running: bool) -> None:
        """Habilita o deshabilita botones mientras el solver está corriendo en WSL2."""
        self.btn_run.setEnabled(not is_running)
        self.btn_abort.setEnabled(is_running)
        self.spin_mach.setEnabled(not is_running)
        self.spin_altitude.setEnabled(not is_running)
        self.combo_mode.setEnabled(not is_running)
