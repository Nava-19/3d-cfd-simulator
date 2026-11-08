"""
Panel Central del Visor 3D Acelerado por Hardware (Fase 8.4).
Integra QVTKRenderWindowInteractor (vía pyvistaqt.QtInteractor) dentro de PyQt6.
Permite superponer en tiempo real:
    - Superficie 3D del vehículo (Temperatura de Pared T, Coeficiente de Presión Cp, Presión p).
    - Plano de Efecto Schlieren Numérico (|∇ρ| y log10(1 + |∇ρ|)) en escala de grises óptica.
    - Generatriz analítica exacta de la onda de choque oblicua de Taylor-Maccoll (β_TM).
    - Líneas de corriente (Streamlines) integradas con RK45 y coloreadas por número de Mach.
    - Marcador espacial 3D del punto de máximo calentamiento aerodinámico (Hotspot / Remanso).
"""

import math
from pathlib import Path
import numpy as np
import pyvista as pv
from pyvistaqt import QtInteractor
from PyQt6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from src.postprocess.aeroheating import AeroheatingMetrics, extract_aeroheating_surface
from src.postprocess.schlieren import compute_numerical_schlieren
from src.postprocess.streamlines import generate_streamlines


def reconstruct_half_cone_from_wedge(wedge_surface: pv.PolyData, n_sectors: int = 36) -> pv.PolyData:
    """
    Si la superficie extraída proviene de una cuña axisimétrica de 5° (benchmark de cono),
    reconstruye un semicono 3D de 180° (en el semiespacio Z <= 0) rotando copias alrededor
    del eje X. Conserva intactos todos los campos escalares físicos (T, p, Cp).
    """
    _, _, ymin, ymax, zmin, zmax = wedge_surface.bounds
    if (zmax - zmin) > 0.25 * max(ymax - ymin, 1e-4):
        return wedge_surface

    sectors: list[pv.PolyData] = []
    angles = np.linspace(-90.0, 90.0, n_sectors)
    for ang in angles:
        rotated = wedge_surface.copy(deep=True)
        rotated.rotate_x(float(ang), inplace=True)
        sectors.append(rotated)

    merged = sectors[0].merge(sectors[1:])
    return merged


class Viewport3DPanel(QWidget):
    """Widget central con el motor de renderizado 3D VTK/PyVista y barra de capas físicas."""

    def __init__(self, parent: QWidget | None = None, off_screen: bool = False) -> None:
        super().__init__(parent)
        self.off_screen = off_screen

        # Caches de geometría y campos procesados del caso activo
        self.surface_mesh: pv.PolyData | None = None
        self.schlieren_slice: pv.DataSet | None = None
        self.streamlines_mesh: pv.PolyData | None = None
        self.aero_metrics: AeroheatingMetrics | None = None
        self.theo_shock_angle_deg: float | None = None
        self.domain_length_m: float = 0.5

        self._build_ui()

    def _build_ui(self) -> None:
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(4)

        # ------------------------------------------------------------------
        # Barra de Herramientas Superior: Capas Físicas y Vistas de Cámara
        # ------------------------------------------------------------------
        toolbar_layout = QHBoxLayout()
        toolbar_layout.setContentsMargins(8, 6, 8, 4)
        toolbar_layout.setSpacing(10)

        lbl_field = QLabel("Campo Superficie:")
        lbl_field.setStyleSheet("font-weight: bold;")
        self.combo_field = QComboBox()
        self.combo_field.addItems([
            "Temperatura Pared T [K] (Aeroheating)",
            "Coeficiente de Presión Cp [-]",
            "Presión Estática p [Pa]",
            "Sólido Metálico Aeroespacial"
        ])

        self.chk_schlieren = QCheckBox("Plano Schlieren (|∇ρ|)")
        self.chk_schlieren.setChecked(True)

        self.chk_theo_shock = QCheckBox("Choque Teórico Taylor-Maccoll")
        self.chk_theo_shock.setChecked(True)

        self.chk_streamlines = QCheckBox("Streamlines (Mach)")
        self.chk_streamlines.setChecked(True)

        self.chk_edges = QCheckBox("Ver Malla")
        self.chk_edges.setChecked(False)

        self.chk_hotspot = QCheckBox("Etiqueta T_max")
        self.chk_hotspot.setChecked(True)

        btn_xy = QPushButton("Vista XY (Choque)")
        btn_iso = QPushButton("Vista 3D Iso")
        btn_front = QPushButton("Vista Frontal")

        toolbar_layout.addWidget(lbl_field)
        toolbar_layout.addWidget(self.combo_field)
        toolbar_layout.addWidget(self.chk_schlieren)
        toolbar_layout.addWidget(self.chk_theo_shock)
        toolbar_layout.addWidget(self.chk_streamlines)
        toolbar_layout.addWidget(self.chk_edges)
        toolbar_layout.addWidget(self.chk_hotspot)
        toolbar_layout.addStretch(1)
        toolbar_layout.addWidget(btn_xy)
        toolbar_layout.addWidget(btn_iso)
        toolbar_layout.addWidget(btn_front)

        main_layout.addLayout(toolbar_layout)

        # ------------------------------------------------------------------
        # Lienzo Interactivo 3D VTK / PyVista
        # ------------------------------------------------------------------
        if self.off_screen:
            self.plotter = pv.Plotter(off_screen=True)
        else:
            self.plotter = QtInteractor(self)
            main_layout.addWidget(self.plotter.interactor, stretch=1)

        self.plotter.set_background("#0F172A", top="#1E293B")
        self.plotter.add_axes(color="white")

        # Conectar eventos de la barra de herramientas
        self.combo_field.currentIndexChanged.connect(self.refresh_scene)
        self.chk_schlieren.toggled.connect(self.refresh_scene)
        self.chk_theo_shock.toggled.connect(self.refresh_scene)
        self.chk_streamlines.toggled.connect(self.refresh_scene)
        self.chk_edges.toggled.connect(self.refresh_scene)
        self.chk_hotspot.toggled.connect(self.refresh_scene)

        btn_xy.clicked.connect(self.set_view_xy)
        btn_iso.clicked.connect(self.set_view_isometric)
        btn_front.clicked.connect(self.set_view_front)

    def load_openfoam_case(
        self,
        case_dir: Path | str,
        t_inf_k: float = 288.15,
        p_inf_pa: float = 101325.0,
        q_inf_pa: float = 283710.0,
        mach_inf: float = 2.0,
        theo_shock_angle_deg: float | None = None
    ) -> AeroheatingMetrics:
        """
        Lee el último directorio temporal de un caso de OpenFOAM mediante POpenFOAMReader,
        ejecuta los algoritmos de post-procesado (Aeroheating, Schlieren y Streamlines)
        y actualiza la escena 3D en la GPU.
        """
        case_path = Path(case_dir)
        foam_file = case_path / "case.foam"
        if not foam_file.exists():
            foam_file.touch()

        reader = pv.POpenFOAMReader(str(foam_file))
        if len(reader.time_values) == 0:
            raise RuntimeError(f"No se encontraron pasos de tiempo en {case_path}")

        reader.set_active_time_value(reader.time_values[-1])
        multiblock = reader.read()

        internal_mesh = multiblock["internalMesh"]
        xmin, xmax, _, _, _, _ = internal_mesh.bounds
        self.domain_length_m = max(xmax, xmax - xmin)
        self.theo_shock_angle_deg = theo_shock_angle_deg

        # 1. Extraer superficie térmica y de presión del vehículo
        raw_surface, self.aero_metrics = extract_aeroheating_surface(
            multiblock_mesh=multiblock,
            t_inf_k=t_inf_k,
            p_inf_pa=p_inf_pa,
            q_inf_pa=q_inf_pa,
            mach_inf=mach_inf
        )
        self.surface_mesh = reconstruct_half_cone_from_wedge(raw_surface)

        # 2. Calcular plano de efecto Schlieren numérico (|∇ρ|)
        self.schlieren_slice = compute_numerical_schlieren(
            internal_mesh=internal_mesh,
            slice_normal=(0.0, 0.0, 1.0),
            slice_origin=(0.0, 0.0, 0.0)
        )

        # 3. Generar líneas de corriente 3D coloreadas por Mach
        self.streamlines_mesh = generate_streamlines(
            internal_mesh=internal_mesh,
            n_seeds=22
        )

        self.refresh_scene()
        self.set_view_xy()
        return self.aero_metrics

    def refresh_scene(self) -> None:
        """Reconstruye los actores 3D en el renderizador VTK según las capas activas."""
        self.plotter.clear()
        self.plotter.add_axes(color="white")

        show_edges = self.chk_edges.isChecked()

        # 1. Plano de Efecto Schlieren Numérico (|∇ρ|)
        if self.chk_schlieren.isChecked() and self.schlieren_slice is not None:
            self.plotter.add_mesh(
                self.schlieren_slice,
                scalars="schlieren_log",
                cmap="gist_yarg",
                opacity=0.92,
                show_edges=show_edges,
                edge_color="#64748B",
                scalar_bar_args={
                    "title": "Schlieren log10(1 + |∇ρ|)",
                    "color": "white",
                    "position_x": 0.04,
                    "position_y": 0.05,
                    "width": 0.25,
                    "height": 0.06,
                }
            )

        # 2. Superficie Sólida del Cohete / Cono
        if self.surface_mesh is not None:
            field_idx = self.combo_field.currentIndex()
            if field_idx == 0 and "T" in self.surface_mesh.point_data:
                self.plotter.add_mesh(
                    self.surface_mesh,
                    scalars="T",
                    cmap="inferno",
                    show_edges=show_edges,
                    specular=0.35,
                    specular_power=25,
                    scalar_bar_args={
                        "title": "Temp. Pared T [K]",
                        "color": "white",
                        "position_x": 0.36,
                        "position_y": 0.05,
                        "width": 0.25,
                        "height": 0.06,
                    }
                )
            elif field_idx == 1 and "Cp" in self.surface_mesh.point_data:
                self.plotter.add_mesh(
                    self.surface_mesh,
                    scalars="Cp",
                    cmap="coolwarm",
                    show_edges=show_edges,
                    specular=0.35,
                    scalar_bar_args={
                        "title": "Coef. Presión Cp [-]",
                        "color": "white",
                        "position_x": 0.36,
                        "position_y": 0.05,
                        "width": 0.25,
                        "height": 0.06,
                    }
                )
            elif field_idx == 2 and "p" in self.surface_mesh.point_data:
                self.plotter.add_mesh(
                    self.surface_mesh,
                    scalars="p",
                    cmap="viridis",
                    show_edges=show_edges,
                    specular=0.35,
                    scalar_bar_args={
                        "title": "Presión p [Pa]",
                        "color": "white",
                        "position_x": 0.36,
                        "position_y": 0.05,
                        "width": 0.25,
                        "height": 0.06,
                    }
                )
            else:
                self.plotter.add_mesh(
                    self.surface_mesh,
                    color="#94A3B8",
                    show_edges=show_edges,
                    specular=0.6,
                    specular_power=40
                )

        # 3. Generatriz Analítica de Onda de Choque de Taylor-Maccoll (β_TM)
        if (
            self.chk_theo_shock.isChecked()
            and self.theo_shock_angle_deg is not None
            and self.theo_shock_angle_deg > 0.0
        ):
            beta_rad = math.radians(self.theo_shock_angle_deg)
            x_end = self.domain_length_m
            y_end = x_end * math.tan(beta_rad)
            shock_line = pv.Line((0.0, 0.0, 0.001), (x_end, y_end, 0.001))
            self.plotter.add_mesh(
                shock_line,
                color="#10B981",
                line_width=3.5,
                label=f"Choque Exacto Taylor-Maccoll (β={self.theo_shock_angle_deg:.2f}°)"
            )

        # 4. Líneas de Corriente (Streamlines coloreadas por Mach)
        if (
            self.chk_streamlines.isChecked()
            and self.streamlines_mesh is not None
            and self.streamlines_mesh.n_points > 0
        ):
            self.plotter.add_mesh(
                self.streamlines_mesh,
                scalars="Mach",
                cmap="plasma",
                line_width=2.2,
                scalar_bar_args={
                    "title": "Mach Local [-]",
                    "color": "white",
                    "position_x": 0.68,
                    "position_y": 0.05,
                    "width": 0.25,
                    "height": 0.06,
                }
            )

        # 5. Etiqueta 3D del Punto Caliente / Remanso
        if self.chk_hotspot.isChecked() and self.aero_metrics is not None:
            pt = np.array([self.aero_metrics.hotspot_coordinates_m])
            label_str = f"T_max = {self.aero_metrics.max_wall_temperature_k:.1f} K"
            self.plotter.add_point_labels(
                pt,
                [label_str],
                point_size=10,
                point_color="#EF4444",
                text_color="white",
                font_size=12,
                shape_color="#1E293B",
                shape_opacity=0.85,
                always_visible=True
            )

        self.plotter.render()

    def set_view_xy(self) -> None:
        """Coloca la cámara perpendicular al plano meridiano X-Y."""
        self.plotter.view_xy()
        self.plotter.reset_camera()
        self.plotter.render()

    def set_view_isometric(self) -> None:
        """Coloca la cámara en perspectiva isométrica 3D."""
        self.plotter.view_isometric()
        self.plotter.reset_camera()
        self.plotter.render()

    def set_view_front(self) -> None:
        """Coloca la cámara mirando de frente a la proa del vehículo (plano Y-Z)."""
        self.plotter.view_yz()
        self.plotter.reset_camera()
        self.plotter.render()
