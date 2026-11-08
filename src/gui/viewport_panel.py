"""
Panel Central del Visor 3D Acelerado por Hardware (Fase 8.4 y 8.5).
Integra QVTKRenderWindowInteractor (vía pyvistaqt.QtInteractor) dentro de PyQt6.
Permite superponer en tiempo real:
    - Cuerpo 3D Completo (360°) o Corte Seccional Meridiano (180° en Z <= 0) con opacidad regulable.
    - Volumen de Control 3D del túnel numérico (cilindro/caja de dominio + ejes graduados en metros).
    - Advección Cinemática de Partículas de Fluido en tiempo real (60 FPS mediante QTimer).
    - Plano de Efecto Schlieren Numérico simétrico (|∇ρ| y log10(1 + |∇ρ|)).
    - Cono de onda de choque analítico exacto de Taylor-Maccoll (±β_TM).
    - Líneas de corriente (Streamlines) integradas con RK45 y coloreadas por número de Mach.
"""

import math
from pathlib import Path
import numpy as np
import pyvista as pv
from pyvistaqt import QtInteractor
from PyQt6.QtCore import QTimer
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
from src.postprocess.particles import LagrangianParticleSystem
from src.postprocess.schlieren import compute_numerical_schlieren
from src.postprocess.streamlines import generate_streamlines


def reconstruct_cone_from_wedge(
    wedge_surface: pv.PolyData,
    full_360: bool = True,
    n_sectors_full: int = 72
) -> pv.PolyData:
    """
    Si la superficie extraída proviene de una cuña axisimétrica de 5° (benchmark de cono),
    reconstruye por revolución alrededor del eje X:
        - Si full_360=True: el cono 3D completo de 360° (-180° a +180°).
        - Si full_360=False: el semicono posterior de 180° (Z <= 0, abarcando de +Y a -Y)
          para dejar despejado el plano meridiano Z=0 sin perder la silueta superior e inferior.
    Conserva intactos todos los campos escalares físicos (T, p, Cp).
    """
    _, _, ymin, ymax, zmin, zmax = wedge_surface.bounds
    if (zmax - zmin) > 0.25 * max(ymax - ymin, 1e-4):
        return wedge_surface

    sectors: list[pv.PolyData] = []
    if full_360:
        angles = np.linspace(0.0, 360.0, n_sectors_full, endpoint=False)
    else:
        # En la cuña original, ángulo 0° está en +Y (Z=0).
        # Rotando de 0° a -180° (en sentido horario visto desde +X) barremos de +Y a -Y por Z <= 0.
        angles = np.linspace(-180.0, 0.0, max(n_sectors_full // 2, 18), endpoint=True)

    for ang in angles:
        rotated = wedge_surface.copy(deep=True)
        rotated.rotate_x(float(ang), inplace=True)
        sectors.append(rotated)

    return sectors[0].merge(sectors[1:])


def build_control_volume_cylinder(bounds: tuple[float, float, float, float, float, float]) -> pv.PolyData:
    """
    Construye la malla envolvente 3D del volumen de control del túnel de viento numérico
    para dominios axisimétricos o prismáticos.
    """
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    length = xmax - xmin
    center_x = 0.5 * (xmin + xmax)
    radius = max(abs(ymin), abs(ymax), abs(zmin), abs(zmax))

    cyl = pv.Cylinder(
        center=(center_x, 0.0, 0.0),
        direction=(1.0, 0.0, 0.0),
        radius=radius,
        height=length,
        resolution=48,
        capping=False
    )
    return cyl


class Viewport3DPanel(QWidget):
    """Widget central con el motor de renderizado 3D VTK/PyVista, partículas en vivo y volumen de control."""

    def __init__(self, parent: QWidget | None = None, off_screen: bool = False) -> None:
        super().__init__(parent)
        self.off_screen = off_screen

        # Caches de geometría y campos procesados del caso activo
        self.raw_surface_mesh: pv.PolyData | None = None
        self.surface_mesh: pv.PolyData | None = None
        self.schlieren_slice: pv.DataSet | None = None
        self.streamlines_mesh: pv.PolyData | None = None
        self.control_volume_mesh: pv.PolyData | None = None
        self.particle_system: LagrangianParticleSystem | None = None
        self.aero_metrics: AeroheatingMetrics | None = None
        self.theo_shock_angle_deg: float | None = None
        self.domain_length_m: float = 0.5
        self.domain_bounds: tuple[float, float, float, float, float, float] | None = None

        # Temporizador para la animación de partículas a ~33 FPS
        self.anim_timer = QTimer(self)
        self.anim_timer.setInterval(30)
        self.anim_timer.timeout.connect(self._on_animation_tick)

        self._build_ui()

    def _build_ui(self) -> None:
        main_layout = QVBoxLayout(self)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(2)

        # ------------------------------------------------------------------
        # Fila 1 de Barra de Herramientas: Superficie, Geometría 360°, Opacidad y Vistas
        # ------------------------------------------------------------------
        toolbar_row1 = QHBoxLayout()
        toolbar_row1.setContentsMargins(8, 6, 8, 2)
        toolbar_row1.setSpacing(8)

        lbl_field = QLabel("Campo:")
        lbl_field.setStyleSheet("font-weight: bold;")
        self.combo_field = QComboBox()
        self.combo_field.addItems([
            "Temperatura Pared T [K] (Aeroheating)",
            "Coeficiente de Presión Cp [-]",
            "Presión Estática p [Pa]",
            "Sólido Metálico Aeroespacial"
        ])

        lbl_body = QLabel("Geometría:")
        lbl_body.setStyleSheet("font-weight: bold;")
        self.combo_body_mode = QComboBox()
        self.combo_body_mode.addItems([
            "Cuerpo Completo 360°",
            "Corte Seccional 180°"
        ])

        lbl_opacity = QLabel("Opacidad:")
        lbl_opacity.setStyleSheet("font-weight: bold;")
        self.combo_opacity = QComboBox()
        self.combo_opacity.addItems([
            "85% (Semi-translúcido)",
            "100% (Sólido Opaco)",
            "50% (Rayos-X Interior)"
        ])

        btn_iso = QPushButton("🎥 Vista 3D Iso")
        btn_xy = QPushButton("📐 Vista XY (Choque)")
        btn_front = QPushButton("🎯 Vista Frontal")

        toolbar_row1.addWidget(lbl_field)
        toolbar_row1.addWidget(self.combo_field)
        toolbar_row1.addWidget(lbl_body)
        toolbar_row1.addWidget(self.combo_body_mode)
        toolbar_row1.addWidget(lbl_opacity)
        toolbar_row1.addWidget(self.combo_opacity)
        toolbar_row1.addStretch(1)
        toolbar_row1.addWidget(btn_iso)
        toolbar_row1.addWidget(btn_xy)
        toolbar_row1.addWidget(btn_front)

        # ------------------------------------------------------------------
        # Fila 2 de Barra de Herramientas: Capas Físicas, Partículas y Volumen de Control
        # ------------------------------------------------------------------
        toolbar_row2 = QHBoxLayout()
        toolbar_row2.setContentsMargins(8, 2, 8, 4)
        toolbar_row2.setSpacing(12)

        self.chk_particles = QCheckBox("✨ Partículas de Fluido (Animación)")
        self.chk_particles.setChecked(True)
        self.chk_particles.setStyleSheet("color: #38BDF8; font-weight: bold;")

        self.chk_control_vol = QCheckBox("📦 Volumen de Control 3D")
        self.chk_control_vol.setChecked(True)

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

        toolbar_row2.addWidget(self.chk_particles)
        toolbar_row2.addWidget(self.chk_control_vol)
        toolbar_row2.addWidget(self.chk_schlieren)
        toolbar_row2.addWidget(self.chk_theo_shock)
        toolbar_row2.addWidget(self.chk_streamlines)
        toolbar_row2.addWidget(self.chk_edges)
        toolbar_row2.addWidget(self.chk_hotspot)
        toolbar_row2.addStretch(1)

        main_layout.addLayout(toolbar_row1)
        main_layout.addLayout(toolbar_row2)

        # ------------------------------------------------------------------
        # Lienzo Interactivo 3D VTK / PyVista
        # ------------------------------------------------------------------
        if self.off_screen:
            self.plotter = pv.Plotter(off_screen=True)
        else:
            self.plotter = QtInteractor(self)
            main_layout.addWidget(self.plotter.interactor, stretch=1)

        self.plotter.set_background("#0B0F19", top="#1E293B")
        self.plotter.add_axes(color="white")

        # Conectar eventos de ambas barras de herramientas
        self.combo_field.currentIndexChanged.connect(self.refresh_scene)
        self.combo_body_mode.currentIndexChanged.connect(self._on_body_mode_changed)
        self.combo_opacity.currentIndexChanged.connect(self.refresh_scene)

        self.chk_particles.toggled.connect(self._on_particles_toggled)
        self.chk_control_vol.toggled.connect(self.refresh_scene)
        self.chk_schlieren.toggled.connect(self.refresh_scene)
        self.chk_theo_shock.toggled.connect(self.refresh_scene)
        self.chk_streamlines.toggled.connect(self.refresh_scene)
        self.chk_edges.toggled.connect(self.refresh_scene)
        self.chk_hotspot.toggled.connect(self.refresh_scene)

        btn_xy.clicked.connect(self.set_view_xy)
        btn_iso.clicked.connect(self.set_view_isometric)
        btn_front.clicked.connect(self.set_view_front)

    def _get_selected_opacity(self) -> float:
        idx = self.combo_opacity.currentIndex()
        if idx == 0:
            return 0.85
        if idx == 1:
            return 1.0
        return 0.50

    def _on_body_mode_changed(self) -> None:
        if self.raw_surface_mesh is not None:
            full_360 = self.combo_body_mode.currentIndex() == 0
            self.surface_mesh = reconstruct_cone_from_wedge(self.raw_surface_mesh, full_360=full_360)
            self.refresh_scene()

    def _on_particles_toggled(self, checked: bool) -> None:
        if checked and not self.off_screen and self.particle_system is not None:
            self.anim_timer.start()
        else:
            self.anim_timer.stop()
        self.refresh_scene()

    def _on_animation_tick(self) -> None:
        """Actualiza in-place las posiciones 3D de las partículas en la GPU sin reconstruir la escena."""
        if self.particle_system is not None and self.chk_particles.isChecked():
            self.particle_system.step(delta_phase=0.012)
            self.plotter.render()

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
        ejecuta los algoritmos de post-procesado (Aeroheating, Schlieren, Streamlines,
        Volumen de Control 3D y Sistema Lagrangiano de Partículas) y renderiza en la GPU.
        """
        self.anim_timer.stop()

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
        xmin, xmax, ymin, ymax, zmin, zmax = internal_mesh.bounds
        self.domain_length_m = max(xmax, xmax - xmin)
        self.theo_shock_angle_deg = theo_shock_angle_deg

        # 1. Extraer superficie térmica y de presión del vehículo y reconstruir 360° / 180°
        self.raw_surface_mesh, self.aero_metrics = extract_aeroheating_surface(
            multiblock_mesh=multiblock,
            t_inf_k=t_inf_k,
            p_inf_pa=p_inf_pa,
            q_inf_pa=q_inf_pa,
            mach_inf=mach_inf
        )
        full_360 = self.combo_body_mode.currentIndex() == 0
        self.surface_mesh = reconstruct_cone_from_wedge(self.raw_surface_mesh, full_360=full_360)

        # 2. Calcular plano de efecto Schlieren numérico simétrico (|∇ρ|)
        self.schlieren_slice = compute_numerical_schlieren(
            internal_mesh=internal_mesh,
            slice_normal=(0.0, 0.0, 1.0),
            slice_origin=(0.0, 0.0, 0.0),
            mirror_symmetry=True
        )

        # 3. Construir el Volumen de Control 3D del túnel numérico
        rmax = max(abs(ymin), abs(ymax))
        self.domain_bounds = (xmin, xmax, -rmax, rmax, -rmax, rmax)
        self.control_volume_mesh = build_control_volume_cylinder(self.domain_bounds)

        # 4. Generar líneas de corriente 3D y Sistema Lagrangiano de Partículas
        self.streamlines_mesh = generate_streamlines(
            internal_mesh=internal_mesh,
            n_seeds=18,
            mirror_symmetry=True
        )
        self.particle_system = LagrangianParticleSystem(
            streamlines_polydata=self.streamlines_mesh,
            n_particles_per_track=8,
            revolve_3d_azimuths=6
        )

        self.refresh_scene()
        self.set_view_isometric()

        if self.chk_particles.isChecked() and not self.off_screen:
            self.anim_timer.start()

        return self.aero_metrics

    def refresh_scene(self) -> None:
        """Reconstruye los actores 3D en el renderizador VTK según las capas activas."""
        self.plotter.clear()
        self.plotter.add_axes(color="white")

        show_edges = self.chk_edges.isChecked()
        body_opacity = self._get_selected_opacity()

        # 1. Volumen de Control 3D (Túnel de Viento Numérico y Caja Graduada)
        if self.chk_control_vol.isChecked() and self.control_volume_mesh is not None:
            self.plotter.add_mesh(
                self.control_volume_mesh,
                style="wireframe",
                color="#38BDF8",
                opacity=0.22,
                line_width=1.0
            )
            self.plotter.show_bounds(
                mesh=self.control_volume_mesh,
                grid="back",
                location="outer",
                color="#94A3B8",
                font_size=9,
                xtitle="X Axial [m]",
                ytitle="Y Vertical [m]",
                ztitle="Z Lateral [m]"
            )

        # 2. Plano de Efecto Schlieren Numérico (|∇ρ|)
        if self.chk_schlieren.isChecked() and self.schlieren_slice is not None:
            self.plotter.add_mesh(
                self.schlieren_slice,
                scalars="schlieren_log",
                cmap="gist_yarg",
                opacity=0.88,
                show_edges=show_edges,
                edge_color="#64748B",
                scalar_bar_args={
                    "title": "Schlieren log10(1+|∇ρ|)",
                    "color": "white",
                    "title_font_size": 11,
                    "label_font_size": 9,
                    "position_x": 0.14,
                    "position_y": 0.02,
                    "width": 0.22,
                    "height": 0.045,
                }
            )

        # 3. Superficie Sólida 3D del Vehículo (Cono 360° / Cohete)
        if self.surface_mesh is not None:
            field_idx = self.combo_field.currentIndex()
            if field_idx == 0 and "T" in self.surface_mesh.point_data:
                self.plotter.add_mesh(
                    self.surface_mesh,
                    scalars="T",
                    cmap="inferno",
                    opacity=body_opacity,
                    show_edges=show_edges,
                    specular=0.45,
                    specular_power=30,
                    scalar_bar_args={
                        "title": "Temp. Pared T [K]",
                        "color": "white",
                        "title_font_size": 11,
                        "label_font_size": 9,
                        "position_x": 0.42,
                        "position_y": 0.02,
                        "width": 0.22,
                        "height": 0.045,
                    }
                )
            elif field_idx == 1 and "Cp" in self.surface_mesh.point_data:
                self.plotter.add_mesh(
                    self.surface_mesh,
                    scalars="Cp",
                    cmap="coolwarm",
                    opacity=body_opacity,
                    show_edges=show_edges,
                    specular=0.45,
                    scalar_bar_args={
                        "title": "Coef. Presión Cp [-]",
                        "color": "white",
                        "title_font_size": 11,
                        "label_font_size": 9,
                        "position_x": 0.42,
                        "position_y": 0.02,
                        "width": 0.22,
                        "height": 0.045,
                    }
                )
            elif field_idx == 2 and "p" in self.surface_mesh.point_data:
                self.plotter.add_mesh(
                    self.surface_mesh,
                    scalars="p",
                    cmap="viridis",
                    opacity=body_opacity,
                    show_edges=show_edges,
                    specular=0.45,
                    scalar_bar_args={
                        "title": "Presión p [Pa]",
                        "color": "white",
                        "title_font_size": 11,
                        "label_font_size": 9,
                        "position_x": 0.42,
                        "position_y": 0.02,
                        "width": 0.22,
                        "height": 0.045,
                    }
                )
            else:
                self.plotter.add_mesh(
                    self.surface_mesh,
                    color="#94A3B8",
                    opacity=body_opacity,
                    show_edges=show_edges,
                    specular=0.6,
                    specular_power=40
                )

        # 4. Generatrices Analíticas Superior e Inferior de Taylor-Maccoll (±β_TM)
        if (
            self.chk_theo_shock.isChecked()
            and self.theo_shock_angle_deg is not None
            and self.theo_shock_angle_deg > 0.0
        ):
            beta_rad = math.radians(self.theo_shock_angle_deg)
            x_end = self.domain_length_m
            y_end = x_end * math.tan(beta_rad)
            shock_top = pv.Line((0.0, 0.0, 0.002), (x_end, y_end, 0.002))
            shock_bot = pv.Line((0.0, 0.0, 0.002), (x_end, -y_end, 0.002))
            self.plotter.add_mesh(shock_top, color="#10B981", line_width=3.5)
            self.plotter.add_mesh(shock_bot, color="#10B981", line_width=3.5)

        # 5. Líneas de Corriente Estáticas (Streamlines coloreadas por Mach)
        if (
            self.chk_streamlines.isChecked()
            and self.streamlines_mesh is not None
            and self.streamlines_mesh.n_points > 0
        ):
            self.plotter.add_mesh(
                self.streamlines_mesh,
                scalars="Mach",
                cmap="plasma",
                line_width=1.8,
                opacity=0.65,
                scalar_bar_args={
                    "title": "Mach Local [-]",
                    "color": "white",
                    "title_font_size": 11,
                    "label_font_size": 9,
                    "position_x": 0.70,
                    "position_y": 0.02,
                    "width": 0.22,
                    "height": 0.045,
                }
            )

        # 6. Nube 3D de Partículas de Fluido en Movimiento (Lagrangian Tracers)
        if (
            self.chk_particles.isChecked()
            and self.particle_system is not None
            and self.particle_system.cloud_polydata.n_points > 1
        ):
            self.plotter.add_mesh(
                self.particle_system.cloud_polydata,
                scalars="Mach",
                cmap="plasma",
                render_points_as_spheres=True,
                point_size=8.5,
                show_scalar_bar=not self.chk_streamlines.isChecked(),
                scalar_bar_args={
                    "title": "Mach Local [-]",
                    "color": "white",
                    "title_font_size": 11,
                    "label_font_size": 9,
                    "position_x": 0.70,
                    "position_y": 0.02,
                    "width": 0.22,
                    "height": 0.045,
                }
            )

        # 7. Etiqueta 3D del Punto Caliente / Remanso
        if self.chk_hotspot.isChecked() and self.aero_metrics is not None:
            pt = np.array([self.aero_metrics.hotspot_coordinates_m])
            label_str = f"T_max = {self.aero_metrics.max_wall_temperature_k:.1f} K"
            self.plotter.add_point_labels(
                pt,
                [label_str],
                point_size=10,
                point_color="#EF4444",
                text_color="white",
                font_size=11,
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
        """Coloca la cámara en perspectiva isométrica 3D aeroespacial."""
        self.plotter.view_isometric()
        self.plotter.reset_camera()
        self.plotter.render()

    def set_view_front(self) -> None:
        """Coloca la cámara mirando de frente a la proa del vehículo (plano Y-Z)."""
        self.plotter.view_yz()
        self.plotter.reset_camera()
        self.plotter.render()
