"""
Generador automatizado de casos de OpenFOAM para rhoCentralFoam.
Utiliza Jinja2 para inyectar condiciones de vuelo calculadas por AtmosphereISA
y parámetros geométricos en los diccionarios de OpenFOAM:
    1. Cono supersónico axisimétrico (Validación Taylor-Maccoll).
    2. Esfera supersónica axisimétrica (Cuerpo romo con Bow Shock desprendido y estela).
    3. Cohete lanzador 3D completo con parches STL y dominio de estela (snappyHexMesh).
"""

import math
from pathlib import Path
import shutil
from jinja2 import Environment, FileSystemLoader
from src.physics.atmosphere import FlightCondition


class OpenFOAMCaseGenerator:
    """Genera la estructura completa de un caso de OpenFOAM a partir de plantillas."""

    def __init__(self, templates_dir: Path | str | None = None) -> None:
        if templates_dir is None:
            self.templates_dir = Path(__file__).resolve().parent.parent.parent / "templates" / "openfoam"
        else:
            self.templates_dir = Path(templates_dir)

        if not self.templates_dir.exists():
            raise FileNotFoundError(f"Directorio de plantillas no encontrado: {self.templates_dir}")

        self.env = Environment(
            loader=FileSystemLoader(str(self.templates_dir)),
            trim_blocks=True,
            lstrip_blocks=True
        )

    def generate_case(
        self,
        output_case_dir: Path | str,
        flight_cond: FlightCondition,
        rocket_patches: list[str] | None = None,
        ref_length: float = 1.0,
        ref_area: float = 0.05,
        end_time: float = 0.05,
        initial_delta_t: float = 1e-7,
        write_interval: int = 500,
        is_axisymmetric: bool = False,
        inviscid: bool = False
    ) -> Path:
        """
        Crea o actualiza el directorio del caso con todos los diccionarios renderizados.
        """
        out_dir = Path(output_case_dir)
        dir_0 = out_dir / "0"
        dir_constant = out_dir / "constant"
        dir_system = out_dir / "system"

        for d in [dir_0, dir_constant, dir_system]:
            d.mkdir(parents=True, exist_ok=True)

        if rocket_patches is None:
            rocket_patches = ["nose", "fuselage", "fins", "base"]

        ux, uy, uz = flight_cond.velocity_vector_m_s

        # Estimación de turbulencia en freestream (Intensidad 1%)
        turb_intensity = 0.01
        u_mag = flight_cond.velocity_magnitude_m_s
        k_inf = 1.5 * (u_mag * turb_intensity) ** 2
        c_mu = 0.09
        turb_length_scale = 0.07 * ref_length
        omega_inf = (k_inf ** 0.5) / ((c_mu ** 0.25) * turb_length_scale)

        # Vector de sustentación y resistencia según alfa
        alpha_rad = math.radians(flight_cond.alpha_deg)
        drag_dir_x = math.cos(alpha_rad)
        drag_dir_y = math.sin(alpha_rad)
        lift_dir_x = -math.sin(alpha_rad)
        lift_dir_y = math.cos(alpha_rad)

        context = {
            "mach": flight_cond.mach,
            "altitude_m": flight_cond.altitude_m,
            "alpha_deg": flight_cond.alpha_deg,
            "velocity_mag": u_mag,
            "ux": ux,
            "uy": uy,
            "uz": uz,
            "p_inf": flight_cond.pressure_pa,
            "t_inf": flight_cond.temperature_k,
            "density_kg_m3": flight_cond.density_kg_m3,
            "k_inf": max(k_inf, 1e-4),
            "omega_inf": max(omega_inf, 1.0),
            "rocket_patches": rocket_patches,
            "ref_length": ref_length,
            "ref_area": ref_area,
            "end_time": end_time,
            "initial_delta_t": initial_delta_t,
            "write_interval": write_interval,
            "cofr_x": 0.5 * ref_length,
            "cofr_y": 0.0,
            "cofr_z": 0.0,
            "drag_dir_x": drag_dir_x,
            "drag_dir_y": drag_dir_y,
            "lift_dir_x": lift_dir_x,
            "lift_dir_y": lift_dir_y,
            "is_axisymmetric": is_axisymmetric,
            "inviscid": inviscid
        }

        # 1. constant/
        self._render_template("constant/thermophysicalProperties.jinja", dir_constant / "thermophysicalProperties", context)
        self._render_template("constant/turbulenceProperties.jinja", dir_constant / "turbulenceProperties", context)

        # 2. system/
        self._render_template("system/fvSchemes.jinja", dir_system / "fvSchemes", context)
        self._render_template("system/fvSolution.jinja", dir_system / "fvSolution", context)
        self._render_template("system/controlDict.jinja", dir_system / "controlDict", context)

        # 3. 0/
        self._render_template("0.orig/U.jinja", dir_0 / "U", context)
        self._render_template("0.orig/p.jinja", dir_0 / "p", context)
        self._render_template("0.orig/T.jinja", dir_0 / "T", context)
        self._render_template("0.orig/k.jinja", dir_0 / "k", context)
        self._render_template("0.orig/omega.jinja", dir_0 / "omega", context)

        return out_dir

    def generate_cone_benchmark_case(
        self,
        output_case_dir: Path | str,
        flight_cond: FlightCondition,
        cone_angle_deg: float = 15.0,
        cone_length_m: float = 0.5,
        nx_cone: int = 120,
        ny: int = 80,
        end_time: float = 0.004
    ) -> Path:
        """
        Genera el caso canónico de validación de un cono supersónico axisimétrico (cuña 3D de 5 grados)
        para comparar resultados de OpenFOAM frente a la solución exacta de Taylor-Maccoll.
        """
        out_dir = self.generate_case(
            output_case_dir=output_case_dir,
            flight_cond=flight_cond,
            rocket_patches=["cone"],
            ref_length=cone_length_m,
            ref_area=math.pi * (cone_length_m * math.tan(math.radians(cone_angle_deg))) ** 2,
            end_time=end_time,
            initial_delta_t=5e-7,
            write_interval=200,
            is_axisymmetric=True,
            inviscid=True
        )

        wedge_half_rad = math.radians(2.5)
        theta_c_rad = math.radians(cone_angle_deg)

        x_inlet = -0.15 * cone_length_m
        x_outlet = cone_length_m
        r_cone = cone_length_m * math.tan(theta_c_rad)
        r_top = 1.1 * cone_length_m

        y_cone = r_cone * math.cos(wedge_half_rad)
        z_cone = r_cone * math.sin(wedge_half_rad)
        y_top = r_top * math.cos(wedge_half_rad)
        z_top = r_top * math.sin(wedge_half_rad)

        nx_up = max(int(nx_cone * abs(x_inlet) / cone_length_m), 20)

        mesh_context = {
            "x_inlet": round(x_inlet, 6),
            "x_outlet": round(x_outlet, 6),
            "y_cone": round(y_cone, 8),
            "z_cone": round(z_cone, 8),
            "y_top": round(y_top, 8),
            "z_top": round(z_top, 8),
            "nx_up": nx_up,
            "nx_cone": nx_cone,
            "ny": ny
        }

        self._render_template(
            "system/blockMeshDict_cone.jinja",
            out_dir / "system" / "blockMeshDict",
            mesh_context
        )

        return out_dir

    def generate_sphere_benchmark_case(
        self,
        output_case_dir: Path | str,
        flight_cond: FlightCondition,
        sphere_radius_m: float = 0.10,
        n_theta_front: int = 60,
        n_theta_wake: int = 60,
        n_radial: int = 70,
        end_time: float = 0.0035
    ) -> Path:
        """
        Genera el caso canónico de una esfera supersónica (Blunt Body) con malla O-Grid polar.
        Captura el choque en arco desprendido (Bow Shock) aguas arriba, la línea sónica
        y el desarrollo de la estela supersónica aguas abajo.
        """
        out_dir = self.generate_case(
            output_case_dir=output_case_dir,
            flight_cond=flight_cond,
            rocket_patches=["sphere"],
            ref_length=2.0 * sphere_radius_m,
            ref_area=math.pi * (sphere_radius_m ** 2),
            end_time=end_time,
            initial_delta_t=3e-7,
            write_interval=250,
            is_axisymmetric=True,
            inviscid=True
        )

        wedge_half_rad = math.radians(2.5)
        cos_w = math.cos(wedge_half_rad)
        sin_w = math.sin(wedge_half_rad)

        # Centro de la esfera en x = R para que el punto de remanso (proa) esté exactamente en x = 0.0
        r_sph = sphere_radius_m
        x_center = r_sph
        x_nose = 0.0
        x_tail = 2.0 * r_sph

        # Dominio exterior: 2.0*R aguas arriba de la proa y 4.5*R aguas abajo en la estela
        r_out_front = 3.0 * r_sph
        r_out_top = 3.8 * r_sph
        r_out_wake = 4.5 * r_sph

        x_inlet = x_center - r_out_front
        x_outlet = x_center + r_out_wake

        # Puntos intermedios a 45° y 135° para los arcos de blockMesh
        c45 = math.cos(math.radians(45.0))
        s45 = math.sin(math.radians(45.0))

        r_mid_front_out = 0.5 * (r_out_front + r_out_top)
        r_mid_rear_out = 0.5 * (r_out_top + r_out_wake)

        mesh_context = {
            "x_nose": round(x_nose, 6),
            "x_inlet": round(x_inlet, 6),
            "x_center": round(x_center, 6),
            "x_tail": round(x_tail, 6),
            "x_outlet": round(x_outlet, 6),
            "y_eq": round(r_sph * cos_w, 8),
            "z_eq": round(r_sph * sin_w, 8),
            "y_top": round(r_out_top * cos_w, 8),
            "z_top": round(r_out_top * sin_w, 8),
            "x_arc_front_wall": round(x_center - r_sph * c45, 8),
            "y_arc_front_wall": round((r_sph * s45) * cos_w, 8),
            "z_arc_front_wall": round((r_sph * s45) * sin_w, 8),
            "x_arc_front_out": round(x_center - r_mid_front_out * c45, 8),
            "y_arc_front_out": round((r_mid_front_out * s45) * cos_w, 8),
            "z_arc_front_out": round((r_mid_front_out * s45) * sin_w, 8),
            "x_arc_rear_wall": round(x_center + r_sph * c45, 8),
            "y_arc_rear_wall": round((r_sph * s45) * cos_w, 8),
            "z_arc_rear_wall": round((r_sph * s45) * sin_w, 8),
            "x_arc_rear_out": round(x_center + r_mid_rear_out * c45, 8),
            "y_arc_rear_out": round((r_mid_rear_out * s45) * cos_w, 8),
            "z_arc_rear_out": round((r_mid_rear_out * s45) * sin_w, 8),
            "n_theta_front": n_theta_front,
            "n_theta_wake": n_theta_wake,
            "n_radial": n_radial,
            "radial_grading": 2.5
        }

        self._render_template(
            "system/blockMeshDict_sphere.jinja",
            out_dir / "system" / "blockMeshDict",
            mesh_context
        )

        return out_dir

    def generate_rocket_stl_case(
        self,
        output_case_dir: Path | str,
        flight_cond: FlightCondition,
        stl_files: list[Path] | None = None,
        rocket_length_m: float = 1.0,
        rocket_radius_m: float = 0.08,
        end_time: float = 0.02
    ) -> Path:
        """
        Prepara el caso 3D completo del cohete lanzador con un amplio dominio aguas arriba (3*L)
        para el Bow Shock de la ojiva redondeada y aguas abajo (8*L) para el desarrollo de la estela.
        Copia los archivos STL a constant/triSurface/ y genera blockMeshDict y snappyHexMeshDict.
        """
        patches = [p.stem for p in stl_files] if stl_files else ["nose", "fuselage", "fins", "base"]

        out_dir = self.generate_case(
            output_case_dir=output_case_dir,
            flight_cond=flight_cond,
            rocket_patches=patches,
            ref_length=rocket_length_m,
            ref_area=math.pi * (rocket_radius_m ** 2),
            end_time=end_time,
            initial_delta_t=1e-7,
            write_interval=250,
            is_axisymmetric=False,
            inviscid=False
        )

        tri_dir = out_dir / "constant" / "triSurface"
        tri_dir.mkdir(parents=True, exist_ok=True)
        if stl_files:
            for stl_path in stl_files:
                if stl_path.exists():
                    shutil.copy2(stl_path, tri_dir / stl_path.name)

        # Túnel numérico con 3*L aguas arriba y 8*L aguas abajo para capturar la estela completa
        x_min = -3.0 * rocket_length_m
        x_max = 8.0 * rocket_length_m
        y_max = 4.0 * rocket_length_m
        z_min = -4.0 * rocket_length_m
        z_max = 4.0 * rocket_length_m

        block_ctx = {
            "x_min": round(x_min, 4),
            "x_max": round(x_max, 4),
            "y_max": round(y_max, 4),
            "z_min": round(z_min, 4),
            "z_max": round(z_max, 4),
            "nx": 88,
            "ny": 32,
            "nz": 64,
        }
        self._render_template(
            "system/blockMeshDict_rocket.jinja",
            out_dir / "system" / "blockMeshDict",
            block_ctx
        )

        snappy_ctx = {
            "rocket_patches": patches,
            "box_x_min": round(-0.5 * rocket_length_m, 4),
            "box_x_max": round(4.5 * rocket_length_m, 4),
            "box_y_max": round(0.8 * rocket_length_m, 4),
            "box_z_min": round(-0.8 * rocket_length_m, 4),
            "box_z_max": round(0.8 * rocket_length_m, 4),
            "loc_x": round(-1.5 * rocket_length_m, 4),
            "loc_y": round(1.5 * rocket_length_m, 4),
            "loc_z": 0.0,
        }
        self._render_template(
            "system/snappyHexMeshDict.jinja",
            out_dir / "system" / "snappyHexMeshDict",
            snappy_ctx
        )

        return out_dir

    def _render_template(self, template_rel_path: str, output_path: Path, context: dict) -> None:
        template = self.env.get_template(template_rel_path)
        content = template.render(**context)
        with open(output_path, "w", encoding="utf-8", newline="\n") as f:
            f.write(content)
