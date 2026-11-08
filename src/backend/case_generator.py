"""
Generador automatizado de casos de OpenFOAM para rhoCentralFoam.
Utiliza Jinja2 para inyectar condiciones de vuelo calculadas por AtmosphereISA
y parámetros geométricos en los diccionarios de OpenFOAM.
"""

import math
from pathlib import Path
from jinja2 import Environment, FileSystemLoader
from src.physics.atmosphere import FlightCondition


class OpenFOAMCaseGenerator:
    """Genera la estructura completa de un caso de OpenFOAM a partir de plantillas."""

    def __init__(self, templates_dir: Path | str | None = None) -> None:
        if templates_dir is None:
            # Directorio por defecto: templates/openfoam relativo a la raíz del proyecto
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

        # Geometría de la cuña axisimétrica (+/- 2.5 grados alrededor del plano X-Y)
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

    def _render_template(self, template_rel_path: str, output_path: Path, context: dict) -> None:
        template = self.env.get_template(template_rel_path)
        content = template.render(**context)
        with open(output_path, "w", encoding="utf-8", newline="\n") as f:
            f.write(content)
