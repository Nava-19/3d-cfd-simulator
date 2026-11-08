"""
Módulo de Validación Física: Cono Supersónico (OpenFOAM rhoCentralFoam vs. Taylor-Maccoll).
Automatiza la generación del caso canónico axisimétrico, la extracción de resultados CFD
mediante PyVista/VTK y el cálculo riguroso de errores relativos frente a la teoría exacta.
"""

from dataclasses import dataclass
import math
from pathlib import Path
import numpy as np
import pyvista as pv

from src.physics.atmosphere import AtmosphereISA, FlightCondition
from src.physics.taylor_maccoll import TaylorMaccollSolver, TaylorMaccollResult
from src.backend.case_generator import OpenFOAMCaseGenerator


@dataclass(frozen=True)
class ValidationComparisonReport:
    """Informe cuantitativo de validación CFD vs. Teoría Analítica de Taylor-Maccoll."""
    mach_inf: float
    cone_angle_deg: float
    # Valores analíticos exactos (Taylor-Maccoll)
    theo_shock_angle_deg: float
    theo_pressure_ratio: float
    theo_mach_surface: float
    # Valores obtenidos por CFD (OpenFOAM rhoCentralFoam)
    cfd_shock_angle_deg: float
    cfd_pressure_ratio: float
    cfd_mach_surface: float
    # Errores relativos porcentuales [%]
    error_shock_angle_pct: float
    error_pressure_ratio_pct: float
    error_mach_surface_pct: float

    def to_markdown_table(self) -> str:
        """Devuelve una tabla en formato Markdown con el desglose de validación."""
        return (
            f"### Validación Cono Supersónico ($\\theta_c = {self.cone_angle_deg:.1f}^\\circ$, $M_\\infty = {self.mach_inf:.2f}$)\n\n"
            "| Magnitud Física | Taylor-Maccoll (Exacto) | OpenFOAM (`rhoCentralFoam`) | Error Relativo (%) |\n"
            "| :--- | :---: | :---: | :---: |\n"
            f"| **Ángulo de onda de choque ($\\beta$)** | `{self.theo_shock_angle_deg:.3f}°` | `{self.cfd_shock_angle_deg:.3f}°` | `{self.error_shock_angle_pct:.2f}%` |\n"
            f"| **Relación de presión ($p_c / p_\\infty$)** | `{self.theo_pressure_ratio:.4f}` | `{self.cfd_pressure_ratio:.4f}` | `{self.error_pressure_ratio_pct:.2f}%` |\n"
            f"| **Mach en superficie ($M_c$)** | `{self.theo_mach_surface:.4f}` | `{self.cfd_mach_surface:.4f}` | `{self.error_mach_surface_pct:.2f}%` |\n"
        )


class SupersonicConeBenchmark:
    """Orquesta el caso de validación del cono supersónico y extrae las métricas físicas."""

    GAMMA: float = 1.4
    R_AIR: float = 287.05287

    def __init__(self) -> None:
        self.isa = AtmosphereISA()
        self.tm_solver = TaylorMaccollSolver(gamma=self.GAMMA)
        self.case_gen = OpenFOAMCaseGenerator()

    def setup_case(
        self,
        case_dir: Path | str,
        mach_inf: float = 2.0,
        cone_angle_deg: float = 15.0,
        cone_length_m: float = 0.5,
        nx_cone: int = 120,
        ny: int = 80
    ) -> tuple[Path, FlightCondition, TaylorMaccollResult]:
        """
        Prepara el directorio del caso de OpenFOAM y calcula la solución analítica de referencia.
        """
        cond = self.isa.compute_flight_condition(mach=mach_inf, altitude_km=0.0, alpha_deg=0.0)
        tm_res = self.tm_solver.solve(mach_inf=mach_inf, cone_angle_deg=cone_angle_deg)

        # Tiempo físico suficiente para que el flujo supersónico recorra ~4 veces el dominio
        flow_transit_time = cone_length_m / cond.velocity_magnitude_m_s
        end_time = round(4.0 * flow_transit_time, 5)

        out_path = self.case_gen.generate_cone_benchmark_case(
            output_case_dir=case_dir,
            flight_cond=cond,
            cone_angle_deg=cone_angle_deg,
            cone_length_m=cone_length_m,
            nx_cone=nx_cone,
            ny=ny,
            end_time=end_time
        )

        return out_path, cond, tm_res

    def evaluate_cfd_results(
        self,
        case_dir: Path | str,
        mach_inf: float = 2.0,
        cone_angle_deg: float = 15.0,
        cone_length_m: float = 0.5
    ) -> ValidationComparisonReport:
        """
        Lee el último paso de tiempo del caso de OpenFOAM usando PyVista/VTK,
        extrae el ángulo de onda de choque, la presión en pared y el Mach superficial,
        y los compara frente a la solución exacta de Taylor-Maccoll.
        """
        case_path = Path(case_dir)
        foam_file = case_path / "case.foam"
        if not foam_file.exists():
            foam_file.touch()

        # 1. Solución analítica exacta (Ground Truth)
        cond = self.isa.compute_flight_condition(mach=mach_inf, altitude_km=0.0, alpha_deg=0.0)
        tm_res = self.tm_solver.solve(mach_inf=mach_inf, cone_angle_deg=cone_angle_deg)
        p_inf = cond.pressure_pa

        # 2. Lectura nativa del caso OpenFOAM con VTK/PyVista
        reader = pv.POpenFOAMReader(str(foam_file))
        if len(reader.time_values) < 2:
            raise RuntimeError(f"El caso en {case_path} no contiene resultados convergidos (time_values={reader.time_values}).")

        reader.set_active_time_value(reader.time_values[-1])
        mesh = reader.read()

        # 3. Extracción en la superficie del cono (parche 'cone')
        cone_patch = mesh["boundary"]["cone"]
        cone_centers = cone_patch.cell_centers().points
        x_cone = cone_centers[:, 0]

        # Filtrar región autosemejante (20% a 90% de la longitud del cono)
        mask_cone = (x_cone >= 0.20 * cone_length_m) & (x_cone <= 0.90 * cone_length_m)
        p_cone_vals = cone_patch.cell_data["p"][mask_cone]
        t_cone_vals = cone_patch.cell_data["T"][mask_cone]
        u_cone_vecs = cone_patch.cell_data["U"][mask_cone]

        cfd_pressure_ratio = float(np.mean(p_cone_vals) / p_inf)

        u_mag_cone = np.linalg.norm(u_cone_vecs, axis=1)
        a_cone = np.sqrt(self.GAMMA * self.R_AIR * t_cone_vals)
        cfd_mach_surface = float(np.mean(u_mag_cone / a_cone))

        # 4. Detección del ángulo de la onda de choque oblicua en la malla interna
        internal_mesh = mesh["internalMesh"]
        cell_pts = internal_mesh.cell_centers().points
        p_internal = internal_mesh.cell_data["p"]

        x_int = cell_pts[:, 0]
        r_int = np.sqrt(cell_pts[:, 1] ** 2 + cell_pts[:, 2] ** 2)

        # Umbral de mitad del salto de Rankine-Hugoniot en la onda de choque cónica
        # En flujo cónico 3D, el salto brusco en el choque (p2) es ~50% del aumento total hasta la pared (pc),
        # por lo que el punto medio de la discontinuidad del choque se sitúa en ~25% de (pc - p_inf).
        p_surface_mean = float(np.mean(p_cone_vals))
        p_shock_threshold = p_inf + 0.25 * (p_surface_mean - p_inf)

        # Muestrear estaciones axiales x_k entre 0.15*L y 0.85*L
        stations = np.linspace(0.18 * cone_length_m, 0.85 * cone_length_m, 25)
        dx_window = 0.015 * cone_length_m
        shock_angles_deg: list[float] = []

        for x_k in stations:
            col_mask = np.abs(x_int - x_k) <= dx_window
            if np.sum(col_mask) < 5:
                continue

            r_col = r_int[col_mask]
            p_col = p_internal[col_mask]
            x_col = x_int[col_mask]

            # Ordenar de menor a mayor radio r
            sort_idx = np.argsort(r_col)
            r_sorted = r_col[sort_idx]
            p_sorted = p_col[sort_idx]
            x_sorted = x_col[sort_idx]

            # Buscar el cruce donde la presión cae por debajo de p_shock_threshold al salir del cono de choque
            crossings = np.where((p_sorted[:-1] >= p_shock_threshold) & (p_sorted[1:] < p_shock_threshold))[0]
            if len(crossings) > 0:
                idx = crossings[-1]
                p1, p2 = p_sorted[idx], p_sorted[idx + 1]
                r1, r2 = r_sorted[idx], r_sorted[idx + 1]
                x_loc = 0.5 * (x_sorted[idx] + x_sorted[idx + 1])

                # Interpolación lineal sub-celda para localizar la posición exacta del choque
                weight = (p_shock_threshold - p1) / (p2 - p1) if abs(p2 - p1) > 1e-6 else 0.5
                r_shock = r1 + weight * (r2 - r1)
                angle_deg = math.degrees(math.atan2(r_shock, x_loc))
                shock_angles_deg.append(angle_deg)

        cfd_shock_angle_deg = float(np.median(shock_angles_deg))

        # 5. Cálculo de errores relativos porcentuales
        err_shock = abs(cfd_shock_angle_deg - tm_res.shock_angle_deg) / tm_res.shock_angle_deg * 100.0
        err_p = abs(cfd_pressure_ratio - tm_res.pressure_ratio_surface) / tm_res.pressure_ratio_surface * 100.0
        err_m = abs(cfd_mach_surface - tm_res.mach_surface) / tm_res.mach_surface * 100.0

        return ValidationComparisonReport(
            mach_inf=mach_inf,
            cone_angle_deg=cone_angle_deg,
            theo_shock_angle_deg=tm_res.shock_angle_deg,
            theo_pressure_ratio=tm_res.pressure_ratio_surface,
            theo_mach_surface=tm_res.mach_surface,
            cfd_shock_angle_deg=cfd_shock_angle_deg,
            cfd_pressure_ratio=cfd_pressure_ratio,
            cfd_mach_surface=cfd_mach_surface,
            error_shock_angle_pct=err_shock,
            error_pressure_ratio_pct=err_p,
            error_mach_surface_pct=err_m
        )
