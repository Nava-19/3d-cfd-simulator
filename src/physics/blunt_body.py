"""
Módulo Analítico de Cuerpos Romos Supersónicos (Blunt Body & Detached Bow Shock).
Implementa:
    1. Relaciones de choque normal de Rankine-Hugoniot y fórmula de tubo de Pitot supersónico
       de Rayleigh para la presión de remanso p_02 y Cp_max tras el Bow Shock.
    2. Correlaciones hiperbólicas de Billig (1967) para la distancia de desprendimiento
       (standoff distance Δ) y la geometría 3D del choque en arco (Bow Shock) sobre una
       esfera u ojiva redondeada de radio R.
"""

from dataclasses import dataclass
import math
import numpy as np
import pyvista as pv


@dataclass(frozen=True)
class BluntBodyShockResult:
    """Solución teórica de choque desprendido (Bow Shock) sobre una esfera/ojiva de radio R."""
    mach_inf: float
    nose_radius_m: float
    standoff_distance_m: float           # Distancia de desprendimiento Δ delante de la proa [m]
    shock_vertex_x_m: float              # Coordenada axial del vértice del choque (-Δ) [m]
    shock_curvature_radius_m: float      # Radio de curvatura del choque en el eje Rc [m]
    mach_after_normal_shock: float       # Número de Mach subsónico M2 tras el choque normal [-]
    static_pressure_ratio_shock: float   # Salto de presión estática p2 / p_inf [-]
    stagnation_pressure_ratio: float     # Relación de Rayleigh Pitot p02 / p_inf en el punto de remanso [-]
    cp_stagnation_max: float             # Coeficiente de presión máximo en remanso Cp,max [-]


class BluntBodySolver:
    """Calcula las propiedades exactas de remanso y la forma 3D del Bow Shock (Billig, 1967)."""

    def __init__(self, gamma: float = 1.4) -> None:
        self.gamma = gamma

    def solve(self, mach_inf: float, nose_radius_m: float = 0.10) -> BluntBodyShockResult:
        if mach_inf <= 1.0:
            raise ValueError(f"El flujo debe ser supersónico (M > 1.0), recibido M={mach_inf}")

        g = self.gamma
        m2_sq = mach_inf ** 2

        # 1. Salto estático tras choque normal (Rankine-Hugoniot)
        p2_pinf = (2.0 * g * m2_sq - (g - 1.0)) / (g + 1.0)
        mach2 = math.sqrt(((g - 1.0) * m2_sq + 2.0) / (2.0 * g * m2_sq - (g - 1.0)))

        # 2. Fórmula supersónica de Rayleigh-Pitot (presión total en el punto de remanso de la nariz)
        term1 = (((g + 1.0) ** 2) * m2_sq / (4.0 * g * m2_sq - 2.0 * (g - 1.0))) ** (g / (g - 1.0))
        term2 = (1.0 - g + 2.0 * g * m2_sq) / (g + 1.0)
        p02_pinf = term1 * term2

        # 3. Coeficiente de presión máximo en el punto de remanso (Cp_max)
        cp_max = (2.0 / (g * m2_sq)) * (p02_pinf - 1.0)

        # 4. Correlación de Billig (1967) para esfera en flujo supersónico
        delta = nose_radius_m * 0.143 * math.exp(3.24 / m2_sq)
        rc = nose_radius_m * 1.143 * math.exp(0.54 / ((mach_inf - 1.0) ** 1.2))

        return BluntBodyShockResult(
            mach_inf=mach_inf,
            nose_radius_m=nose_radius_m,
            standoff_distance_m=delta,
            shock_vertex_x_m=-delta,
            shock_curvature_radius_m=rc,
            mach_after_normal_shock=mach2,
            static_pressure_ratio_shock=p2_pinf,
            stagnation_pressure_ratio=p02_pinf,
            cp_stagnation_max=cp_max
        )

    def build_3d_bow_shock_surface(
        self,
        mach_inf: float,
        nose_radius_m: float = 0.10,
        x_max_m: float = 0.45,
        n_axial: int = 50,
        n_theta: int = 72
    ) -> pv.PolyData:
        """
        Construye la superficie hiperboloide 3D de revolución del choque desprendido (Bow Shock)
        según la ecuación analítica de Billig (1967), con el punto de remanso de la esfera en x = 0.
        """
        res = self.solve(mach_inf=mach_inf, nose_radius_m=nose_radius_m)
        mu_rad = math.asin(1.0 / mach_inf)
        tan_mu = math.tan(mu_rad)
        cot_mu_sq = 1.0 / (tan_mu ** 2)

        delta = res.standoff_distance_m
        rc = res.shock_curvature_radius_m

        # Parametrizar x desde el vértice del choque (-delta) hasta x_max_m
        x_vals = np.linspace(-delta, x_max_m, n_axial)
        # Invertir la ecuación hiperbólica x(r) para obtener r(x):
        # x + delta = rc * cot^2(mu) * (sqrt(1 + r^2 * tan^2(mu) / rc^2) - 1)
        term = 1.0 + (x_vals + delta) / (rc * cot_mu_sq)
        r_vals = (rc / tan_mu) * np.sqrt(np.maximum(term ** 2 - 1.0, 0.0))

        # Generar malla estructurada de revolución 3D (360 grados alrededor del eje X)
        phis = np.linspace(0.0, 2.0 * math.pi, n_theta)
        xx, pp = np.meshgrid(x_vals, phis, indexing="ij")
        rr = np.repeat(r_vals[:, None], n_theta, axis=1)

        yy = rr * np.cos(pp)
        zz = rr * np.sin(pp)

        grid = pv.StructuredGrid(xx, yy, zz)
        return grid.extract_surface(algorithm="dataset_surface")


def build_3d_conical_shock_surface(
    shock_angle_deg: float,
    x_end_m: float = 0.50,
    n_axial: int = 40,
    n_theta: int = 72
) -> pv.PolyData:
    """
    Construye la superficie cónica 3D completa (360°) de la onda de choque oblicua
    de Taylor-Maccoll con vértice en (0, 0, 0) y semiángulo β_TM.
    """
    beta_rad = math.radians(shock_angle_deg)
    tan_b = math.tan(beta_rad)

    x_vals = np.linspace(0.0, x_end_m, n_axial)
    r_vals = x_vals * tan_b

    phis = np.linspace(0.0, 2.0 * math.pi, n_theta)
    xx, pp = np.meshgrid(x_vals, phis, indexing="ij")
    rr = np.repeat(r_vals[:, None], n_theta, axis=1)

    yy = rr * np.cos(pp)
    zz = rr * np.sin(pp)

    grid = pv.StructuredGrid(xx, yy, zz)
    return grid.extract_surface(algorithm="dataset_surface")
