"""
Módulo de Post-procesado Aerotermodinámico y Calentamiento por Fricción (Fase 7.2).
Extrae las superficies sólidas del cohete/cono, mapea el campo de temperatura de pared T_w
y el coeficiente de presión Cp, y localiza el punto de máximo calentamiento térmico (remanso).

Fundamento físico:
1. Temperatura de estancamiento (remanso isentrópico / conservación de entalpía total h_0):
       T_0 = T_∞ * (1 + (γ - 1)/2 * M_∞²)
2. Temperatura adiabática de pared o de recuperación (fricción viscosa en capa límite):
       T_rec = T_e * (1 + r * (γ - 1)/2 * M_e²)
   donde el factor de recuperación térmica r depende del número de Prandtl (Pr ≈ 0.71 en aire):
       - Capa límite laminar:    r = Pr^(1/2) ≈ 0.843
       - Capa límite turbulenta: r = Pr^(1/3) ≈ 0.892
"""

from dataclasses import dataclass
import numpy as np
import pyvista as pv


GAMMA: float = 1.4
PRANDTL_AIR: float = 0.71


@dataclass(frozen=True)
class AeroheatingMetrics:
    """Resumen cuantitativo del análisis térmico superficial."""
    max_wall_temperature_k: float
    min_wall_temperature_k: float
    mean_wall_temperature_k: float
    hotspot_coordinates_m: tuple[float, float, float]
    stagnation_point_coordinates_m: tuple[float, float, float]
    theoretical_stagnation_temp_k: float
    theoretical_recovery_temp_turb_k: float
    theoretical_recovery_temp_lam_k: float


def compute_theoretical_thermal_limits(
    t_inf_k: float,
    mach_inf: float,
    gamma: float = GAMMA,
    prandtl: float = PRANDTL_AIR
) -> tuple[float, float, float]:
    """
    Calcula los límites teóricos de temperatura [K] para flujo compresible:
    Retorna:
        (T_stagnation, T_recovery_turbulent, T_recovery_laminar)
    """
    kinetic_factor = 0.5 * (gamma - 1.0) * (mach_inf ** 2)
    r_turb = prandtl ** (1.0 / 3.0)
    r_lam = prandtl ** 0.5

    t_stag = t_inf_k * (1.0 + kinetic_factor)
    t_rec_turb = t_inf_k * (1.0 + r_turb * kinetic_factor)
    t_rec_lam = t_inf_k * (1.0 + r_lam * kinetic_factor)
    return t_stag, t_rec_turb, t_rec_lam


def extract_aeroheating_surface(
    multiblock_mesh: pv.MultiBlock,
    target_patches: list[str] | None = None,
    t_inf_k: float = 288.15,
    p_inf_pa: float = 101325.0,
    q_inf_pa: float = 1.0,
    mach_inf: float = 2.0
) -> tuple[pv.PolyData, AeroheatingMetrics]:
    """
    Extrae y combina los parches de pared del vehículo espacial desde el bloque 'boundary'
    de OpenFOAM, añade el coeficiente de presión Cp y calcula las métricas térmicas.
    """
    if target_patches is None:
        target_patches = ["nose", "fuselage", "fins", "base", "cone", "sphere"]

    if "boundary" not in multiblock_mesh.keys():
        raise ValueError("El objeto MultiBlock no contiene el bloque 'boundary' de OpenFOAM.")

    boundary_block = multiblock_mesh["boundary"]
    collected_meshes: list[pv.DataSet] = []

    for patch_name in target_patches:
        if patch_name in boundary_block.keys():
            patch_ds = boundary_block[patch_name]
            if patch_ds is not None and patch_ds.n_cells > 0:
                collected_meshes.append(patch_ds.extract_surface(algorithm="dataset_surface"))

    if not collected_meshes:
        raise ValueError(
            f"No se encontró ninguno de los parches solicitados {target_patches} "
            f"en boundary (disponibles: {list(boundary_block.keys())})."
        )

    if len(collected_meshes) == 1:
        surface = collected_meshes[0].copy()
    else:
        surface = collected_meshes[0].merge(collected_meshes[1:])

    # Interpolar variables de celda a punto para sombreado Gouraud/Phong suave en GPU
    if "T" in surface.cell_data and "T" not in surface.point_data:
        surface = surface.cell_data_to_point_data(pass_cell_data=True)

    # Calcular Coeficiente de Presión adimensional Cp = (p - p_inf) / q_inf
    q_safe = max(q_inf_pa, 1e-6)
    if "p" in surface.point_data:
        p_pts = np.asarray(surface.point_data["p"], dtype=float)
        surface.point_data["Cp"] = (p_pts - p_inf_pa) / q_safe
    if "p" in surface.cell_data:
        p_cells = np.asarray(surface.cell_data["p"], dtype=float)
        surface.cell_data["Cp"] = (p_cells - p_inf_pa) / q_safe

    # Extraer temperaturas de pared y localizar punto caliente y proa (remanso)
    centers = surface.cell_centers().points
    t_vals = np.asarray(surface.cell_data["T"], dtype=float)

    idx_max_t = int(np.argmax(t_vals))
    idx_nose_tip = int(np.argmin(centers[:, 0]))

    hotspot_pt = tuple(float(c) for c in centers[idx_max_t])
    stag_pt = tuple(float(c) for c in centers[idx_nose_tip])

    t_stag, t_rec_turb, t_rec_lam = compute_theoretical_thermal_limits(
        t_inf_k=t_inf_k,
        mach_inf=mach_inf
    )

    metrics = AeroheatingMetrics(
        max_wall_temperature_k=float(np.max(t_vals)),
        min_wall_temperature_k=float(np.min(t_vals)),
        mean_wall_temperature_k=float(np.mean(t_vals)),
        hotspot_coordinates_m=(hotspot_pt[0], hotspot_pt[1], hotspot_pt[2]),
        stagnation_point_coordinates_m=(stag_pt[0], stag_pt[1], stag_pt[2]),
        theoretical_stagnation_temp_k=t_stag,
        theoretical_recovery_temp_turb_k=t_rec_turb,
        theoretical_recovery_temp_lam_k=t_rec_lam
    )

    return surface, metrics
