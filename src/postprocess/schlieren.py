"""
Módulo de Efecto Schlieren Numérico (Fase 7.1).
Calcula el módulo del gradiente espacial de densidad |∇ρ| en mallas 3D de OpenFOAM
para visualizar ondas de choque oblicuas, ondas de expansión de Prandtl-Meyer y estelas.

Fundamento físico (Relación de Gladstone-Dale):
    n - 1 = K_GD * ρ
La deflexión angular de un haz de luz en un túnel de viento óptico es proporcional a ∇n,
y por tanto directamente proporcional al gradiente de densidad:
    |∇ρ| = sqrt((∂ρ/∂x)² + (∂ρ/∂y)² + (∂ρ/∂z)²)
"""

import numpy as np
import pyvista as pv


R_AIR: float = 287.05287  # Constante específica del aire [J/(kg·K)]
GAMMA: float = 1.4        # Coeficiente adiabático del aire


def ensure_thermodynamic_fields(dataset: pv.DataSet) -> pv.DataSet:
    """
    Garantiza que el dataset de PyVista disponga de los campos 'rho' (densidad [kg/m³])
    y 'Mach' (número de Mach local [-]) tanto en cell_data como en point_data.
    Si 'rho' no fue exportado directamente por el solver, lo reconstruye mediante
    la ecuación de estado de gas ideal: ρ = p / (R_air * T).
    """
    # Convertir cell_data a point_data si aún no están interpolados a los nodos
    # (necesario para filtros continuos de derivadas e isocontornos en VTK)
    if "p" in dataset.cell_data and "p" not in dataset.point_data:
        dataset = dataset.cell_data_to_point_data(pass_cell_data=True)

    # 1. Densidad en datos de celda (Cell Data)
    if "rho" not in dataset.cell_data and "p" in dataset.cell_data and "T" in dataset.cell_data:
        p_cell = np.asarray(dataset.cell_data["p"], dtype=float)
        t_cell = np.maximum(np.asarray(dataset.cell_data["T"], dtype=float), 1.0)
        dataset.cell_data["rho"] = p_cell / (R_AIR * t_cell)

    # 2. Densidad en datos de nodo (Point Data)
    if "rho" not in dataset.point_data and "p" in dataset.point_data and "T" in dataset.point_data:
        p_pt = np.asarray(dataset.point_data["p"], dtype=float)
        t_pt = np.maximum(np.asarray(dataset.point_data["T"], dtype=float), 1.0)
        dataset.point_data["rho"] = p_pt / (R_AIR * t_pt)

    # 3. Número de Mach local M = |U| / sqrt(gamma * R * T)
    if "U" in dataset.point_data and "T" in dataset.point_data and "Mach" not in dataset.point_data:
        u_vec = np.asarray(dataset.point_data["U"], dtype=float)
        u_mag = np.linalg.norm(u_vec, axis=1)
        t_pt = np.maximum(np.asarray(dataset.point_data["T"], dtype=float), 1.0)
        a_local = np.sqrt(GAMMA * R_AIR * t_pt)
        dataset.point_data["Mach"] = u_mag / a_local

    if "U" in dataset.cell_data and "T" in dataset.cell_data and "Mach" not in dataset.cell_data:
        u_vec = np.asarray(dataset.cell_data["U"], dtype=float)
        u_mag = np.linalg.norm(u_vec, axis=1)
        t_cell = np.maximum(np.asarray(dataset.cell_data["T"], dtype=float), 1.0)
        a_local = np.sqrt(GAMMA * R_AIR * t_cell)
        dataset.cell_data["Mach"] = u_mag / a_local

    return dataset


def compute_numerical_schlieren(
    internal_mesh: pv.DataSet,
    slice_normal: tuple[float, float, float] | None = (0.0, 0.0, 1.0),
    slice_origin: tuple[float, float, float] | None = (0.0, 0.0, 0.0)
) -> pv.DataSet:
    """
    Calcula el campo de Schlieren numérico sobre la malla interna 3D y opcionalmente
    extrae el plano de simetría longitudinal (por defecto Z=0, plano X-Y).

    Añade los campos escalares:
        - 'grad_rho_mag': |∇ρ| = sqrt((∂ρ/∂x)² + (∂ρ/∂y)² + (∂ρ/∂z)²) [kg/m⁴]
        - 'schlieren_log': log10(1 + |∇ρ|) para realzar expansiones débiles junto a choques fuertes.
    """
    mesh_with_rho = ensure_thermodynamic_fields(internal_mesh.copy())

    # Calcular el vector gradiente ∇ρ en cada vértice mediante funciones de forma de elementos finitos VTK
    deriv_mesh = mesh_with_rho.compute_derivative(
        scalars="rho",
        gradient="grad_rho",
        preference="point"
    )

    grad_vec = np.asarray(deriv_mesh.point_data["grad_rho"], dtype=float)
    grad_mag = np.linalg.norm(grad_vec, axis=1)

    deriv_mesh.point_data["grad_rho_mag"] = grad_mag
    deriv_mesh.point_data["schlieren_log"] = np.log10(1.0 + grad_mag)

    if slice_normal is not None:
        # Si la malla es una cuña axisimétrica delgada (ej. benchmark cono +/- 2.5 grados),
        # el corte por z=0 es un plano meridiano exacto.
        sliced = deriv_mesh.slice(normal=slice_normal, origin=slice_origin)
        if sliced.n_points > 0:
            return sliced

    return deriv_mesh
