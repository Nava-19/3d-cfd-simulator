"""
Módulo de Generación de Líneas de Corriente 3D (Streamlines) (Fase 7.3 y 8.5).
Integra numéricamente las trayectorias tangentes al campo vectorial de velocidad U(x, y, z)
mediante Runge-Kutta de orden 4/5 (VTK StreamTracer) y calcula el número de Mach local.

Ecuación diferencial de las líneas de corriente:
    dx / ds = U(x(s))   <=>   dx / Ux = dy / Uy = dz / Uz
"""

import numpy as np
import pyvista as pv

from src.postprocess.schlieren import ensure_thermodynamic_fields


def generate_streamlines(
    internal_mesh: pv.DataSet,
    n_seeds: int = 25,
    seed_point_a: tuple[float, float, float] | None = None,
    seed_point_b: tuple[float, float, float] | None = None,
    max_length: float | None = None,
    mirror_symmetry: bool = True
) -> pv.PolyData:
    """
    Traza líneas de corriente sobre el dominio fluido a partir de un peine (rake)
    de semillas situado aguas arriba del cuerpo.
    Si el dominio es axisimétrico (Y >= 0), refleja las líneas respecto a Y=0 para
    cubrir simétricamente tanto la parte superior como inferior del vehículo.
    """
    mesh = ensure_thermodynamic_fields(internal_mesh.copy())
    xmin, xmax, ymin, ymax, zmin, zmax = mesh.bounds

    dx = xmax - xmin
    dy = ymax - ymin
    dz = zmax - zmin

    is_thin_wedge = dz < 0.15 * max(dx, dy)

    if seed_point_a is None or seed_point_b is None:
        x_seed = xmin + 0.02 * dx
        y_low = ymin + 0.03 * dy
        y_high = ymin + 0.88 * dy
        z_mid = 0.0 if (zmin <= 0.0 <= zmax) else 0.5 * (zmin + zmax)
        seed_point_a = (x_seed, y_low, z_mid)
        seed_point_b = (x_seed, y_high, z_mid)

    if max_length is None:
        max_length = 2.5 * dx

    if is_thin_wedge:
        slice_plane = mesh.slice(normal=(0.0, 0.0, 1.0), origin=(0.0, 0.0, 0.0))
        if slice_plane.n_points > 0 and "U" in slice_plane.point_data:
            u_plane = np.asarray(slice_plane.point_data["U"], dtype=float).copy()
            u_plane[:, 2] = 0.0
            slice_plane.point_data["U_planar"] = u_plane
            streams = slice_plane.streamlines(
                vectors="U_planar",
                pointa=seed_point_a,
                pointb=seed_point_b,
                n_points=n_seeds,
                max_length=max_length,
                surface_streamlines=True
            )
            streams = ensure_thermodynamic_fields(streams)
            if mirror_symmetry and ymin >= -1e-4 and streams.n_points > 0:
                streams_bot = streams.copy(deep=True)
                streams_bot.points[:, 1] *= -1.0
                if "U" in streams_bot.point_data:
                    u_bot = np.asarray(streams_bot.point_data["U"], dtype=float).copy()
                    u_bot[:, 1] *= -1.0
                    streams_bot.point_data["U"] = u_bot
                return streams.merge(streams_bot)
            return streams

    streams = mesh.streamlines(
        vectors="U",
        pointa=seed_point_a,
        pointb=seed_point_b,
        n_points=n_seeds,
        max_length=max_length,
        integrator_type=45
    )

    return ensure_thermodynamic_fields(streams)
