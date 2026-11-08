"""
Módulo de Advección Cinemática Lagrangiana de Partículas de Fluido (Fase 8.5).
Parametriza las líneas de corriente por tiempo físico de residencia real:
    t(s) = ∫ ds / |U(x(s))|
y actualiza en tiempo real (60 FPS) una nube 3D de partículas trazadoras mediante
interpolación vectorizada en NumPy para actualizar el Vertex Buffer Object (VBO) de la GPU.
"""

from dataclasses import dataclass
import math
import numpy as np
import pyvista as pv


@dataclass
class TrajectoryTrack:
    """Trayectoria lagrangiana individual parametrizada por tiempo normalizado tau en [0, 1]."""
    tau: np.ndarray        # Tiempo adimensional normalizado [0.0 .. 1.0], shape (M,)
    points: np.ndarray     # Coordenadas espaciales 3D (x, y, z) [m], shape (M, 3)
    mach: np.ndarray       # Número de Mach local en cada punto [-], shape (M,)


class LagrangianParticleSystem:
    """
    Motor cinemático que gestiona cientos de partículas de fluido moviéndose a la
    velocidad física local U(x, y, z) a través de ondas de choque y expansiones.
    """

    def __init__(
        self,
        streamlines_polydata: pv.PolyData,
        n_particles_per_track: int = 8,
        revolve_3d_azimuths: int = 6
    ) -> None:
        self.tracks: list[TrajectoryTrack] = []
        self.n_particles_per_track = max(n_particles_per_track, 1)
        self.phase: float = 0.0

        self._build_tracks_from_streamlines(streamlines_polydata, revolve_3d_azimuths)
        self.cloud_polydata = self._create_initial_polydata()

    def _build_tracks_from_streamlines(
        self,
        streamlines: pv.PolyData,
        revolve_3d_azimuths: int
    ) -> None:
        if streamlines is None or streamlines.n_cells == 0:
            return

        mach_all = (
            np.asarray(streamlines.point_data["Mach"], dtype=float)
            if "Mach" in streamlines.point_data
            else np.ones(streamlines.n_points, dtype=float)
        )
        u_all = (
            np.asarray(streamlines.point_data["U"], dtype=float)
            if "U" in streamlines.point_data
            else None
        )

        # Extraer cada polilínea individual del objeto PolyData de VTK
        base_tracks: list[TrajectoryTrack] = []
        for cell_idx in range(streamlines.n_cells):
            cell = streamlines.get_cell(cell_idx)
            pt_ids = cell.point_ids
            if len(pt_ids) < 4:
                continue

            pts = np.asarray(streamlines.points[pt_ids], dtype=float)
            mach_seg = mach_all[pt_ids]

            # Asegurar que la trayectoria está ordenada de aguas arriba (menor x) a aguas abajo (mayor x)
            if pts[0, 0] > pts[-1, 0]:
                pts = pts[::-1].copy()
                mach_seg = mach_seg[::-1].copy()

            # Calcular longitudes de arco ds entre puntos consecutivos
            diffs = np.diff(pts, axis=0)
            ds = np.linalg.norm(diffs, axis=1)

            # Filtrar puntos duplicados
            valid_mask = np.concatenate(([True], ds > 1e-7))
            pts = pts[valid_mask]
            mach_seg = mach_seg[valid_mask]
            if len(pts) < 4:
                continue

            diffs = np.diff(pts, axis=0)
            ds = np.linalg.norm(diffs, axis=1)

            # Velocidad local |U| para integrar el tiempo físico real dt = ds / |U|
            if u_all is not None:
                u_seg = np.linalg.norm(u_all[pt_ids][valid_mask], axis=1)
                u_mid = np.maximum(0.5 * (u_seg[:-1] + u_seg[1:]), 1.0)
            else:
                u_mid = np.maximum(0.5 * (mach_seg[:-1] + mach_seg[1:]), 0.1)

            dt = ds / u_mid
            t_cum = np.concatenate(([0.0], np.cumsum(dt)))
            t_total = float(t_cum[-1])
            if t_total <= 1e-12:
                continue

            tau = t_cum / t_total
            base_tracks.append(TrajectoryTrack(tau=tau, points=pts, mach=mach_seg))

        if not base_tracks:
            return

        # Distribuir las trayectorias en el volumen 3D (planos azimutales alrededor del eje X)
        azimuth_angles_deg = np.linspace(0.0, 360.0, revolve_3d_azimuths, endpoint=False)
        for track in base_tracks:
            for phi_deg in azimuth_angles_deg:
                phi_rad = math.radians(float(phi_deg))
                cos_p = math.cos(phi_rad)
                sin_p = math.sin(phi_rad)

                pts_rot = track.points.copy()
                y_orig = pts_rot[:, 1].copy()
                z_orig = pts_rot[:, 2].copy()
                pts_rot[:, 1] = y_orig * cos_p - z_orig * sin_p
                pts_rot[:, 2] = y_orig * sin_p + z_orig * cos_p

                self.tracks.append(
                    TrajectoryTrack(
                        tau=track.tau,
                        points=pts_rot,
                        mach=track.mach.copy()
                    )
                )

    def _create_initial_polydata(self) -> pv.PolyData:
        pts, mach = self.evaluate_at_phase(0.0)
        if len(pts) == 0:
            pts = np.zeros((1, 3), dtype=float)
            mach = np.zeros(1, dtype=float)
        cloud = pv.PolyData(pts)
        cloud.point_data["Mach"] = mach
        return cloud

    def evaluate_at_phase(self, phase: float) -> tuple[np.ndarray, np.ndarray]:
        """
        Evalúa la posición 3D y el Mach local de todas las partículas en el instante de fase [0, 1].
        """
        if not self.tracks:
            return np.empty((0, 3), dtype=float), np.empty(0, dtype=float)

        n_tracks = len(self.tracks)
        k = self.n_particles_per_track
        total_particles = n_tracks * k

        all_pts = np.empty((total_particles, 3), dtype=float)
        all_mach = np.empty(total_particles, dtype=float)

        base_offsets = np.linspace(0.0, 1.0, k, endpoint=False)

        idx = 0
        for t_idx, track in enumerate(self.tracks):
            # Desfase determinista suave entre líneas adyacentes para evitar frentes artificiales rígidos
            stagger = (t_idx * 0.173) % (1.0 / k)
            taus = (base_offsets + phase + stagger) % 1.0

            all_pts[idx:idx + k, 0] = np.interp(taus, track.tau, track.points[:, 0])
            all_pts[idx:idx + k, 1] = np.interp(taus, track.tau, track.points[:, 1])
            all_pts[idx:idx + k, 2] = np.interp(taus, track.tau, track.points[:, 2])
            all_mach[idx:idx + k] = np.interp(taus, track.tau, track.mach)
            idx += k

        return all_pts, all_mach

    def step(self, delta_phase: float = 0.012) -> pv.PolyData:
        """
        Avanza la fase temporal de las partículas y actualiza in-place el PolyData en la GPU.
        """
        self.phase = (self.phase + delta_phase) % 1.0
        pts, mach = self.evaluate_at_phase(self.phase)
        if len(pts) > 0 and self.cloud_polydata.n_points == len(pts):
            self.cloud_polydata.points = pts
            self.cloud_polydata.point_data["Mach"] = mach
        return self.cloud_polydata
