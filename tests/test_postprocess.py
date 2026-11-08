"""
Tests unitarios y de integración para los algoritmos de Post-procesado 3D (Fase 7 y 8.5):
    - Efecto Schlieren numérico simétrico (|∇ρ|)
    - Extracción aerotérmica superficial (Aeroheating y Cp)
    - Trazado de líneas de corriente (Streamlines RK45)
    - Sistema Lagrangiano de Partículas de Fluido (LagrangianParticleSystem)
    - Exportador de Base de Datos Aerodinámica (AeroDB CSV/JSON)
"""

from pathlib import Path
import math
import numpy as np
import pyvista as pv

from src.postprocess.aeroheating import (
    compute_theoretical_thermal_limits,
    extract_aeroheating_surface,
)
from src.postprocess.aerodb_exporter import AeroDBExporter, AeroDBRecord
from src.postprocess.particles import LagrangianParticleSystem
from src.postprocess.schlieren import compute_numerical_schlieren
from src.postprocess.streamlines import generate_streamlines


BENCHMARK_CASE_DIR = Path(__file__).resolve().parent.parent / "cases" / "benchmark_cone_m2"


def _load_benchmark_multiblock() -> pv.MultiBlock:
    foam_file = BENCHMARK_CASE_DIR / "case.foam"
    if not foam_file.exists():
        foam_file.touch()
    reader = pv.POpenFOAMReader(str(foam_file))
    reader.set_active_time_value(reader.time_values[-1])
    return reader.read()


def test_theoretical_thermal_limits():
    t_stag, t_rec_turb, t_rec_lam = compute_theoretical_thermal_limits(
        t_inf_k=288.15, mach_inf=2.0
    )
    assert math.isclose(t_stag, 518.67, rel_tol=1e-4)
    assert t_rec_lam < t_rec_turb < t_stag


def test_numerical_schlieren_captures_symmetric_shock_gradient():
    mb = _load_benchmark_multiblock()
    internal = mb["internalMesh"]

    schlieren_slice = compute_numerical_schlieren(
        internal_mesh=internal,
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        mirror_symmetry=True
    )

    assert schlieren_slice.n_points > 200
    assert "grad_rho_mag" in schlieren_slice.point_data
    assert "schlieren_log" in schlieren_slice.point_data
    assert float(schlieren_slice.point_data["grad_rho_mag"].max()) > 10.0
    # Verificar que el plano reflejado cubre tanto Y positivo como Y negativo
    _, _, ymin, ymax, _, _ = schlieren_slice.bounds
    assert ymin < -0.2 and ymax > 0.2


def test_aeroheating_surface_extraction():
    mb = _load_benchmark_multiblock()
    surface, metrics = extract_aeroheating_surface(
        multiblock_mesh=mb,
        target_patches=["cone"],
        t_inf_k=288.15,
        p_inf_pa=101325.0,
        q_inf_pa=283710.0,
        mach_inf=2.0
    )

    assert surface.n_cells > 50
    assert "Cp" in surface.point_data
    assert "T" in surface.point_data
    assert metrics.max_wall_temperature_k > 320.0
    assert float(surface.point_data["Cp"].mean()) > 0.15


def test_streamlines_and_lagrangian_particles():
    mb = _load_benchmark_multiblock()
    internal = mb["internalMesh"]

    streams = generate_streamlines(internal_mesh=internal, n_seeds=12, mirror_symmetry=True)
    assert streams.n_points > 40
    assert "Mach" in streams.point_data
    assert float(streams.point_data["Mach"].max()) >= 1.9

    # Verificar el sistema cinemático de partículas 3D
    psys = LagrangianParticleSystem(
        streamlines_polydata=streams,
        n_particles_per_track=6,
        revolve_3d_azimuths=4
    )
    assert psys.cloud_polydata.n_points > 50
    pts_t0 = psys.cloud_polydata.points.copy()

    # Avanzar un paso de tiempo y comprobar que las partículas se desplazan en el espacio 3D
    psys.step(delta_phase=0.05)
    pts_t1 = psys.cloud_polydata.points.copy()
    displacement = float(np.linalg.norm(pts_t1 - pts_t0))
    assert displacement > 1e-3


def test_aerodb_exporter_csv_and_json(tmp_path: Path):
    x_cp = AeroDBExporter.compute_center_of_pressure(
        cd=0.25, cl=0.50, cm_pitch=-0.05, alpha_deg=5.0, x_ref_m=0.5, l_ref_m=1.0
    )
    assert x_cp > 0.5

    rec = AeroDBRecord(
        case_name="test_case_m2",
        mach=2.0,
        altitude_km=5.0,
        alpha_deg=2.0,
        p_inf_pa=54000.0,
        t_inf_k=255.7,
        u_inf_m_s=641.0,
        cd=0.21,
        cl=0.12,
        cm_pitch=-0.03,
        cd_pressure=0.18,
        cd_friction=0.03,
        x_cp_m=x_cp,
        max_wall_temp_k=445.0,
        stagnation_temp_theo_k=460.2
    )

    csv_file = AeroDBExporter.append_csv(rec, tmp_path / "aerodb.csv")
    json_file = AeroDBExporter.export_json(rec, tmp_path / "result.json")

    assert csv_file.exists()
    assert json_file.exists()
    assert "test_case_m2" in csv_file.read_text(encoding="utf-8")
