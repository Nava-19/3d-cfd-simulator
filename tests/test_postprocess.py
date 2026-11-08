"""
Tests unitarios y de integración para los algoritmos de Post-procesado 3D (Fase 7):
    - Efecto Schlieren numérico (|∇ρ|)
    - Extracción aerotérmica superficial (Aeroheating y Cp)
    - Trazado de líneas de corriente (Streamlines RK45)
    - Exportador de Base de Datos Aerodinámica (AeroDB CSV/JSON)
"""

from pathlib import Path
import math
import pyvista as pv

from src.postprocess.aeroheating import (
    compute_theoretical_thermal_limits,
    extract_aeroheating_surface,
)
from src.postprocess.aerodb_exporter import AeroDBExporter, AeroDBRecord
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
    # A M=2.0, T0 = 288.15 * (1 + 0.2 * 4) = 288.15 * 1.8 = 518.67 K
    assert math.isclose(t_stag, 518.67, rel_tol=1e-4)
    # T_rec_turb > T_rec_lam y ambos menores que T_stag (pues Pr = 0.71 < 1)
    assert t_rec_lam < t_rec_turb < t_stag


def test_numerical_schlieren_captures_shock_gradient():
    mb = _load_benchmark_multiblock()
    internal = mb["internalMesh"]

    schlieren_slice = compute_numerical_schlieren(
        internal_mesh=internal,
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0)
    )

    assert schlieren_slice.n_points > 100
    assert "grad_rho_mag" in schlieren_slice.point_data
    assert "schlieren_log" in schlieren_slice.point_data
    # En una onda de choque supersónica a M=2.0, el gradiente de densidad |∇ρ| supera ampliamente 10 kg/m^4
    assert float(schlieren_slice.point_data["grad_rho_mag"].max()) > 10.0


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
    # Tras la onda de choque cónica de 15° a M=2.0, la temperatura estática superficial sube por encima de 320 K
    assert metrics.max_wall_temperature_k > 320.0
    assert float(surface.point_data["Cp"].mean()) > 0.15


def test_streamlines_generation():
    mb = _load_benchmark_multiblock()
    internal = mb["internalMesh"]

    streams = generate_streamlines(internal_mesh=internal, n_seeds=15)
    assert streams.n_points > 20
    assert "Mach" in streams.point_data
    assert float(streams.point_data["Mach"].max()) >= 1.9


def test_aerodb_exporter_csv_and_json(tmp_path: Path):
    x_cp = AeroDBExporter.compute_center_of_pressure(
        cd=0.25, cl=0.50, cm_pitch=-0.05, alpha_deg=5.0, x_ref_m=0.5, l_ref_m=1.0
    )
    assert x_cp > 0.5  # Momento picado negativo -> Xcp detrás de Xref (estabilidad estática positiva)

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
