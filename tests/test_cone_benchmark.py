"""
Tests de integración y validación para el benchmark del cono supersónico
comparando OpenFOAM (rhoCentralFoam) frente a la solución exacta de Taylor-Maccoll.
"""

from pathlib import Path
import pytest
from src.validation.cone_benchmark import SupersonicConeBenchmark


def test_cone_benchmark_setup(tmp_path: Path):
    """Verifica que el caso canónico axisimétrico se genera con blockMeshDict y BCs de cuña."""
    bench = SupersonicConeBenchmark()
    case_dir, cond, tm_res = bench.setup_case(
        case_dir=tmp_path / "cone_setup_test",
        mach_inf=2.0,
        cone_angle_deg=15.0
    )

    assert (case_dir / "system" / "blockMeshDict").exists()
    assert (case_dir / "0" / "U").exists()
    u_text = (case_dir / "0" / "U").read_text(encoding="utf-8")
    assert "wedge_front" in u_text
    assert "slip" in u_text
    assert 33.5 < tm_res.shock_angle_deg < 34.5


def test_cone_benchmark_cfd_vs_taylor_maccoll():
    """
    Valida los resultados reales de OpenFOAM (rhoCentralFoam) del cono de 15° a Mach 2.0
    frente a la solución analítica exacta de Taylor-Maccoll (error exigido < 2%).
    """
    case_dir = Path(__file__).resolve().parent.parent / "cases" / "benchmark_cone_m2"
    if not case_dir.exists():
        pytest.skip("El caso benchmark_cone_m2 aún no ha sido ejecutado.")

    bench = SupersonicConeBenchmark()
    report = bench.evaluate_cfd_results(
        case_dir=case_dir,
        mach_inf=2.0,
        cone_angle_deg=15.0,
        cone_length_m=0.5
    )

    # Error relativo en salto de presión en pared < 1.0%
    assert report.error_pressure_ratio_pct < 1.0

    # Error relativo en número de Mach superficial < 1.5%
    assert report.error_mach_surface_pct < 1.5

    # Error relativo en ángulo de onda de choque oblicua < 2.0%
    assert report.error_shock_angle_pct < 2.0
