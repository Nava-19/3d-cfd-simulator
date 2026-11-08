"""
Tests unitarios para el generador de casos OpenFOAM.
Verifica que los archivos se renderizan con los valores calculados por ISA
para el cono, la esfera supersónica y el cohete 3D con dominio de estela.
"""

from pathlib import Path
from src.physics.atmosphere import AtmosphereISA
from src.backend.case_generator import OpenFOAMCaseGenerator


def test_case_generator_renders_case(tmp_path: Path):
    isa = AtmosphereISA()
    cond = isa.compute_flight_condition(mach=2.0, altitude_km=10.0, alpha_deg=2.0)

    generator = OpenFOAMCaseGenerator()
    case_dir = generator.generate_case(
        output_case_dir=tmp_path / "test_rocket_case",
        flight_cond=cond,
        rocket_patches=["nose", "fuselage", "fins"]
    )

    assert (case_dir / "0" / "U").exists()
    assert (case_dir / "0" / "p").exists()
    assert (case_dir / "0" / "T").exists()
    assert (case_dir / "constant" / "thermophysicalProperties").exists()
    assert (case_dir / "system" / "controlDict").exists()

    u_content = (case_dir / "0" / "U").read_text(encoding="utf-8")
    ux, uy, uz = cond.velocity_vector_m_s
    assert f"({ux} {uy} {uz})" in u_content
    assert "nose" in u_content
    assert "fuselage" in u_content
    assert "fins" in u_content

    ctrl_content = (case_dir / "system" / "controlDict").read_text(encoding="utf-8")
    assert "application         rhoCentralFoam;" in ctrl_content
    assert "forces_breakdown" in ctrl_content


def test_case_generator_sphere_and_rocket_stl_with_wake(tmp_path: Path):
    isa = AtmosphereISA()
    cond = isa.compute_flight_condition(mach=2.0, altitude_km=0.0, alpha_deg=0.0)
    generator = OpenFOAMCaseGenerator()

    sphere_dir = generator.generate_sphere_benchmark_case(
        output_case_dir=tmp_path / "sphere_case",
        flight_cond=cond
    )
    assert (sphere_dir / "system" / "blockMeshDict").exists()
    assert "sphere" in (sphere_dir / "system" / "blockMeshDict").read_text(encoding="utf-8")

    rocket_dir = generator.generate_rocket_stl_case(
        output_case_dir=tmp_path / "rocket_wake_case",
        flight_cond=cond,
        rocket_length_m=1.0
    )
    assert (rocket_dir / "system" / "blockMeshDict").exists()
    assert (rocket_dir / "system" / "snappyHexMeshDict").exists()
    bm_text = (rocket_dir / "system" / "blockMeshDict").read_text(encoding="utf-8")
    # Verificar que el túnel deja -3*L aguas arriba y +8*L aguas abajo para la estela
    assert "-3.0" in bm_text and "8.0" in bm_text
