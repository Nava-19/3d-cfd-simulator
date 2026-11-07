"""
Tests unitarios para el generador de casos OpenFOAM.
Verifica que los archivos se renderizan con los valores calculados por ISA.
"""

from pathlib import Path
import pytest
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

    # Verificar que los directorios existen
    assert (case_dir / "0" / "U").exists()
    assert (case_dir / "0" / "p").exists()
    assert (case_dir / "0" / "T").exists()
    assert (case_dir / "constant" / "thermophysicalProperties").exists()
    assert (case_dir / "system" / "controlDict").exists()

    # Verificar contenido renderizado en U
    u_content = (case_dir / "0" / "U").read_text(encoding="utf-8")
    ux, uy, uz = cond.velocity_vector_m_s
    assert f"({ux} {uy} {uz})" in u_content
    assert "nose" in u_content
    assert "fuselage" in u_content
    assert "fins" in u_content

    # Verificar contenido renderizado en controlDict
    ctrl_content = (case_dir / "system" / "controlDict").read_text(encoding="utf-8")
    assert "application         rhoCentralFoam;" in ctrl_content
    assert "forces_breakdown" in ctrl_content
