"""
Tests unitarios para el runner de OpenFOAM en WSL2.
"""

from pathlib import Path
from src.backend.runner import windows_to_wsl_path, OpenFOAMRunner


def test_windows_to_wsl_path():
    win_path = Path("C:/Users/34647/Desktop/CFD/cases/test_case")
    wsl_path = windows_to_wsl_path(win_path)
    assert wsl_path == "/mnt/c/Users/34647/Desktop/CFD/cases/test_case"


def test_openfoam_runner_instantiation(tmp_path: Path):
    runner = OpenFOAMRunner(case_dir_windows=tmp_path, num_cores=6)
    assert runner.wsl_distro == "Ubuntu-22.04"
    assert runner.num_cores == 6
    assert runner.case_dir_wsl == windows_to_wsl_path(tmp_path)
