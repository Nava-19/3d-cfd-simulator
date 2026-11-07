"""
Tests unitarios para el parser de logs de rhoCentralFoam.
"""

from src.backend.log_parser import OpenFOAMLogParser


SAMPLE_LOG_CHUNK = """
Time = 1e-05

Courant Number mean: 0.123 max: 0.654
diagonal:  Solving for rho, Initial residual = 0, Final residual = 0, No Iterations 0
smoothSolver:  Solving for Ux, Initial residual = 0.054321, Final residual = 0.00012, No Iterations 4
smoothSolver:  Solving for Uy, Initial residual = 0.043210, Final residual = 0.00015, No Iterations 4
smoothSolver:  Solving for e, Initial residual = 0.087654, Final residual = 0.00021, No Iterations 5
smoothSolver:  Solving for k, Initial residual = 0.012345, Final residual = 0.00001, No Iterations 2
smoothSolver:  Solving for omega, Initial residual = 0.023456, Final residual = 0.00002, No Iterations 2
ExecutionTime = 0.45 s  ClockTime = 1 s

Time = 2e-05
"""

SAMPLE_FORCE_LINE = "1e-05    0.4521    0.0000    0.0823    0.0000    -0.0215    0.0000    0.0521    0.4000    0.0010    0.0813"


def test_log_parser_extracts_residuals():
    parser = OpenFOAMLogParser()
    completed_iter = None

    for line in SAMPLE_LOG_CHUNK.strip().split("\n"):
        res = parser.parse_line(line)
        if res is not None:
            completed_iter = res

    assert completed_iter is not None
    assert completed_iter.time == 1e-05
    assert completed_iter.courant_max == 0.654
    assert completed_iter.residuals["Ux"] == 0.054321
    assert completed_iter.residuals["Uy"] == 0.043210
    assert completed_iter.residuals["e"] == 0.087654
    assert completed_iter.residuals["k"] == 0.012345


def test_log_parser_extracts_forces():
    parser = OpenFOAMLogParser()
    forces = parser.parse_force_coeffs_line(SAMPLE_FORCE_LINE)

    assert forces is not None
    assert forces["time"] == 1e-05
    assert forces["cd"] == 0.4521
    assert forces["cl"] == 0.0823
    assert forces["cm_pitch"] == -0.0215
    assert forces["cd_pressure"] == 0.4000
    assert forces["cd_friction"] == 0.0521
