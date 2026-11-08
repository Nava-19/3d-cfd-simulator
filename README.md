# 3D Rocket Aerodynamics & Aeroheating CFD Simulator

[![Python](https://img.shields.io/badge/Python-3.13-blue.svg)](https://www.python.org/)
[![OpenFOAM](https://img.shields.io/badge/OpenFOAM-v2312%20(OpenCFD)-green.svg)](https://www.openfoam.com/)
[![GUI](https://img.shields.io/badge/GUI-PyQt6%20%7C%20PyVista%20%7C%20VTK-orange.svg)](https://pyvista.org/)
[![Tests](https://img.shields.io/badge/Tests-16%20Passed-brightgreen.svg)](tests/)
[![Architecture](https://img.shields.io/badge/Architecture-Asynchronous%20QThread-purple.svg)]()

Industrial-grade, automated computational fluid dynamics (CFD) pipeline designed for 3D rocket aerodynamics and aerodynamic heating across transonic and supersonic regimes ($\text{Mach } 0.2 \to 3.5+$).

Orchestrates **OpenFOAM (`rhoCentralFoam`)** in WSL2 through a decoupled Python backend, providing real-time convergence telemetry via `PyQtGraph` and GPU-accelerated 3D rendering with `PyVista`/`VTK`.

---

## 🚀 Key Engineering Capabilities

1. **International Standard Atmosphere (ISA 1976) Engine**
   - Automated conversion from `(Mach, Altitude, Angle of Attack)` to exact freestream thermodynamics ($T_\infty$, $p_\infty$, $\rho_\infty$, $\mu_\infty$, $a_\infty$, $\vec{U}_\infty$).
   - Sutherland's dynamic viscosity law and isentropic stagnation recovery temperatures ($T_{0}$).
   - Validated against official NASA/ICAO atmospheric tables.

2. **Analytical Ground Truth (Taylor-Maccoll Exact Solver)**
   - Custom non-linear ODE solver for supersonic axisymmetric cone flow using 4th-order Runge-Kutta (`scipy.integrate.solve_ivp`) and bracketed root-finding (`scipy.optimize.brentq`) for the weak shock solution.
   - Enables automated benchmark validation: extracts shock angles ($\theta_s$) and surface pressure ratios ($p_c / p_\infty$) to compute exact CFD deviation percentages.

3. **Headless OpenFOAM Orchestration (`rhoCentralFoam`)**
   - Dynamic Jinja2 dictionary generation for OpenCFD v2312: `thermophysicalProperties`, `turbulenceProperties` ($k$-$\omega$ SST), `fvSchemes` (Kurganov central-upwind shock capturing), `fvSolution`, and boundary conditions (`0/`).
   - Component-level aerodynamic forces breakdown (`nose`, `fuselage`, `fins`, `base`) via `forceCoeffs` for stability margin and center-of-pressure ($x_{cp}$) tracking.

4. **Asynchronous Execution & Live Telemetry**
   - Non-blocking `QThread` worker managing WSL2 subprocess execution with streaming stdout parsing.
   - Regex-powered real-time residual extraction ($p$, $U_x$, $U_y$, $e$, $k$, $\omega$, $\text{Courant}_{max}$) feeding live high-FPS plots.

5. **Advanced Post-Processing & 3D Visualization**
   - Numerical Schlieren imaging ($|\nabla \rho|$) for oblique shock wave and expansion fan identification.
   - 3D surface aeroheating maps ($T_w = T_{rec}$) on rocket skin.
   - 3D streamtracer particle lines colored by local Mach number.

---

## 🏗️ Architecture Overview

```
                      +---------------------------------------+
                      |         PyQt6 GUI / Dashboard         |
                      |   (Flight Inputs, Telemetry, 3D View) |
                      +-------------------+-------------------+
                                          |
                      +-------------------v-------------------+
                      |      Asynchronous QThread Worker      |
                      |       (src/backend/runner.py)         |
                      +---------+-------------------+---------+
                                |                   |
        +-----------------------v---+   +-----------v-----------------------+
        |   OpenFOAM Case Generator |   |     Real-Time Stream Parser       |
        |  (Jinja2 + ISA 1976 Cond) |   | (Regex: Residuals & Aero Coeffs)  |
        +-----------------------+---+   +-----------+-----------------------+
                                |                   |
                      +---------v-------------------v---------+
                      |     WSL2 Subprocess Execution         |
                      |   Ubuntu 22.04 LTS / OpenFOAM v2312   |
                      | (blockMesh -> snappyHexMesh -> solver)|
                      +-------------------+-------------------+
                                          |
                      +-------------------v-------------------+
                      |    PyVista / VTK 3D Engine (GPU)      |
                      | (Schlieren, Thermal Maps, Streamlines)|
                      +---------------------------------------+
```

---

## 📊 Physical Validation: OpenFOAM (`rhoCentralFoam`) vs. Exact Taylor-Maccoll Theory

Canonical axisymmetric supersonic cone benchmark ($\theta_c = 15.0^\circ$, $M_\infty = 2.00$, $11,200$-cell structured wedge mesh):

| Physical Quantity | Taylor-Maccoll (Analytical) | OpenFOAM (`rhoCentralFoam`) | Relative Error (%) |
| :--- | :---: | :---: | :---: |
| **Oblique Shock Wave Angle ($\beta$)** | `33.915°` | `34.387°` | **`1.39%`** |
| **Surface Pressure Ratio ($p_c / p_\infty$)** | `1.5699` | `1.5701` | **`0.01%`** |
| **Surface Mach Number ($M_c$)** | `1.7069` | `1.7146` | **`0.45%`** |

---

## 🧪 Verification & Unit Testing

The repository maintains automated test coverage across physics, templating, execution, and CFD benchmark validation:

```bash
# Run the complete test suite
pytest -v
```

```text
tests/test_atmosphere.py::test_sea_level_conditions PASSED               [  6%]
tests/test_atmosphere.py::test_troposphere_5000m PASSED                  [ 12%]
tests/test_atmosphere.py::test_tropopause_11000m PASSED                  [ 18%]
tests/test_atmosphere.py::test_lower_stratosphere_20000m PASSED          [ 25%]
tests/test_atmosphere.py::test_speed_of_sound_sea_level PASSED           [ 31%]
tests/test_atmosphere.py::test_flight_condition_supersonic_with_aoa PASSED [ 37%]
tests/test_atmosphere.py::test_altitude_out_of_bounds PASSED             [ 43%]
tests/test_case_generator.py::test_case_generator_renders_case PASSED    [ 50%]
tests/test_cone_benchmark.py::test_cone_benchmark_setup PASSED           [ 56%]
tests/test_cone_benchmark.py::test_cone_benchmark_cfd_vs_taylor_maccoll PASSED [ 62%]
tests/test_log_parser.py::test_log_parser_extracts_residuals PASSED      [ 68%]
tests/test_log_parser.py::test_log_parser_extracts_forces PASSED         [ 75%]
tests/test_runner.py::test_windows_to_wsl_path PASSED                    [ 81%]
tests/test_runner.py::test_openfoam_runner_instantiation PASSED          [ 87%]
tests/test_taylor_maccoll.py::test_taylor_maccoll_cone_15deg_mach_2 PASSED [ 93%]
tests/test_taylor_maccoll.py::test_taylor_maccoll_cone_20deg_mach_3 PASSED [100%]
============================= 16 passed in 2.41s ==============================
```

---

## 🛠️ Setup & Requirements

### Host Environment (Windows 11)
- Python 3.10+ (Tested on 3.13)
- Windows Subsystem for Linux (WSL2) with **Ubuntu 22.04 LTS**
- NVIDIA GPU (RTX 4060 Laptop or equivalent) for VTK rendering acceleration

### WSL2 Environment
- **OpenFOAM v2312** (OpenCFD / ESI distribution)
- OpenMPI

### Python Installation
```bash
# Clone the repository
git clone git@github.com:Nava-19/3d-cfd-simulator.git
cd 3d-cfd-simulator

# Create virtual environment & install dependencies
python -m venv venv
venv\Scripts\activate
pip install -r requirements.txt
```

---

## 📄 License
MIT License. Developed as a final-year aerospace engineering capstone project.
