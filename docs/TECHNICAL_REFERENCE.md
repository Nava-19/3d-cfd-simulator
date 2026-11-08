# Technical Reference — 3D Rocket Aerodynamics & Aeroheating CFD Simulator

> **Version 1.0 — November 2026**
>
> This document collects, in one place, every governing partial differential
> equation, analytical ground-truth formulation, thermodynamic model, coordinate
> frame, numerical discretization scheme, meshing topology, and deliberate
> simplification used by the simulator. Its purpose is to make the project
> auditable as a rigorous aerospace engineering artefact — not just "code that
> runs" — and to provide a single source of truth for anyone reviewing,
> extending, or validating the simulation.
>
> This system generates comprehensive Mach-, altitude-, and angle-of-attack-dependent
> aerodynamic and thermal databases (`AeroDB`), replacing idealized drag approximations
> with high-fidelity compressible Navier-Stokes solutions across the launch envelope.

---

## Table of contents

1. [Coordinate systems and computational domains](#1-coordinate-systems-and-computational-domains)
2. [Units and physical constants](#2-units-and-physical-constants)
3. [Governing equations of compressible flow](#3-governing-equations-of-compressible-flow)
4. [Atmospheric model (ISA 1976) and freestream kinematics](#4-atmospheric-model-isa-1976-and-freestream-kinematics)
5. [Shock-capturing numerical formulation (`rhoCentralFoam`)](#5-shock-capturing-numerical-formulation-rhocentralfoam)
6. [Turbulence closure (k-ω SST) and boundary conditions](#6-turbulence-closure-k-ω-sst-and-boundary-conditions)
7. [Analytical ground truth I — Supersonic cone (Taylor-Maccoll ODE)](#7-analytical-ground-truth-i--supersonic-cone-taylor-maccoll-ode)
8. [Analytical ground truth II — Blunt body bow shock (Rankine-Hugoniot, Rayleigh-Pitot, Billig)](#8-analytical-ground-truth-ii--blunt-body-bow-shock-rankine-hugoniot-rayleigh-pitot-billig)
9. [Aerodynamic heating and recovery temperature model](#9-aerodynamic-heating-and-recovery-temperature-model)
10. [Meshing topologies and boundary layer resolution](#10-meshing-topologies-and-boundary-layer-resolution)
11. [Aerodynamic forces, center of pressure, and GNC export (`AeroDB`)](#11-aerodynamic-forces-center-of-pressure-and-gnc-export-aerodb)
12. [Post-processing: Numerical Schlieren and Lagrangian particle tracking](#12-post-processing-numerical-schlieren-and-lagrangian-particle-tracking)
13. [Asynchronous execution architecture and real-time telemetry](#13-asynchronous-execution-architecture-and-real-time-telemetry)
14. [Validation summary](#14-validation-summary)
15. [Known limitations and deliberate simplifications](#15-known-limitations-and-deliberate-simplifications)
16. [References](#16-references)

---

## 1. Coordinate systems and computational domains

The simulator operates in a right-handed **three-dimensional body-fixed Cartesian coordinate frame** aligned with the vehicle's longitudinal axis.

| Property | Convention |
|---|---|
| Origin $(0, 0, 0)$ | Nose tip apex (cones and rockets) or geometric center (sphere benchmark) |
| $+X$ axis (Axial / Roll) | Directed aft along the vehicle centerline (streamwise direction at $\alpha = 0^\circ$) |
| $+Y$ axis (Normal / Yaw-Pitch plane) | Directed upward in the longitudinal symmetry plane (normal force direction) |
| $+Z$ axis (Lateral / Pitch axis) | Completes the right-handed triad ($\hat{x} \times \hat{y} = \hat{z}$); pitch moment $M_z$ axis |
| State variables | $[\rho, u_x, u_y, u_z, p, T, e, k, \omega]$ |

**Freestream velocity decomposition** (from `src/physics/atmosphere.py`):

For a given flight Mach number $M_\infty$, altitude $h$, and angle of attack $\alpha$ in the $X$-$Y$ pitch plane (with sideslip $\beta_{slip} = 0^\circ$):

$$\mathbf{U}_\infty = \begin{bmatrix} U_\infty \cos\alpha \\ U_\infty \sin\alpha \\ 0 \end{bmatrix}, \quad U_\infty = M_\infty \, a_\infty$$

### 1.1 Computational domain topologies

Three distinct domain topologies are implemented via Jinja2 dictionary templates (`templates/openfoam/system/`):

| Mode | Template | Geometry Type | Domain Extent | Symmetry Assumption |
|---|---|---|---|---|
| **Mode 1: Supersonic Cone** | `blockMeshDict_cone.jinja` | $5^\circ$ axisymmetric wedge ($\theta_w = \pm 2.5^\circ$) | $x \in [-0.1, 0.5]$ m, $r \in [0, 0.4]$ m | Rotational periodicity ($\partial/\partial\phi = 0$); exact for $\alpha = 0^\circ$ |
| **Mode 2: Blunt Sphere** | `blockMeshDict_sphere.jinja` | $5^\circ$ polar O-Grid wedge around $R = 0.10$ m | Upstream $-3R$, radial $4R$, downstream wake $+4.5R$ | Axisymmetric detached bow shock and toroidal recirculation bubble |
| **Mode 3: 3D Rocket STL** | `blockMeshDict_rocket.jinja` | Full 3D half-tunnel ($z \ge 0$) or full-tunnel around STL | $x \in [-3L, +8L]$, $y \in [-4L, +4L]$, $z \in [0, 4L]$ | Longitudinal pitch symmetry across $Z = 0$ (`symmetryPlane`); preserves 3D turbulent wake asymmetry |

---

## 2. Units and physical constants

All internal calculations in `src/` and OpenFOAM dictionaries use **SI base units** (metres, kilograms, seconds, kelvin, pascals).

| Symbol | Value | Unit | Description | Source |
|---|---|---|---|---|
| $g_0$ | 9.80665 | m/s² | Standard sea-level gravitational acceleration | ISO 2533:1975 |
| $R_{\mathrm{air}}$ | 287.05287 | J/(kg·K) | Specific gas constant for dry air | ICAO / ISA 1976 |
| $M_{\mathrm{mol}}$ | 28.96 | kg/kmol | Mean molecular weight of dry air | OpenFOAM `specie` |
| $\gamma$ | 1.4 | — | Ratio of specific heats ($c_p / c_v$) for diatomic air | Calorically perfect gas |
| $c_p$ | 1004.5 | J/(kg·K) | Specific heat at constant pressure | $\gamma R_{\mathrm{air}} / (\gamma - 1)$ |
| $\mathrm{Pr}$ | 0.71 | — | Molecular Prandtl number for air | Standard air tables |
| $T_0$ | 288.15 | K | ISA sea-level standard temperature (15.0 °C) | ISO 2533:1975 |
| $p_0$ | 101 325.0 | Pa | ISA sea-level standard static pressure | ISO 2533:1975 |
| $\rho_0$ | 1.2250 | kg/m³ | ISA sea-level standard density | ISO 2533:1975 |
| $L_{\mathrm{tropo}}$ | $-0.0065$ | K/m | Tropospheric temperature lapse rate ($0\text{--}11$ km) | ISO 2533:1975 |
| $L_{\mathrm{strato}}$ | $+0.0010$ | K/m | Mid-stratospheric temperature lapse rate ($20\text{--}32$ km) | U.S. Standard Atmosphere 1976 |
| $\mu_{\mathrm{ref}}$ | $1.716 \times 10^{-5}$ | Pa·s | Sutherland reference dynamic viscosity at $T_{\mathrm{ref}}$ | White, *Viscous Fluid Flow* |
| $T_{\mathrm{ref}}$ | 273.15 | K | Sutherland reference temperature | White, *Viscous Fluid Flow* |
| $S$ | 110.4 | K | Sutherland effective temperature constant | White, *Viscous Fluid Flow* |
| $A_s$ | $1.458 \times 10^{-6}$ | kg/(m·s·K⁰·⁵) | OpenFOAM Sutherland transport coefficient | $\mu_{\mathrm{ref}} (T_{\mathrm{ref}} + S) / T_{\mathrm{ref}}^{1.5}$ |

---

## 3. Governing equations of compressible flow

The flow field is governed by the three-dimensional, unsteady, **compressible Navier-Stokes equations** in conservative integral form over a control volume $\Omega$ bounded by surface $\partial\Omega$:

$$\frac{\partial}{\partial t} \int_\Omega \mathbf{W} \, d\Omega + \oint_{\partial\Omega} \left( \mathbf{F}_c - \mathbf{F}_v \right) \cdot \hat{n} \, dS = 0$$

where the conservative state vector $\mathbf{W}$, convective flux tensor $\mathbf{F}_c$, and viscous flux tensor $\mathbf{F}_v$ are defined as:

$$\mathbf{W} = \begin{bmatrix} \rho \\ \rho \mathbf{U} \\ \rho E \end{bmatrix}, \quad \mathbf{F}_c = \begin{bmatrix} \rho \mathbf{U} \\ \rho \mathbf{U} \otimes \mathbf{U} + p \mathbf{I} \\ (\rho E + p)\mathbf{U} \end{bmatrix}, \quad \mathbf{F}_v = \begin{bmatrix} 0 \\ \boldsymbol{\tau} \\ \boldsymbol{\tau} \cdot \mathbf{U} - \mathbf{q} \end{bmatrix}$$

### 3.1 Thermodynamic closure (`hePsiThermo`)

The gas is modeled as a calorically perfect ideal gas (`perfectGas` + `hConst` in `thermophysicalProperties`):

$$p = \rho R_{\mathrm{air}} T = \frac{\rho}{\psi}, \quad \psi = \frac{1}{R_{\mathrm{air}} T}$$

where $\psi$ is the isothermal compressibility. Total specific energy $E$ relates to sensible internal energy $e = c_v T$ and kinetic energy via:

$$E = e + \frac{1}{2}|\mathbf{U}|^2, \quad h = e + \frac{p}{\rho} = c_p T$$

### 3.2 Viscous stress tensor and heat flux

Assuming a Newtonian fluid obeying Stokes' hypothesis, the viscous stress tensor $\boldsymbol{\tau}$ and Fourier heat flux $\mathbf{q}$ are:

$$\boldsymbol{\tau} = \mu_{\mathrm{eff}} \left[ \nabla\mathbf{U} + (\nabla\mathbf{U})^T - \frac{2}{3}(\nabla \cdot \mathbf{U})\mathbf{I} \right], \quad \mathbf{q} = -\kappa_{\mathrm{eff}} \nabla T$$

where $\mu_{\mathrm{eff}} = \mu + \mu_t$ and $\kappa_{\mathrm{eff}} = c_p \left( \frac{\mu}{\mathrm{Pr}} + \frac{\mu_t}{\mathrm{Pr}_t} \right)$ incorporate molecular and turbulent Eddy diffusivity.

---

## 4. Atmospheric model (ISA 1976) and freestream kinematics

Implemented in [`src/physics/atmosphere.py`](../src/physics/atmosphere.py) (`AtmosphereISA`), providing exact freestream boundary conditions from sea level ($h = 0$ m) up to the mid-stratosphere ($h = 32\,000$ m) across three geopotential layers:

### 4.1 Hydrostatic pressure and temperature integration

1. **Troposphere** ($0 \le h \le 11\,000$ m) — linear thermal decay ($L_{\mathrm{tropo}} = -0.0065$ K/m):

$$T(h) = T_0 + L_{\mathrm{tropo}} \, h, \quad p(h) = p_0 \left( \frac{T(h)}{T_0} \right)^{-\frac{g_0}{L_{\mathrm{tropo}} R_{\mathrm{air}}}}$$

2. **Lower Stratosphere / Tropopause** ($11\,000 < h \le 20\,000$ m) — isothermal layer ($T_{\mathrm{tropo}} = 216.65$ K):

$$T(h) = T_{\mathrm{tropo}}, \quad p(h) = p_{\mathrm{tropo}} \exp\left( -\frac{g_0 (h - 11\,000)}{R_{\mathrm{air}} T_{\mathrm{tropo}}} \right)$$

3. **Mid-Stratosphere** ($20\,000 < h \le 32\,000$ m) — positive thermal gradient ($L_{\mathrm{strato}} = +0.0010$ K/m):

$$T(h) = T_{\mathrm{tropo}} + L_{\mathrm{strato}} (h - 20\,000), \quad p(h) = p_{\mathrm{strato,1}} \left( \frac{T(h)}{T_{\mathrm{tropo}}} \right)^{-\frac{g_0}{L_{\mathrm{strato}} R_{\mathrm{air}}}}$$

### 4.2 Sutherland's viscosity law and derived freestream parameters

Molecular dynamic viscosity $\mu_\infty$ is evaluated as an explicit function of static temperature:

$$\mu(T) = \mu_{\mathrm{ref}} \left( \frac{T}{T_{\mathrm{ref}}} \right)^{3/2} \frac{T_{\mathrm{ref}} + S}{T + S} = \frac{A_s \, T^{3/2}}{T + T_s}$$

Acoustic speed, dynamic pressure, and isentropic stagnation temperature are then computed directly:

$$a_\infty = \sqrt{\gamma R_{\mathrm{air}} T_\infty}, \quad q_\infty = \frac{1}{2}\rho_\infty U_\infty^2 = \frac{1}{2}\gamma p_\infty M_\infty^2, \quad T_0 = T_\infty \left( 1 + \frac{\gamma - 1}{2}M_\infty^2 \right)$$

---

## 5. Shock-capturing numerical formulation (`rhoCentralFoam`)

Supersonic and transonic flows feature strong discontinuities (oblique shocks, detached bow shocks, contact surfaces) where classical pressure-based PISO/SIMPLE algorithms fail. The simulator orchestrates OpenFOAM's density-based solver **`rhoCentralFoam`**, which implements the **Kurganov and Tadmor (2000) central-upwind scheme**.

### 5.1 Kurganov central-upwind convective flux discretization

Unlike Riemann-solver-based Godunov methods (e.g., Roe or HLLC) that require costly characteristic matrix decompositions at every cell face, the Kurganov scheme integrates over the Riemann fan using exact one-sided local speeds of propagation.

At each internal cell face $f$ separating owner cell $P$ (left state $L$) and neighbour cell $N$ (right state $R$), the volumetric acoustic propagation speeds are:

$$a_{f,L} = \sqrt{\gamma R_{\mathrm{air}} T_L}, \quad a_{f,R} = \sqrt{\gamma R_{\mathrm{air}} T_R}$$

The directional maximum outward and inward wave speeds at the face are bounded by:

$$\phi_+ = \max\left( \mathbf{U}_L \cdot \mathbf{S}_f + a_{f,L}|\mathbf{S}_f|, \, \mathbf{U}_R \cdot \mathbf{S}_f + a_{f,R}|\mathbf{S}_f|, \, 0 \right)$$

$$\phi_- = \min\left( \mathbf{U}_L \cdot \mathbf{S}_f - a_{f,L}|\mathbf{S}_f|, \, \mathbf{U}_R \cdot \mathbf{S}_f - a_{f,R}|\mathbf{S}_f|, \, 0 \right)$$

Defining the weighting factor $\alpha_f = \frac{\phi_+}{\phi_+ - \phi_-}$, the numerical convective flux across face $f$ is evaluated as:

$$\mathbf{F}_{c,f} = \alpha_f \mathbf{F}_c(\mathbf{W}_L) + (1 - \alpha_f)\mathbf{F}_c(\mathbf{W}_R) + \omega_f \left( \mathbf{W}_R - \mathbf{W}_L \right)$$

where $\omega_f = \alpha_f (1 - \alpha_f)(\phi_+ - \phi_-)$ is the **artificial numerical diffusion coefficient** that automatically activates across steep shock gradients to prevent spurious Gibbs oscillations while vanishing in smooth regions.

### 5.2 TVD MUSCL reconstruction (`vanLeer`)

To achieve second-order spatial accuracy away from discontinuities while preserving monotonicity (Total Variation Diminishing property), cell-centered primitive variables $(\rho, \mathbf{U}, T, k, \omega)$ are interpolated to cell faces using the **van Leer flux limiter** (`fvSchemes.jinja`):

$$\Psi(r) = \frac{r + |r|}{1 + |r|}, \quad r = \frac{\phi_P - \phi_U}{\phi_N - \phi_P}$$

### 5.3 Acoustic Courant-Friedrichs-Lewy (CFL) condition

Explicit convective time-stepping (`ddtSchemes: Euler`) is stabilized by dynamically scaling $\Delta t$ (`adjustTimeStep: yes` in `controlDict.jinja`) to satisfy the compressible acoustic CFL bound:

$$\mathrm{Co} = \frac{\Delta t}{2 V_P} \sum_{f} \max\left(|\phi_+|, |\phi_-|\right) \le \mathrm{Co}_{\max} = 0.3$$

---

## 6. Turbulence closure (k-ω SST) and boundary conditions

For viscous rocket and sphere simulations, Reynolds-Averaged (Favre-averaged) Navier-Stokes (RANS) closure is achieved via **Menter's k-ω Shear Stress Transport (SST)** model (`turbulenceProperties.jinja`).

### 6.1 Transport equations

The SST formulation blends the robust near-wall $k$-$\omega$ model with the freestream-insensitive $k$-$\varepsilon$ model in the outer wake using a blending function $F_1$:

$$\frac{\partial (\rho k)}{\partial t} + \nabla \cdot (\rho \mathbf{U} k) = \tilde{P}_k - \beta^* \rho k \omega + \nabla \cdot \left[ (\mu + \sigma_k \mu_t)\nabla k \right]$$

$$\frac{\partial (\rho \omega)}{\partial t} + \nabla \cdot (\rho \mathbf{U} \omega) = \alpha \frac{\rho}{\mu_t}\tilde{P}_k - \beta \rho \omega^2 + \nabla \cdot \left[ (\mu + \sigma_\omega \mu_t)\nabla \omega \right] + 2(1 - F_1)\rho \sigma_{\omega 2}\frac{1}{\omega}\nabla k \cdot \nabla \omega$$

Crucially for supersonic shock-boundary-layer interaction (SBLI), the turbulent eddy viscosity $\mu_t$ is bounded by Bradshaw's structural stress limiter ($a_1 = 0.31$):

$$\mu_t = \frac{\rho a_1 k}{\max(a_1 \omega, \, S F_2)}$$

This prevents unphysical over-prediction of turbulent kinetic energy $k$ across normal and oblique shock waves.

### 6.2 Initial turbulent seeding and boundary conditions (`0/`)

Freestream turbulence intensity is initialized at $I = 1\%$ with turbulent viscosity ratio $\mu_t / \mu_\infty = 1$:

$$k_\infty = \frac{3}{2}(U_\infty I)^2, \quad \omega_\infty = \frac{\rho_\infty k_\infty}{\mu_\infty}$$

| Boundary Patch | Velocity `U` | Pressure `p` | Temperature `T` | Turbulence `k`, `omega`, `nut` |
|---|---|---|---|---|
| `inlet` / `freestream` | `fixedValue` ($\mathbf{U}_\infty$) | `fixedValue` ($p_\infty$) | `fixedValue` ($T_\infty$) | `fixedValue` ($k_\infty, \omega_\infty$) |
| `outlet` | `zeroGradient` / `inletOutlet` | `waveTransmissive` / `zeroGradient` | `zeroGradient` | `zeroGradient` |
| `body` (`nose`, `fuselage`, `fins`, `base`, `sphere`) | `noSlip` ($\mathbf{U} = \mathbf{0}$) or `slip` (inviscid) | `zeroGradient` ($\partial p/\partial n = 0$) | `zeroGradient` (adiabatic wall $T_w = T_{\mathrm{rec}}$) | `kqRWallFunction`, `omegaWallFunction`, `nutkWallFunction` |
| `wedge_front` / `wedge_back` | `wedge` | `wedge` | `wedge` | `wedge` |
| `symmetry` | `symmetryPlane` | `symmetryPlane` | `symmetryPlane` | `symmetryPlane` |

---

## 7. Analytical ground truth I — Supersonic cone (Taylor-Maccoll ODE)

Implemented in [`src/physics/taylor_maccoll.py`](../src/physics/taylor_maccoll.py) (`TaylorMaccollSolver`), providing an exact analytical benchmark for axisymmetric supersonic flow over a sharp cone of half-angle $\theta_c$ at $\alpha = 0^\circ$.

### 7.1 Conical self-similarity and governing ODE

In steady, inviscid, irrotational, axisymmetric supersonic flow over a semi-infinite cone, fluid properties are constant along rays emanating from the apex ($\partial/\partial r = 0$) and depend solely on the polar angle $\theta$. Normalizing velocities by the theoretical maximum adiabatic expansion speed $V_{\max} = \sqrt{2 c_p T_0}$:

$$V' = \frac{V}{V_{\max}} = \left( 1 + \frac{2}{(\gamma - 1)M^2} \right)^{-1/2}$$

Combining continuity, zero vorticity ($\nabla \times \mathbf{V} = \mathbf{0} \implies \frac{d V'_r}{d\theta} = V'_\theta$), and Euler momentum yields the non-linear second-order **Taylor-Maccoll ordinary differential equation**:

$$\frac{\gamma - 1}{2}\left( 1 - V_r^{\prime 2} - V_\theta^{\prime 2} \right)\left( 2 V'_r + V'_\theta \cot\theta + \frac{d V'_\theta}{d\theta} \right) - V_\theta^{\prime 2}\left( V'_r + \frac{d V'_\theta}{d\theta} \right) = 0$$

### 7.2 First-order state-space system and Rankine-Hugoniot initial values

Defining the state vector $\mathbf{y}(\theta) = [V_r', V_\theta']^T$ and dimensionless sound speed squared $\tilde{a}^2$:

$$\tilde{a}^2 = \frac{\gamma - 1}{2}\left(1 - V_r^{\prime 2} - V_\theta^{\prime 2}\right)$$

the system solved by `scipy.integrate.solve_ivp` (RK45) is:

$$\frac{d}{d\theta}\begin{bmatrix} V_r' \\ V_\theta' \end{bmatrix} = \begin{bmatrix} V_\theta' \\ \dfrac{V_\theta^{\prime 2} V_r' - \tilde{a}^2 (2 V_r' + V_\theta' \cot\theta)}{\tilde{a}^2 - V_\theta^{\prime 2}} \end{bmatrix}$$

Integration proceeds inward from a candidate oblique shock angle $\theta = \beta$ where the immediate post-shock state is fixed by the oblique Rankine-Hugoniot relations ($M_{n1} = M_\infty \sin\beta > 1$):

$$\tan\delta_s = 2\cot\beta \frac{M_{n1}^2 - 1}{M_\infty^2(\gamma + \cos 2\beta) + 2}, \quad M_{n2}^2 = \frac{1 + \frac{\gamma - 1}{2}M_{n1}^2}{\gamma M_{n1}^2 - \frac{\gamma - 1}{2}}, \quad M_2 = \frac{M_{n2}}{\sin(\beta - \delta_s)}$$

$$V'_r(\beta) = V'_2 \cos(\beta - \delta_s), \quad V'_\theta(\beta) = -V'_2 \sin(\beta - \delta_s)$$

### 7.3 Shooting method and Brent root-finding

The physical cone surface corresponds to the polar angle $\theta = \theta_c$ where flow tangency is satisfied ($V'_\theta(\theta_c) = 0$). A terminal event function detects the zero-crossing of $V'_\theta$, and **Brent's bracketed root-finding algorithm** (`scipy.optimize.brentq`, tolerance $10^{-6}$ rad) adjusts $\beta$ on the weak-shock branch until $\theta_{c,\mathrm{calc}}(\beta) - \theta_{c,\mathrm{target}} = 0$.

Once converged, the surface Mach number $M_c$ and exact surface-to-freestream static pressure ratio $p_c / p_\infty$ are recovered via isentropic relations behind the conical shock (since entropy is uniform behind a straight conical shock):

$$\frac{p_c}{p_\infty} = \left( \frac{2\gamma M_{n1}^2 - (\gamma - 1)}{\gamma + 1} \right) \left[ \frac{1 + \frac{\gamma - 1}{2}M_2^2}{1 + \frac{\gamma - 1}{2}M_c^2} \right]^{\frac{\gamma}{\gamma - 1}}$$

---

## 8. Analytical ground truth II — Blunt body bow shock (Rankine-Hugoniot, Rayleigh-Pitot, Billig)

Implemented in [`src/physics/blunt_body.py`](../src/physics/blunt_body.py) (`BluntBodySolver`), modeling supersonic flow over a sphere or hemispherical rocket nose cap of radius $R$ where the shock wave detaches upstream by a standoff distance $\Delta$.

### 8.1 Normal shock and Rayleigh supersonic Pitot formula

Along the stagnation streamline ($r = 0$), the detached bow shock is locally normal ($\beta = 90^\circ$). Across the normal shock, static pressure jumps and Mach number drops subsonically ($M_2 < 1$):

$$\frac{p_2}{p_\infty} = \frac{2\gamma M_\infty^2 - (\gamma - 1)}{\gamma + 1}, \quad M_2^2 = \frac{(\gamma - 1)M_\infty^2 + 2}{2\gamma M_\infty^2 - (\gamma - 1)}$$

The subsonic stream then decelerates isentropically from $M_2$ to zero velocity at the nose stagnation point. Combining the normal shock pressure jump with the post-shock isentropic compression yields the **Rayleigh Supersonic Pitot Formula** for the stagnation pressure ratio $p_{02}/p_\infty$:

$$\frac{p_{02}}{p_\infty} = \left[ \frac{(\gamma + 1)^2 M_\infty^2}{4\gamma M_\infty^2 - 2(\gamma - 1)} \right]^{\frac{\gamma}{\gamma - 1}} \left( \frac{1 - \gamma + 2\gamma M_\infty^2}{\gamma + 1} \right)$$

The theoretical maximum stagnation pressure coefficient $C_{p,\max}$ at the blunt nose tip is therefore:

$$C_{p,\max} = \frac{p_{02} - p_\infty}{\frac{1}{2}\gamma p_\infty M_\infty^2} = \frac{2}{\gamma M_\infty^2}\left( \frac{p_{02}}{p_\infty} - 1 \right)$$

At $M_\infty = 2.00$ ($\gamma = 1.4$), $p_{02}/p_\infty = 5.6404$ and $C_{p,\max} = 1.6573$.

### 8.2 Billig (1967) 3D hyperboloid bow shock correlation

The standoff distance $\Delta$ ahead of the nose tip and the shock vertex radius of curvature $R_c$ are computed via **Billig's experimental correlations** for spherical blunt bodies:

$$\frac{\Delta}{R} = 0.143 \exp\left( \frac{3.24}{M_\infty^2} \right), \quad \frac{R_c}{R} = 1.143 \exp\left( \frac{0.54}{(M_\infty - 1)^{1.2}} \right)$$

The three-dimensional bow shock surface of revolution is a hyperboloid asymptotic to the freestream Mach cone ($\mu_\infty = \arcsin(1/M_\infty)$):

$$x(r) = -\Delta + R_c \cot^2\mu_\infty \left( \sqrt{1 + \frac{r^2 \tan^2\mu_\infty}{R_c^2}} - 1 \right), \quad r = \sqrt{y^2 + z^2}$$

This analytical surface is meshed into a 360° `pyvista.PolyData` shell (`build_3d_bow_shock_surface`) and rendered in the 3D viewport as a translucent green reference wave overlaid directly onto the CFD solution.

---

## 9. Aerodynamic heating and recovery temperature model

Implemented in [`src/postprocess/aeroheating.py`](../src/postprocess/aeroheating.py). In high-speed compressible flight, kinetic energy dissipation via viscous shear inside the boundary layer heats the vehicle skin even under adiabatic wall conditions ($\partial T / \partial n|_w = 0$).

### 9.1 Stagnation vs. adiabatic recovery temperature

At the stagnation point ($U_e = 0$), the gas reaches the full **isentropic stagnation temperature** $T_0$ (since a stationary normal shock is isoenergetic, $h_{02} = h_{0\infty} \implies T_{02} = T_0$ for a calorically perfect gas):

$$T_0 = T_\infty \left( 1 + \frac{\gamma - 1}{2}M_\infty^2 \right)$$

Along the downstream conical or cylindrical surfaces where the boundary-layer edge Mach number is $M_e > 0$, imperfect balance between viscous shear work and thermal conduction ($\mathrm{Pr} = 0.71 < 1$) results in an **adiabatic wall recovery temperature** $T_{\mathrm{rec}} < T_0$:

$$T_{\mathrm{rec}} = T_e \left( 1 + r \frac{\gamma - 1}{2}M_e^2 \right)$$

where the **thermal recovery factor** $r$ depends on the boundary-layer regime (Eckert / van Driest):
- **Laminar boundary layer:** $r_{\mathrm{lam}} = \mathrm{Pr}^{1/2} = \sqrt{0.71} \approx 0.8426$
- **Turbulent boundary layer:** $r_{\mathrm{turb}} = \mathrm{Pr}^{1/3} = (0.71)^{1/3} \approx 0.8921$

At sea level ($T_\infty = 288.15$ K) and $M_\infty = 2.00$:
- $T_0 = 518.67\text{ K}$ (245.5 °C)
- $T_{\mathrm{rec,turb}} = 493.80\text{ K}$ (220.7 °C)
- $T_{\mathrm{rec,lam}} = 482.39\text{ K}$ (209.2 °C)

---

## 10. Meshing topologies and boundary layer resolution

### 10.1 Polar O-Grid mesh (`blockMeshDict_sphere.jinja`)

For the supersonic sphere benchmark ($R = 0.10$ m), a multi-block structured **polar O-Grid wedge mesh** ($8\,400$ hexahedral cells) is generated around the circular arc:
- **Upstream bow shock block:** $60 \times 70$ radial-polar cells graded toward the wall (`simpleGrading (4 1 1)`) to capture both the detached shock at $x \approx -0.032$ m and the stagnation boundary layer.
- **Downstream wake block:** $70 \times 70$ cells extending to $x = +0.45$ m ($4.5R$) to resolve the supersonic expansion fan, boundary-layer separation, and recompression wake shock.

### 10.2 Full 3D wind tunnel and `snappyHexMesh` (`snappyHexMeshDict.jinja`)

For complex 3D rocket geometries (`nose.stl`, `fuselage.stl`, `fins.stl`, `base.stl`), an octree-based hex-dominant mesh is constructed in three automated stages:

1. **Castellated Mesh (`castellatedMesh: true`):**
   - Feature edge extraction (`surfaceFeatureExtractDict.jinja`) at $150^\circ$ included angle.
   - Surface refinement levels $(3, 4)$ on fuselage/base and $(4, 5)$ on nose tip and fin leading/trailing edges.
   - Volumetric refinement box (`wakeAndShockBox`, level 3) spanning $x \in [-0.5L, +4.5L]$ to preserve shock sharpness and capture asymmetric turbulent vortex shedding in the rocket wake.
2. **Surface Snapping (`snap: true`):**
   - Projects jagged staircase vertices onto the triangulated STL patches (`nSmoothPatch: 3`, `nRelaxIter: 5`).
3. **Boundary Layer Inflation (`addLayers: true`):**
   - Extrudes $N_{\mathrm{layers}} = 4$ prism/hexahedral inflation layers normal to all rocket walls with expansion ratio $1.25$ and first-layer relative thickness $0.3$, ensuring strict quality thresholds (`maxNonOrtho: 65`, `maxBoundarySkewness: 20`).

---

## 11. Aerodynamic forces, center of pressure, and GNC export (`AeroDB`)

### 11.1 OpenFOAM `forceCoeffs` integration

Aerodynamic forces $\mathbf{F} = \mathbf{F}_p + \mathbf{F}_v$ and moments $\mathbf{M}$ about the reference center $\mathbf{x}_{\mathrm{ref}} = (X_{\mathrm{ref}}, 0, 0)$ are integrated at every solver timestep over the vehicle patches:

$$\mathbf{F} = \oint_{S_{\mathrm{body}}} \left( (p - p_\infty)\hat{n} - \boldsymbol{\tau} \cdot \hat{n} \right) dS, \quad \mathbf{M} = \oint_{S_{\mathrm{body}}} (\mathbf{x} - \mathbf{x}_{\mathrm{ref}}) \times \left( (p - p_\infty)\hat{n} - \boldsymbol{\tau} \cdot \hat{n} \right) dS$$

Non-dimensional aerodynamic coefficients are normalized by freestream dynamic pressure $q_\infty$, reference cross-sectional area $A_{\mathrm{ref}}$, and reference length $L_{\mathrm{ref}}$:

$$C_D = \frac{\mathbf{F} \cdot \hat{e}_{\mathrm{drag}}}{q_\infty A_{\mathrm{ref}}}, \quad C_L = \frac{\mathbf{F} \cdot \hat{e}_{\mathrm{lift}}}{q_\infty A_{\mathrm{ref}}}, \quad C_{m,\mathrm{pitch}} = \frac{M_z}{q_\infty A_{\mathrm{ref}} L_{\mathrm{ref}}}$$

where the wind-axis unit vectors rotate with the angle of attack $\alpha$:

$$\hat{e}_{\mathrm{drag}} = \begin{bmatrix} \cos\alpha \\ \sin\alpha \\ 0 \end{bmatrix}, \quad \hat{e}_{\mathrm{lift}} = \begin{bmatrix} -\sin\alpha \\ \cos\alpha \\ 0 \end{bmatrix}$$

### 11.2 Center of pressure ($X_{cp}$) for GNC stability

Implemented in [`src/postprocess/aerodb_exporter.py`](../src/postprocess/aerodb_exporter.py). Transforming wind-axis coefficients $(C_L, C_D)$ into the body-normal force coefficient $C_N$:

$$C_N = C_L \cos\alpha + C_D \sin\alpha$$

The longitudinal location of the **center of pressure** $X_{cp}$ measured aft from the nose tip ($x = 0$) is obtained from moment equilibrium:

$$X_{cp} = X_{\mathrm{ref}} - \frac{C_{m,\mathrm{pitch}}}{C_N} L_{\mathrm{ref}}$$

A rocket is aerodynamically statically stable in pitch when $X_{cp} > X_{\mathrm{cg}}$, giving a positive static margin $\mathrm{SM} = (X_{cp} - X_{\mathrm{cg}}) / D_{\mathrm{ref}} > 0$. Converged flight points are exported to structured JSON and appended to the master CSV database (`AeroDB`) for direct ingestion by 6-DOF flight dynamics and trajectory simulation tools.

---

## 12. Post-processing: Numerical Schlieren and Lagrangian particle tracking

### 12.1 Numerical Schlieren imaging (`src/postprocess/schlieren.py`)

Optical wind-tunnel Schlieren photography visualizes refractive index variations proportional to the first spatial derivative of gas density (Gladstone-Dale relation). Numerically, the simulator computes the spatial gradient of the cell-centered density field $\nabla\rho = (\partial\rho/\partial x, \partial\rho/\partial y, \partial\rho/\partial z)$ using VTK's green-gauss derivative filter (`compute_derivative`), and evaluates both the Euclidean magnitude and a logarithmic contrast enhancement field:

$$S(\mathbf{x}) = |\nabla\rho| = \sqrt{\left(\frac{\partial\rho}{\partial x}\right)^2 + \left(\frac{\partial\rho}{\partial y}\right)^2 + \left(\frac{\partial\rho}{\partial z}\right)^2}, \quad S_{\log}(\mathbf{x}) = \log_{10}\left( 1 + |\nabla\rho| \right)$$

For axisymmetric wedge cases, the $Z = 0$ slice is mirrored across the longitudinal symmetry axis ($Y \to -Y$) to render a full-field bilateral Schlieren photograph.

### 12.2 Residence-time Lagrangian particle system (`src/postprocess/particles.py`)

To visualize real-time deceleration across shock waves and recirculation inside wake vortices at 60 FPS without per-frame VTK actor reconstruction overhead:

1. **RK45 Streamline Integration:** Seed points are integrated through the steady velocity field $\mathbf{U}(\mathbf{x})$ using 4th/5th-order Runge-Kutta (`streamlines_from_source`).
2. **Physical Residence Time Parameterization:** Along each trajectory of arc-length $s$, the physical residence time $t(s)$ is integrated from the local speed $|\mathbf{U}(s)|$:

$$t(s) = \int_0^s \frac{ds'}{\max\left(|\mathbf{U}(s')|, U_{\min}\right)}$$

3. **360° Volumetric Revolution & In-Place GPU Update:** Trajectories are replicated across 16 diametral planes (32 radial half-planes around the $X$ axis). At each 60 Hz animation tick ($\Delta\tau$), particle phase offsets advance uniformly in physical time $t$, and the single `pyvista.PolyData` vertex buffer (`cloud_polydata.points`) and scalar Mach buffer (`cloud_polydata.point_data["Mach"]`) are updated in-place in GPU memory. Particles naturally bunch up and slow down as they traverse the detached bow shock and base wake.

---

## 13. Asynchronous execution architecture and real-time telemetry

To prevent GUI freezing during multi-minute OpenFOAM runs in WSL2, the software enforces a strict **Model-View-Worker concurrency architecture**:

1. **Subprocess Bridge (`src/backend/runner.py`):** Translates Windows paths (`C:\...`) to POSIX WSL2 mount points (`/mnt/c/...`) and spawns `wsl.exe -d Ubuntu-22.04` executing `blockMesh`, `surfaceFeatureExtract`, `snappyHexMesh`, and `rhoCentralFoam` with unbuffered stdout (`PYTHONUNBUFFERED=1`).
2. **Streaming Regex Parser (`src/backend/log_parser.py`):** Parses solver stdout line-by-line in real time, extracting:
   - Time step $t$ (`^Time = ([0-9.eE+-]+)`)
   - Mean and maximum Courant numbers (`Courant Number mean: ... max: ...`)
   - Initial and final linear-solver residuals for $p, U_x, U_y, U_z, e, k, \omega$
   - Integrated aerodynamic coefficients $C_D, C_L, C_{m,\mathrm{pitch}}$
3. **Qt Signal Telemetry (`PyQt6` + `PyQtGraph`):** Emits thread-safe Qt signals (`residual_parsed`, `forces_parsed`) to update hardware-accelerated `PyQtGraph` plots in real time, as well as reconstructing complete historical convergence curves from `log.rhoCentralFoam` and `postProcessing/forceCoeffs/0/coefficient.dat` when loading pre-computed benchmark cases.

---

## 14. Validation summary

All physical modules, analytical solvers, and OpenFOAM CFD benchmarks are verified by an automated test suite (`25 passed` in `tests/`):

| Physical Quantity / Benchmark | Simulator / CFD Result | Reference / Analytical Value | Source | Relative Error |
|---|---|---|---|---|
| **ISA Sea-Level State** ($h = 0$ m) | $T = 288.15$ K, $p = 101\,325$ Pa, $\rho = 1.2250$ kg/m³ | Exact match | ISO 2533:1975 | `< 0.001%` |
| **ISA Tropopause State** ($h = 11$ km) | $T = 216.65$ K, $p = 22\,632.1$ Pa, $\rho = 0.3639$ kg/m³ | $216.65$ K, $22\,632$ Pa, $0.3639$ kg/m³ | ICAO Standard Atmosphere | `< 0.01%` |
| **ISA Mid-Stratosphere** ($h = 20$ km) | $T = 216.65$ K, $p = 5\,474.9$ Pa, $\rho = 0.0880$ kg/m³ | $216.65$ K, $5\,474.9$ Pa, $0.0880$ kg/m³ | U.S. Standard Atmosphere 1976 | `< 0.01%` |
| **Taylor-Maccoll Solver** ($\theta_c = 15^\circ, M_\infty = 2.0$) | $\beta = 33.915^\circ$, $p_c/p_\infty = 1.5699$, $M_c = 1.7069$ | $\beta = 33.92^\circ$, $p_c/p_\infty = 1.570$ | Sims (1964), NASA SP-3004 | `< 0.02%` |
| **Taylor-Maccoll Solver** ($\theta_c = 20^\circ, M_\infty = 3.0$) | $\beta = 29.882^\circ$, $p_c/p_\infty = 2.6177$, $M_c = 2.2661$ | $\beta = 29.88^\circ$, $p_c/p_\infty = 2.618$ | Anderson, *Mod. Comp. Flow* | `< 0.02%` |
| **OpenFOAM Cone Shock Angle** ($\theta_c = 15^\circ, M_\infty = 2.0$) | $\beta_{\mathrm{CFD}} = 34.387^\circ$ | $\beta_{\mathrm{TM}} = 33.915^\circ$ | Exact Taylor-Maccoll ODE | **`1.39%`** |
| **OpenFOAM Cone Surface Pressure** ($p_c / p_\infty$) | $1.5701$ | $1.5699$ | Exact Taylor-Maccoll ODE | **`0.01%`** |
| **OpenFOAM Cone Surface Mach** ($M_c$) | $1.7146$ | $1.7069$ | Exact Taylor-Maccoll ODE | **`0.45%`** |
| **Blunt Body Rayleigh-Pitot** ($M_\infty = 2.0$) | $p_{02}/p_\infty = 5.6404$, $C_{p,\max} = 1.6573$ | $p_{02}/p_\infty = 5.640$, $C_{p,\max} = 1.657$ | NACA Report 1135 | `< 0.01%` |
| **Billig Sphere Standoff** ($M_\infty = 2.0, R = 0.10$ m) | $\Delta = 0.03215$ m ($\Delta/R = 0.3215$) | $\Delta/R \approx 0.32$ | Billig (1967), *J. Spacecraft* | Consistent |
| **OpenFOAM Sphere Drag** ($M_\infty = 2.0$, viscous) | $C_D \approx 0.96$ (converged steady state) | $C_D \approx 0.92\text{ to }1.00$ | Hoerner, *Fluid-Dynamic Drag* | Within experimental band |

---

## 15. Known limitations and deliberate simplifications

| Limitation | Physical Impact | Engineering Rationale |
|---|---|---|
| **Calorically perfect gas** ($\gamma = 1.4, c_p = \text{const}$) | Over-predicts stagnation temperature $T_0$ at hypersonic speeds ($M_\infty > 5$) where vibrational excitation and $\mathrm{O}_2/\mathrm{N}_2$ dissociation absorb thermal energy | Valid and exact across the target subsonic, transonic, and supersonic launch envelope ($M_\infty \in [0.2, 3.5+]$) where $T_0 < 800$ K |
| **5° axisymmetric wedge for Cone & Sphere benchmarks** | Enforces circumferential symmetry ($\partial/\partial\phi = 0, u_\phi = 0$); prevents 3D helical vortex shedding in the sphere wake | Reduces cell count from $\sim 1.5\times 10^6$ to $\sim 10^4$ cells, enabling sub-minute benchmark verification against analytical 2D/axisymmetric theory |
| **Adiabatic wall thermal boundary condition (`zeroGradient`)** | Computes maximum steady-state recovery temperature $T_{\mathrm{rec}}$ rather than transient conductive heat transfer into the rocket skin | Provides a conservative upper-bound thermal design envelope without requiring structural material thickness or internal thermal capacitance specs |
| **Cold-wall base flow (no rocket exhaust plume)** | Base drag ($C_{D,\mathrm{base}}$) is evaluated for an unpowered coast / power-off condition | Power-on plume interaction requires multi-species combustion products; power-off base drag represents the worst-case coasting drag post-MECO |
| **RANS turbulence closure** ($k$-$\omega$ SST) | Models ensemble-averaged Reynolds stresses rather than resolving unsteady large-scale turbulent eddies (LES/DES) | Standard industrial aerospace methodology for aerodynamic database generation, balancing boundary-layer separation accuracy with computational tractability |

---

## 16. References

1. **Anderson, J. D.** *Modern Compressible Flow: With Historical Perspective*, 4th ed. McGraw-Hill, 2021. — Oblique shock relations, Taylor-Maccoll conical flow ODE, normal shock and Rayleigh-Pitot derivations.
2. **Kurganov, A., and Tadmor, E.** "New High-Resolution Central Schemes for Nonlinear Conservation Laws and Convection-Diffusion Equations," *Journal of Computational Physics*, Vol. 160, No. 1, 2000, pp. 241–282. — Mathematical foundation of OpenFOAM's `rhoCentralFoam` solver.
3. **Greenshields, C. J., Weller, H. G., Gasparini, L., and Reese, J. M.** "Implementation of Semi-Discrete, Non-Staggered Central Schemes in a Colocated, Polyhedral, Finite Volume Framework, for High-Speed Viscous Flows," *International Journal for Numerical Methods in Fluids*, Vol. 63, No. 1, 2010, pp. 1–21. — Official validation paper of `rhoCentralFoam`.
4. **Billig, F. S.** "Shock-Wave Shapes Around Spherical- and Cylindrical-Nosed Bodies," *Journal of Spacecraft and Rockets*, Vol. 4, No. 6, 1967, pp. 822–823. — Empirical standoff distance $\Delta$ and hyperboloid bow shock equations.
5. **Menter, F. R.** "Two-Equation Eddy-Viscosity Turbulence Models for Engineering Applications," *AIAA Journal*, Vol. 32, No. 8, 1994, pp. 1598–1605. — Formulation of the $k$-$\omega$ SST turbulence model.
6. **White, F. M.** *Viscous Fluid Flow*, 3rd ed. McGraw-Hill, 2006. — Sutherland's viscosity law, compressible boundary-layer recovery factor ($r = \mathrm{Pr}^{1/2}, \mathrm{Pr}^{1/3}$), and aerodynamic heating.
7. **Sims, J. L.** *Tables for Supersonic Flow Around Right Circular Cones at Zero Angle of Attack*, NASA SP-3004, 1964. — Reference benchmark tables for Taylor-Maccoll validation.
8. **Ames Research Staff.** *Equations, Tables, and Charts for Compressible Flow*, NACA Report 1135, 1953. — Normal shock, oblique shock, and Rayleigh-Pitot reference data.
9. **International Organization for Standardization.** *Standard Atmosphere*, ISO 2533:1975 / U.S. Standard Atmosphere, 1976. — Geopotential altitude layers, lapse rates, and sea-level reference constants.
