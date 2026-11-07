# ROADMAP Y TO-DO LIST OPERATIVO: SIMULADOR CFD DE COHETES

Este documento establece la ruta crítica secuencial de desarrollo. Para maximizar la probabilidad de éxito y evitar reescrituras de código, se aplica el principio de ingeniería de software **"Headless-First, Physics-First"**: la física y la automatización por consola deben funcionar y estar validadas al 100% antes de construir la interfaz gráfica.

---

## 🗺️ Mapa de Fases (Ruta Crítica)

```
[Fase 1: Infraestructura] ──> [Fase 2: Física ISA y Analítica] ──> [Fase 3: Backend OpenFOAM Headless]
                                                                                │
                                                                                ▼
[Fase 6: Orquestador Asíncrono] <── [Fase 5: Geometría Cohete y Malla] <── [Fase 4: Benchmark Cono Validado]
          │
          ▼
[Fase 7: Algoritmos Postproceso] ──> [Fase 8: GUI PyQt6 + PyVista] ──> [Fase 9: AeroDB y Portfolio]
```

---

## 📋 Lista Maestra de Tareas (To-Do Checklist)

### Fase 1: Infraestructura y Blindaje del Hardware
> **Objetivo:** Preparar WSL2 y el entorno Python para evitar cuelgues por falta de RAM (OOM Killer) y habilitar aceleración GPU.

- [x] `1.1` **Crear archivo `.wslconfig` en Windows (`%USERPROFILE%\.wslconfig`)**
  - Asignar límite de memoria RAM a WSL2: `memory=10GB`, `swap=8GB` y `processors=12`.
- [x] `1.2` **Instalar OpenFOAM v2312 en Ubuntu WSL2**
  - Distribución oficial OpenCFD instalada (`openfoam2312-default`).
  - Entorno configurado en `/root/.bashrc` y comprobado con `rhoCentralFoam -help`.
  - Verificados binarios clave: `blockMesh`, `snappyHexMesh`, `surfaceFeatureExtract`, `checkMesh` y `mpirun`.
- [x] `1.3` **Configurar entorno virtual Python**
  - Entorno virtual creado en `venv/` con Python 3.13 y `.gitignore` configurado.
- [x] `1.4` **Instalar dependencias de cálculo y renderizado**
  - Paquetes instalados: `PyQt6`, `PyVista 0.49`, `VTK 9.7`, `PyQtGraph`, `Jinja2`, `SciPy`, `NumPy`, `Matplotlib`, `pytest`.
- [x] `1.5` **Verificar aceleración por hardware de la GPU con PyVista**
  - Verificado: PyVista detecta aceleración gráfica. (Nota: por defecto toma Intel UHD; configuraremos variable para forzar la RTX 4060).

---

### Fase 2: Módulo de Física de Vuelo (Atmósfera ISA 1976)
> **Objetivo:** Crear el motor matemático que transforma `(Mach, Altitud, Alfa)` en variables termodinámicas y cinemáticas exactas.

- [x] `2.1` **Implementar ecuaciones atmosféricas ISA 1976 (`src/physics/atmosphere.py`)**
  - Programados gradientes para troposfera (0 a 11 km) y estratosfera (11 a 32 km).
  - Cálculo de presión $p_\infty$, temperatura $T_\infty$ y densidad $\rho_\infty$.
- [x] `2.2` **Calcular viscosidad y velocidad del sonido**
  - Ley de Sutherland para $\mu_\infty(T)$.
  - Cálculo de $a = \sqrt{\gamma R T_\infty}$ y velocidad de flujo $U_\infty = \text{Mach} \cdot a$.
- [x] `2.3` **Proyectar el vector velocidad según ángulo de ataque $\alpha$**
  - Componentes $U_x = U_\infty \cos(\alpha)$, $U_y = U_\infty \sin(\alpha)$, $U_z = 0$.
- [x] `2.4` **Crear suite de pruebas unitarias (`tests/test_atmosphere.py`)**
  - Suite de 7 tests pasando al 100% validados frente a datos oficiales de la NASA/OACI.

---

### Fase 3: Automatización de OpenFOAM (Backend Headless)
> **Objetivo:** Generar casos de simulación compresibles de forma 100% programática usando plantillas Jinja2.

- [x] `3.1` **Crear estructura de carpetas para plantillas (`templates/openfoam/`)**
  - Creadas carpetas `0.orig/`, `constant/` y `system/`.
- [x] `3.2` **Escribir plantilla Jinja2 de `thermophysicalProperties`**
  - Configurado `hePsiThermo`, mezcla pura `pureMixture`, gas perfecto `hConst` y Sutherland para aire v2312.
- [x] `3.3` **Escribir plantilla Jinja2 de `turbulenceProperties`**
  - Configurado `simulationType RAS;` con modelo `kOmegaSST`.
- [x] `3.4` **Escribir plantilla Jinja2 de `fvSchemes` y `fvSolution` para `rhoCentralFoam`**
  - Esquemas Kurganov con limitadores de flujo `vanLeer` para captura de choques.
  - Solvers `diagonal` para variables hiperbólicas y `smoothSolver` para transporte.
- [x] `3.5` **Escribir plantilla Jinja2 de `controlDict`**
  - Control de paso de tiempo variable con CFL < 0.75 y función `forceCoeffs` desglosada por parches.
- [x] `3.6` **Escribir plantilla Jinja2 de condiciones de contorno (`0.orig/`)**
  - Plantillas para `U`, `p`, `T`, `k`, `omega` con pared adiabática (`zeroGradient` para recuperación de calor $T_{rec}$).
- [x] `3.7` **Crear generador de casos en Python (`src/backend/case_generator.py`)**
  - Clase `OpenFOAMCaseGenerator` testeada y validada con 1 test unitario passing.

---

### Fase 4: Modo de Validación Física (Benchmark Taylor-Maccoll)
> **Objetivo:** Probar el solver y el generador de casos contra una solución exacta antes de simular el cohete complejo.

- [x] `4.1` **Implementar solver analítico de choque cónico exacto (`src/physics/taylor_maccoll.py`)**
  - EDOs no lineales de Taylor-Maccoll resueltas con `solve_ivp` y búsqueda acotada de la solución de choque débil con `brentq`.
  - Tests unitarios pasando al 100% (`tests/test_taylor_maccoll.py`) con ángulos de choque y saltos de presión exactos de NACA 1135.
- [ ] `4.2` **Crear caso benchmark de cono supersónico en OpenFOAM**
  - Crear cono analítico en STL o con `blockMesh`.
  - Mallar con refinamiento zonal.
- [ ] `4.3` **Ejecutar simulación del cono con `rhoCentralFoam`**
  - Correr caso en WSL2 hasta convergencia de residuos.
- [ ] `4.4` **Crear script de extracción y comparación de choque (`src/validation/compare_benchmark.py`)**
  - Extraer gradiente de densidad $|\nabla \rho|$ y presión en pared de los resultados CFD.
  - Calcular la desviación porcentual respecto a la teoría de Taylor-Maccoll (debe ser $< 3\%$).

---

### Fase 5: Preparación de Geometría del Cohete y Mallado 3D
> **Objetivo:** Importar el diseño de SolidWorks por piezas y generar una malla simétrica de alta calidad en menos de 2.5 millones de celdas.

- [ ] `5.1` **Modelar cohete en SolidWorks y aplicar corte longitudinal**
  - Generar el modelo del cohete y seccionar longitudinalmente por el plano $Y=0$ para obtener exactamente la mitad ($1/2$).
- [ ] `5.2` **Exportar parches STL independientes en formato ASCII**
  - Exportar por separado: `nose.stl`, `fuselage.stl`, `fins.stl`, `base.stl`.
  - Verificar que las coordenadas coinciden en el mismo origen espacial y la escala está en metros.
- [ ] `5.3` **Diseñar dominio exterior en `blockMeshDict`**
  - Crear caja de túnel numérico: entrada a $5L$ aguas arriba, salida a $10L$ aguas abajo, radio $5L$.
  - Cara longitudinal configurada como parche tipo `symmetryPlane`.
- [ ] `5.4` **Configurar `snappyHexMeshDict` con refinamiento selectivo**
  - Definir caja de refinamiento (refinement box) envolviendo el cohete y la estela.
  - Refinamiento nivel alto en `nose` y bordes de ataque de `fins`.
  - Refinamiento medio en `fuselage`.
  - Configurar 3 a 5 capas prismáticas de inflación (`addLayers`) con factor de expansión $1.2$.
- [ ] `5.5` **Validar calidad de malla (`checkMesh`)**
  - Verificar que el número de celdas se mantiene entre $1.5\text{M}$ y $2.8\text{M}$ para respetar los 16 GB de RAM.
  - Comprobar que `Max skewness < 4` y no hay caras no ortogonales graves.

---

### Fase 6: Orquestador Asíncrono y Streaming de Residuos
> **Objetivo:** Conectar Python con WSL2 sin bloquear la interfaz, leyendo la salida del solver en tiempo real.

- [x] `6.1` **Implementar worker en segundo plano (`src/backend/runner.py`)**
  - Clase `OpenFOAMRunner` basada en `QThread` que ejecuta comandos en WSL2 (`wsl.exe -d Ubuntu-22.04`).
  - Redirección de stdout en tiempo real y soporte para cancelación segura (`cancel()`).
  - Conversor bidireccional de rutas Windows/WSL (`windows_to_wsl_path`).
- [x] `6.2` **Crear parser de logs en tiempo real (`src/backend/log_parser.py`)**
  - Implementado con Regex optimizadas para capturar `Time`, residuos (`Ux`, `Uy`, `e`, `k`, `omega`), `Courant max` y datos de `forceCoeffs` (`Cd`, `Cl`, `Cm`).
  - Suite de tests unitarios pasando al 100% (`tests/test_log_parser.py`).
- [x] `6.3` **Diseñar sistema de Señales y Slots en el Worker**
  - Señales Qt: `status_changed`, `log_line_received`, `iteration_updated(IterationMetrics)`, `force_updated(dict)`, `progress_changed(float)`, `finished_success(float)`, `error_occurred(str)`.
  - Verificado con tests unitarios (`tests/test_runner.py`).

---

### Fase 7: Algoritmos de Post-procesado y Renderizado 3D
> **Objetivo:** Extraer los campos físicos clave (ondas de choque, temperatura de fricción y líneas de corriente).

- [ ] `7.1` **Implementar filtro de efecto Schlieren numérico**
  - En PyVista/VTK, calcular el gradiente del campo escalar de densidad: $|\nabla \rho| = \sqrt{(\partial\rho/\partial x)^2 + (\partial\rho/\partial y)^2 + (\partial\rho/\partial z)^2}$.
  - Renderizar en escala de grises invertida para visualizar ondas de choque oblicuas y expansiones de Prandtl-Meyer.
- [ ] `7.2` **Mapear temperatura de pared (Aeroheating)**
  - Extraer el parche de superficie del cohete (`nose`, `fuselage`, `fins`).
  - Aplicar mapa de color térmico (inferno / coolwarm) mostrando la temperatura de recuperación $T_{rec}$ en Kelvin.
  - Localizar y etiquetar el punto de remanso (stagnation point) en la punta de la ojiva ($T_{stagnation} = T_\infty(1 + \frac{\gamma-1}{2}M^2)$).
- [ ] `7.3` **Generar líneas de corriente 3D (Streamlines)**
  - Configurar semillero de partículas en la entrada (`source_rake` o plano).
  - Trazar integrales del campo vectorial de velocidad $U$ coloreadas por número de Mach local.
- [ ] `7.4` **Crear generador de informes AeroDB (`src/postprocess/aerodb_exporter.py`)**
  - Exportar resultados convergidos a un archivo `results.json` y `aerodb.csv` con: Mach, Altitud, $\alpha$, $C_D$, $C_L$, $C_m$, $x_{cp}$ (centro de presiones) y $T_{max,nose}$.

---

### Fase 8: Interfaz Gráfica de Usuario (PyQt6 + PyVista)
> **Objetivo:** Construir una aplicación de escritorio profesional con visor 3D interactivo y telemetría de convergencia.

- [ ] `8.1` **Crear estructura modular de la GUI (`src/gui/`)**
  - `main_window.py`: Ventana principal, menús y distribución de paneles.
  - `control_panel.py`: Barra lateral izquierda de configuración.
  - `plots_panel.py`: Panel inferior con pestañas para gráficas en vivo.
  - `viewport_panel.py`: Área central con el visor 3D.
- [ ] `8.2` **Construir el panel de control de vuelo**
  - Spinboxes para `Mach` (rango 0.2 a 4.0), `Altitud (km)` (0 a 30 km) y `Ángulo de ataque` ($0^\circ$ a $15^\circ$).
  - Selector de archivos para importar los STL del cohete.
  - Botón selector entre **"Simulación Completa"** y **"Modo Validación (Cono Teórico)"**.
  - Display informativo con las variables ISA calculadas en tiempo real al mover los controles.
  - Botón "Iniciar Simulación" y botón "Abortar".
- [ ] `8.3` **Integrar gráficas de convergencia en vivo con PyQtGraph**
  - Gráfica 1: Residuos logarítmicos vs. Iteraciones/Tiempo físico.
  - Gráfica 2: Coeficiente de resistencia ($C_D$) y sustentación ($C_L$) estabilizándose.
- [ ] `8.4` **Embeber el visor 3D de PyVista (`QVTKRenderWindowInteractor`)**
  - Integrar el widget interactivo en el centro de la interfaz.
  - Configurar iluminación especular para materiales aerodinámicos.
  - Añadir barra de herramientas 3D: conmutador para ver Malla, Geometría sólida, Presión en pared, Temperatura de recuperación, Schlieren en plano simétrico y Líneas de corriente.
  - Añadir leyenda de colores interactiva y ejes coordenados 3D.

---

### Fase 9: Pulido Industrial, Validación Final y Presentación
> **Objetivo:** Generar el caso de estudio definitivo para captar la atención de reclutadores en el sector aeroespacial.

- [ ] `9.1` **Correr una campaña de vuelo completa (Mach Sweep)**
  - Ejecutar la simulación para $M = 0.5, 0.9, 1.2, 2.0, 3.0$.
  - Generar la curva transónica crítica de resistencia ($C_D$ vs. Mach) demostrando el aumento brusco en la barrera del sonido (Mach drag divergence).
- [ ] `9.2` **Automatizar capturas de alta definición con `pvpython`**
  - Script que guarde renders automáticos 4K de las ondas de choque y mapas de temperatura para informes técnicos.
- [ ] `9.3` **Redactar documentación técnica y README del repositorio**
  - Explicar la arquitectura software (diagrama de bloques).
  - Incluir la tabla de validación analítica frente a Taylor-Maccoll (% de error).
  - Añadir GIFs/vídeos de la interfaz funcionando con el cohete a Mach 3.0.

---

## 🎯 Próximo Paso Inmediato

Para arrancar hoy mismo sin bloqueos, la primera tarea física a ejecutar es:
👉 **Tarea 1.1 y 1.2: Configurar `.wslconfig` e instalar OpenFOAM en WSL2.**

¿Quieres que te redacte el contenido exacto de `.wslconfig` y los comandos de terminal de WSL2 para completar la **Fase 1**, o prefieres empezar programando la clase del **Módulo ISA 1976 (`src/physics/atmosphere.py`)** con sus tests unitarios?
