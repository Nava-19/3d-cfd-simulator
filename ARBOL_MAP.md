# Árbol de Descomposición del Problema (Metodología MAP)
## Proyecto: Simulador Aerodinámico y Térmico de Cohetes (OpenFOAM + Python + PyQt6 + PyVista)

```
simulador-cfd-cohetes/
├── configuracion-entorno/
│   ├── crear-archivo-wslconfig.txt                 # limitar RAM a 12GB y habilitar swap en Windows
│   ├── instalar-openfoam-v2312-wsl.txt             # distribución OpenCFD/ESI en Ubuntu WSL2
│   ├── configurar-entorno-virtual-python.txt       # Python 3.10+, pip venv
│   ├── instalar-dependencias-gui-y-calculo.txt     # pyqt6, pyvista, pyqtgraph, jinja2, scipy, numpy
│   └── verificar-aceleracion-gpu-pyvista.txt       # validar renderizado con RTX 4060 en WSL/Windows
├── modulo-atmosfera-isa/
│   ├── implementar-ecuaciones-troposfera-estratosfera.txt # modelo ISA 1976 hasta 30 km
│   ├── calcular-propiedades-termodinamicas.txt     # T_inf, p_inf, rho_inf, mu_inf según altitud
│   ├── calcular-velocidad-sonido-y-flujo.txt       # a = sqrt(gamma*R*T) y U_inf = M * a
│   ├── proyectar-vector-velocidad-incidencia.txt   # componentes Ux, Uy, Uz según ángulo de ataque alpha
│   └── crear-tests-unitarios-isa.txt               # validar salidas frente a tablas estándar ISA
├── preparacion-cad-y-geometria/
│   ├── disenar-modelo-cohete-solidworks.txt        # ojiva, fuselaje, aletas y base
│   ├── aplicar-corte-simetria-longitudinal.txt     # plano y=0 para simular medio cohete (guiñada nula)
│   ├── exportar-parches-stl-independientes.txt     # nose.stl, fuselage.stl, fins.stl, base.stl en ASCII
│   └── verificar-estanqueidad-y-escalado.txt       # comprobar escala en metros y orientación de normales
├── automatizacion-openfoam/
│   ├── plantillas-diccionarios/
│   │   ├── crear-plantilla-blockmeshdict.txt       # dominio exterior con caras farfield y symmetryPlane
│   │   ├── crear-plantilla-snappyhexmeshdict.txt   # refinamiento zonal alto en ojiva/aletas y prism layers
│   │   ├── crear-plantilla-thermophysicalproperties.txt # hePsiThermo, pureMixture, gas perfecto (v2312)
│   │   ├── crear-plantilla-turbulenceproperties.txt # modelo RAS kOmegaSST
│   │   ├── crear-plantilla-fvschemes.txt           # esquemas Kurganov/Tadmor para choques en rhoCentralFoam
│   │   ├── crear-plantilla-fvsolution.txt          # solvers lineales, tolerancias y límites térmicos
│   │   ├── crear-plantilla-controldict.txt         # control CFL < 0.8 y función forceCoeffs por componente
│   │   └── crear-plantilla-condiciones-iniciales-0.txt # campos U, p, T, k, omega con pared adiabática
│   └── orquestador-backend/
│       ├── crear-generador-caso-jinja2.txt         # inyectar Mach, altitud y geometrías en plantillas
│       ├── implementar-hilo-worker-qthread.txt     # subproceso asíncrono para ejecutar comandos en WSL2
│       ├── automatizar-mallado-en-paralelo.txt     # llamadas a surfaceFeatures, blockMesh y snappyHexMesh
│       ├── gestionar-ejecucion-rhocentralfoam.txt  # lanzar solver con redirección a log.rhoCentralFoam
│       └── parsear-log-residuos-y-fuerzas.txt      # regex en streaming para extraer residuos y coeficientes
├── modulo-validacion-teorica/
│   ├── resolver-taylor-maccoll-analitico.txt       # solver ODE en Python para choque cónico exacto
│   ├── crear-caso-benchmark-cono-supersonico.txt   # geometría cónica estándar mallada con snappyHexMesh
│   ├── ejecutar-simulacion-benchmark.txt           # correr cono a Mach 2.0 con rhoCentralFoam
│   ├── extraer-angulo-choque-y-presion.txt         # medir discontinuidad de presión y ángulo en CFD
│   └── generar-tabla-comparativa-error.txt         # calcular % de desviación entre OpenFOAM y teoría
├── interfaz-grafica-pyqt6/
│   ├── arquitectura-base/
│   │   ├── estructurar-proyecto-mvc.txt            # desacoplar modelos de datos, vistas y controladores
│   │   ├── disenar-layout-principal-pyqt6.txt      # barra lateral de control, visor central y panel inferior
│   │   └── conectar-signals-slots-worker.txt       # comunicar worker QThread con elementos de la GUI
│   ├── panel-condiciones-vuelo/
│   │   ├── crear-inputs-mach-altitud-alfa.txt      # spinboxes numéricos con validación de rangos
│   │   ├── crear-selector-archivos-stl.txt         # explorador para asignar parches CAD independientes
│   │   └── mostrar-display-condiciones-isa.txt     # visualizador de T_inf, p_inf y U_inf calculadas
│   ├── panel-convergencia-tiempo-real/
│   │   ├── integrar-widget-pyqtgraph.txt           # canvas de alto rendimiento para refresco continuo
│   │   ├── plotear-curvas-residuos.txt             # actualización en vivo de residuos (p, U, T, k, omega)
│   │   └── plotear-evolucion-fuerzas.txt           # curvas de convergencia de Cd, Cl y Cm
│   └── visor-3d-pyvista/
│       ├── embeber-qvtkwidgewt-en-pyqt6.txt        # integrar viewport VTK acelerado por RTX 4060
│       ├── cargar-geometria-cohete-3d.txt          # renderizado de parches STL con paleta de colores
│       ├── configurar-camara-y-controles.txt       # rotación orbital, paneo, zoom y vistas predefinidas
│       └── superponer-plano-simetria-2d.txt        # plano de corte interactivo para campos de flujo
└── postprocesado-y-visualizacion-avanzada/
    ├── algoritmo-gradiente-densidad-schlieren.txt  # calcular |grad(rho)| para visualizar ondas de choque
    ├── mapear-temperatura-recuperacion-pared.txt   # extraer campo T en superficie 3D (aeroheating)
    ├── generar-lineas-de-corriente-3d.txt          # trazado de partículas de flujo alrededor del cuerpo
    ├── automatizar-capturas-pvpython.txt           # exportación desatendida de renders en alta resolución
    └── exportar-aerodb-csv-json.txt                # exportar matriz de coeficientes aerodinámicos y centro de presiones
```
