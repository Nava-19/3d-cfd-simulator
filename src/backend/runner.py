"""
Módulo de ejecución asíncrona de OpenFOAM en WSL2 mediante QThread.
Permite ejecutar blockMesh, snappyHexMesh y rhoCentralFoam en segundo plano sin congelar la GUI,
capturando la salida en tiempo real y emitiendo señales Qt con métricas para visualización.
"""

import subprocess
import time
from pathlib import Path
from PyQt6.QtCore import QThread, pyqtSignal

from src.backend.log_parser import OpenFOAMLogParser


def windows_to_wsl_path(win_path: Path | str) -> str:
    """Convierte una ruta de Windows (ej. C:\\Users\\...) en ruta WSL (/mnt/c/Users/...)."""
    p = Path(win_path).resolve()
    drive = p.drive.replace(":", "").lower()
    path_without_drive = p.as_posix().split(":", 1)[-1].lstrip("/")
    return f"/mnt/{drive}/{path_without_drive}"


class OpenFOAMRunner(QThread):
    """Worker asíncrono para orquestar la ejecución de casos de OpenFOAM en WSL2."""

    # Señales Qt para comunicar el progreso a la interfaz gráfica
    status_changed = pyqtSignal(str)              # Mensaje de estado ("Mallando...", etc.)
    log_line_received = pyqtSignal(str)           # Línea de texto sin procesar
    iteration_updated = pyqtSignal(object)        # Objeto IterationMetrics procesado
    force_updated = pyqtSignal(dict)              # Diccionario con Cd, Cl, Cm
    progress_changed = pyqtSignal(float)          # Progreso porcentual (0.0 a 100.0)
    finished_success = pyqtSignal(float)          # Tiempo total de ejecución en segundos
    error_occurred = pyqtSignal(str)              # Mensaje de error si falla

    def __init__(
        self,
        case_dir_windows: Path | str,
        wsl_distro: str = "Ubuntu-22.04",
        run_parallel: bool = False,
        num_cores: int = 6,
        run_blockmesh: bool = False,
        end_time: float | None = None
    ) -> None:
        super().__init__()
        self.case_dir_win = Path(case_dir_windows).resolve()
        self.case_dir_wsl = windows_to_wsl_path(self.case_dir_win)
        self.wsl_distro = wsl_distro
        self.run_parallel = run_parallel
        self.num_cores = num_cores
        self.run_blockmesh = run_blockmesh
        self.end_time = end_time

        self._is_cancelled = False
        self._current_process: subprocess.Popen | None = None
        self.parser = OpenFOAMLogParser()

    def run(self) -> None:
        """Flujo de trabajo completo del hilo."""
        start_time = time.time()
        self._is_cancelled = False

        try:
            self.status_changed.emit("Iniciando pipeline de simulación en WSL2...")

            # 0. Opcionalmente generar malla con blockMesh
            if self.run_blockmesh:
                self.status_changed.emit("Generando malla estructurada (blockMesh)...")
                ok_mesh = self._execute_wsl_command("blockMesh", log_filename="log.blockMesh")
                if not ok_mesh or self._is_cancelled:
                    return

            # 1. Ejecutar solver
            self.status_changed.emit("Ejecutando solver compresible rhoCentralFoam...")
            solver_cmd = "rhoCentralFoam"
            if self.run_parallel:
                solver_cmd = f"mpirun -np {self.num_cores} rhoCentralFoam -parallel"

            success = self._execute_wsl_command(solver_cmd, log_filename="log.rhoCentralFoam")
            if not success or self._is_cancelled:
                if self._is_cancelled:
                    self.status_changed.emit("Simulación abortada por el usuario.")
                return

            elapsed = time.time() - start_time
            self.status_changed.emit(f"Simulación finalizada con éxito en {elapsed:.2f} s")
            self.progress_changed.emit(100.0)
            self.finished_success.emit(elapsed)

        except Exception as e:
            self.error_occurred.emit(str(e))
            self.status_changed.emit(f"Error: {e}")

    def _execute_wsl_command(self, bash_cmd: str, log_filename: str | None = None) -> bool:
        """Lanza un comando bash en WSL2 capturando stdout en streaming."""
        full_bash = f"source /usr/lib/openfoam/openfoam2312/etc/bashrc && cd '{self.case_dir_wsl}' && {bash_cmd}"
        cmd_args = ["wsl.exe", "-d", self.wsl_distro, "-u", "root", "-e", "bash", "-c", full_bash]

        log_file = None
        if log_filename:
            log_path = self.case_dir_win / log_filename
            log_file = open(log_path, "w", encoding="utf-8")

        self._current_process = subprocess.Popen(
            cmd_args,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True
        )

        try:
            assert self._current_process.stdout is not None
            for raw_line in iter(self._current_process.stdout.readline, ""):
                if self._is_cancelled:
                    self._current_process.terminate()
                    break

                line = raw_line.rstrip()
                if log_file:
                    log_file.write(raw_line)
                    log_file.flush()

                self.log_line_received.emit(line)

                # Parsear métricas de la simulación
                metrics = self.parser.parse_line(line)
                if metrics:
                    self.iteration_updated.emit(metrics)
                    if self.end_time and self.end_time > 0:
                        pct = min(max(metrics.time / self.end_time * 100.0, 0.0), 99.9)
                        self.progress_changed.emit(pct)

            self._current_process.stdout.close()
            ret_code = self._current_process.wait()

            if ret_code != 0 and not self._is_cancelled:
                self.error_occurred.emit(f"El comando '{bash_cmd}' falló con código {ret_code}.")
                return False

            return True

        finally:
            if log_file:
                log_file.close()
            self._current_process = None

    def cancel(self) -> None:
        """Detiene la simulación de forma segura."""
        self._is_cancelled = True
        if self._current_process:
            try:
                self._current_process.terminate()
            except Exception:
                pass
