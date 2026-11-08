"""
Punto de Entrada Principal del Simulador CFD 3D de Cohetes.
Configura el entorno gráfico para alto DPI y aceleración en GPU dedicada (NVIDIA RTX 4060)
antes de iniciar el bucle de eventos de PyQt6.
"""

import os
import sys

# Forzar vinculación de pyvistaqt con PyQt6 y preferencia de GPU de alto rendimiento
os.environ["QT_API"] = "pyqt6"
os.environ["SHIM_MCCOMPAT"] = "0x800000001"

from PyQt6.QtWidgets import QApplication
from src.gui.main_window import CFDMainWindow


def main() -> int:
    app = QApplication(sys.argv)
    app.setApplicationName("RocketCFD3D")
    window = CFDMainWindow(off_screen=False)
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
