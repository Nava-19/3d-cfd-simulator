"""
Módulo Exportador de Base de Datos Aerodinámica (AeroDB) (Fase 7.4).
Exporta los resultados estacionarios convergidos (Mach, Altitud, α, CD, CL, Cm, Xcp, T_max)
a formatos estándar de ingeniería de guiado, navegación y control (GNC): JSON y CSV.

Cálculo del Centro de Presiones (X_cp):
    A partir del coeficiente de momento de cabeceo Cm respecto al centro de referencia X_ref
    y el coeficiente de fuerza normal CN ≈ CL * cos(α) + CD * sin(α):
        X_cp = X_ref - (Cm / CN) * L_ref
"""

import csv
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path


@dataclass(frozen=True)
class AeroDBRecord:
    """Registro consolidado de un punto de la envolvente de vuelo para GNC."""
    case_name: str
    mach: float
    altitude_km: float
    alpha_deg: float
    p_inf_pa: float
    t_inf_k: float
    u_inf_m_s: float
    cd: float
    cl: float
    cm_pitch: float
    cd_pressure: float
    cd_friction: float
    x_cp_m: float
    max_wall_temp_k: float
    stagnation_temp_theo_k: float


class AeroDBExporter:
    """Gestiona la exportación y acumulación de resultados aerodinámicos."""

    @staticmethod
    def compute_center_of_pressure(
        cd: float,
        cl: float,
        cm_pitch: float,
        alpha_deg: float,
        x_ref_m: float = 0.5,
        l_ref_m: float = 1.0
    ) -> float:
        """
        Calcula la posición longitudinal del centro de presiones X_cp [m] desde la punta (x=0).
        Si el ángulo de ataque y la fuerza normal son nulos, retorna x_ref_m por simetría axial.
        """
        alpha_rad = math.radians(alpha_deg)
        cn = cl * math.cos(alpha_rad) + cd * math.sin(alpha_rad)
        if abs(cn) < 1e-5:
            return x_ref_m
        return x_ref_m - (cm_pitch / cn) * l_ref_m

    @staticmethod
    def export_json(record: AeroDBRecord, output_json_path: Path | str) -> Path:
        """Guarda el resumen de la simulación en un archivo JSON estructurado."""
        out_path = Path(output_json_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(asdict(record), f, indent=4, ensure_ascii=False)
        return out_path

    @staticmethod
    def append_csv(record: AeroDBRecord, output_csv_path: Path | str) -> Path:
        """Añade una fila a la tabla maestra CSV de la base de datos aerodinámica (AeroDB)."""
        out_path = Path(output_csv_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        row_dict = asdict(record)
        file_exists = out_path.exists() and out_path.stat().st_size > 0

        with open(out_path, "a", encoding="utf-8", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(row_dict.keys()))
            if not file_exists:
                writer.writeheader()
            writer.writerow(row_dict)

        return out_path
