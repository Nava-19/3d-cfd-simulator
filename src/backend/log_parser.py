"""
Parser de logs en streaming para rhoCentralFoam y forceCoeffs.
Extrae residuos, número de Courant y coeficientes aerodinámicos (Cd, Cl, Cm)
mediante expresiones regulares optimizadas para alimentar las gráficas en tiempo real.
"""

import re
from dataclasses import dataclass, field


@dataclass
class IterationMetrics:
    """Métricas numéricas de una iteración / paso de tiempo en OpenFOAM."""
    time: float = 0.0
    delta_t: float = 0.0
    courant_max: float = 0.0
    residuals: dict[str, float] = field(default_factory=dict)
    cd: float | None = None
    cl: float | None = None
    cm_pitch: float | None = None
    cd_pressure: float | None = None
    cd_friction: float | None = None


class OpenFOAMLogParser:
    """Parsea líneas de log de OpenFOAM en tiempo real."""

    # Regex para paso de tiempo y Courant (soporta formato estándar y formato rhoCentralFoam)
    RE_TIME = re.compile(r"^Time\s*=\s*([0-9\.eE\+\-]+)")
    RE_COURANT = re.compile(r"Courant Number mean:\s*([0-9\.eE\+\-]+)\s*max:\s*([0-9\.eE\+\-]+)")
    RE_COURANT_CENTRAL = re.compile(r"Mean and max Courant Numbers\s*=\s*([0-9\.eE\+\-]+)\s+([0-9\.eE\+\-]+)")
    RE_DELTAT = re.compile(r"^deltaT\s*=\s*([0-9\.eE\+\-]+)")

    # Regex para residuos: Solving for <field>, Initial residual = <val>, ...
    RE_RESIDUAL = re.compile(r"Solving for\s+([a-zA-Z0-9_]+),\s+Initial residual\s*=\s*([0-9\.eE\+\-]+)")

    # Regex para salida en vivo de forceCoeffs en stdout (v2312):
    # Cd:   <Total>   <Pressure>   <Viscous>   <Internal>
    RE_STDOUT_CD = re.compile(
        r"^Cd:\s+([0-9\.eE\+\-]+)\s+([0-9\.eE\+\-]+)\s+([0-9\.eE\+\-]+)"
    )
    RE_STDOUT_CL = re.compile(
        r"^Cl:\s+([0-9\.eE\+\-]+)\s+([0-9\.eE\+\-]+)\s+([0-9\.eE\+\-]+)"
    )
    RE_STDOUT_CM = re.compile(
        r"^CmPitch:\s+([0-9\.eE\+\-]+)"
    )

    # Regex para forceCoeffs en archivo coefficient.dat
    RE_FORCE_DAT = re.compile(
        r"^\s*([0-9\.eE\+\-]+)\s+([0-9\.eE\+\-]+)\s+([0-9\.eE\+\-]+)\s+([0-9\.eE\+\-]+)\s+"
        r"([0-9\.eE\+\-]+)\s+([0-9\.eE\+\-]+)\s+([0-9\.eE\+\-]+)\s+([0-9\.eE\+\-]+)\s+([0-9\.eE\+\-]+)"
    )

    def __init__(self) -> None:
        self.current_metrics = IterationMetrics()
        self.history: list[IterationMetrics] = []
        self._pending_courant_max: float = 0.0
        self._pending_delta_t: float = 0.0

    def parse_line(self, line: str) -> IterationMetrics | None:
        """
        Procesa una línea de salida de rhoCentralFoam.
        Retorna la métrica completa de la iteración anterior cuando cambia el Time.
        """
        stripped = line.strip()

        # deltaT (en rhoCentralFoam se imprime justo antes de Time =)
        m_dt = self.RE_DELTAT.match(stripped)
        if m_dt:
            self._pending_delta_t = float(m_dt.group(1))
            self.current_metrics.delta_t = self._pending_delta_t
            return None

        # Courant Max (formato rhoCentralFoam antes de Time =)
        m_co_central = self.RE_COURANT_CENTRAL.search(stripped)
        if m_co_central:
            self._pending_courant_max = float(m_co_central.group(2))
            self.current_metrics.courant_max = self._pending_courant_max
            return None

        # Nuevo Time step detectado
        m_time = self.RE_TIME.match(stripped)
        if m_time:
            new_time = float(m_time.group(1))
            completed_metric = None
            if self.current_metrics.residuals or self.current_metrics.cd is not None:
                completed_metric = self.current_metrics
                self.history.append(completed_metric)

            self.current_metrics = IterationMetrics(
                time=new_time,
                delta_t=self._pending_delta_t,
                courant_max=self._pending_courant_max
            )
            return completed_metric

        # Courant Max (formato tradicional dentro del paso de tiempo)
        m_co = self.RE_COURANT.search(stripped)
        if m_co:
            self.current_metrics.courant_max = float(m_co.group(2))
            return None

        # Residuos
        m_res = self.RE_RESIDUAL.search(stripped)
        if m_res:
            field_name = m_res.group(1)
            residual_val = float(m_res.group(2))
            self.current_metrics.residuals[field_name] = residual_val
            return None

        # Coeficientes aerodinámicos impresos en vivo en stdout
        m_cd = self.RE_STDOUT_CD.match(stripped)
        if m_cd:
            self.current_metrics.cd = float(m_cd.group(1))
            self.current_metrics.cd_pressure = float(m_cd.group(2))
            self.current_metrics.cd_friction = float(m_cd.group(3))
            return None

        m_cl = self.RE_STDOUT_CL.match(stripped)
        if m_cl:
            self.current_metrics.cl = float(m_cl.group(1))
            return None

        m_cm = self.RE_STDOUT_CM.match(stripped)
        if m_cm:
            self.current_metrics.cm_pitch = float(m_cm.group(1))
            return None

        return None

    def parse_force_coeffs_line(self, line: str) -> dict[str, float] | None:
        """Parsea una línea del archivo coefficient.dat generado por forceCoeffs."""
        stripped = line.strip()
        if stripped.startswith("#") or not stripped:
            return None

        m = self.RE_FORCE_DAT.match(stripped)
        if m:
            t = float(m.group(1))
            cd = float(m.group(2))
            cl = float(m.group(4))
            cm_pitch = float(m.group(6))
            cdf = float(m.group(8))
            cdp = float(m.group(9))

            data = {
                "time": t,
                "cd": cd,
                "cl": cl,
                "cm_pitch": cm_pitch,
                "cd_friction": cdf,
                "cd_pressure": cdp
            }

            # Sincronizar con la métrica actual
            self.current_metrics.cd = cd
            self.current_metrics.cl = cl
            self.current_metrics.cm_pitch = cm_pitch
            self.current_metrics.cd_friction = cdf
            self.current_metrics.cd_pressure = cdp
            return data

        return None
