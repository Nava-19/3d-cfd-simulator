"""
Solver analítico exacto de Taylor-Maccoll para flujo supersónico sobre conos simétricos.
Resuelve la EDO no lineal de flujo cónico axisimétrico mediante integración numérica (SciPy)
para obtener el ángulo de onda de choque oblicua y la relación de presiones en la superficie.
"""

import math
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


@dataclass(frozen=True)
class TaylorMaccollResult:
    """Resultados de la solución analítica de Taylor-Maccoll."""
    mach_inf: float
    cone_angle_deg: float
    shock_angle_deg: float
    pressure_ratio_surface: float   # p_cone / p_inf
    mach_surface: float
    theta_array_deg: np.ndarray
    pressure_ratio_array: np.ndarray


class TaylorMaccollSolver:
    """
    Resuelve el flujo cónico supersónico no viscoso sobre un cono sin guiñada (alfa = 0).
    Ecuaciones estándar de Taylor-Maccoll normalizadas con la velocidad máxima V_max.
    """

    def __init__(self, gamma: float = 1.4) -> None:
        self.gamma = gamma
        self.gm1_over_2 = (gamma - 1.0) / 2.0

    def _taylor_maccoll_ode(self, theta: float, y: list[float]) -> list[float]:
        """
        Sistema de EDOs:
        y[0] = v_r (velocidad radial / V_max)
        y[1] = v_theta (velocidad azimutal / V_max)
        """
        vr, vtheta = y[0], y[1]
        a2_over_vmax2 = self.gm1_over_2 * (1.0 - vr**2 - vtheta**2)

        # Evitar singularidades
        if abs(math.sin(theta)) < 1e-9:
            cot_theta = 0.0
        else:
            cot_theta = math.cos(theta) / math.sin(theta)

        denominator = a2_over_vmax2 - vtheta**2
        if abs(denominator) < 1e-12:
            denominator = 1e-12 if denominator >= 0 else -1e-12

        dvr = vtheta
        dvtheta = (vtheta**2 * vr - a2_over_vmax2 * (2.0 * vr + vtheta * cot_theta)) / denominator

        return [dvr, dvtheta]

    def _integrate_from_shock(
        self,
        mach_inf: float,
        shock_angle_rad: float
    ) -> tuple[float, np.ndarray, np.ndarray]:
        """
        Integra hacia atrás desde la onda de choque (theta_s) hacia la superficie del cono
        hasta que v_theta(theta_c) == 0.
        """
        m_inf = mach_inf
        beta = shock_angle_rad
        g = self.gamma

        # Componente normal del Mach incidente
        mn1 = m_inf * math.sin(beta)
        if mn1 <= 1.0:
            raise ValueError(f"Mach normal mn1={mn1:.3f} debe ser > 1 para que exista onda de choque.")

        # Salto de velocidad y ángulo de deflexión delta_s tras el choque oblicuo
        tan_delta = 2.0 * (1.0 / math.tan(beta)) * (mn1**2 - 1.0) / (m_inf**2 * (g + math.cos(2.0 * beta)) + 2.0)
        delta_s = math.atan(tan_delta)

        # Mach normal aguas abajo del choque mn2
        mn2_sq = (1.0 + self.gm1_over_2 * mn1**2) / (g * mn1**2 - self.gm1_over_2)
        m2 = math.sqrt(mn2_sq) / math.sin(beta - delta_s)

        # Relación de velocidad V2 / V_max
        v2_ratio = math.sqrt(1.0 / (1.0 / (self.gm1_over_2 * m2**2) + 1.0))

        # En coordenadas esféricas en la onda de choque theta = beta:
        vr_shock = v2_ratio * math.cos(beta - delta_s)
        vtheta_shock = -v2_ratio * math.sin(beta - delta_s)

        # Integrar desde beta hacia cero buscando el cruce v_theta = 0
        def zero_vtheta_event(theta, y):
            return y[1]
        zero_vtheta_event.terminal = True
        zero_vtheta_event.direction = 1  # Cruce de negativo a cero

        sol = solve_ivp(
            self._taylor_maccoll_ode,
            t_span=(beta, 0.01),
            y0=[vr_shock, vtheta_shock],
            events=zero_vtheta_event,
            max_step=0.002,
            rtol=1e-7,
            atol=1e-9
        )

        if len(sol.t_events[0]) == 0:
            cone_angle = sol.t[-1]
        else:
            cone_angle = sol.t_events[0][0]

        return cone_angle, sol.t, sol.y

    def solve(
        self,
        mach_inf: float,
        cone_angle_deg: float
    ) -> TaylorMaccollResult:
        """
        Calcula la solución exacta de onda de choque débil para un semiángulo de cono y Mach dados.
        """
        if mach_inf <= 1.0:
            raise ValueError("El número de Mach debe ser estrictamente supersónico (M > 1.0).")

        theta_c_target = math.radians(cone_angle_deg)

        # Límite inferior físico: onda de Mach asintótica
        mu_inf = math.asin(1.0 / mach_inf)
        b_curr = max(mu_inf + math.radians(0.2), theta_c_target + math.radians(0.1))
        step = math.radians(0.5)
        max_search_angle = math.radians(75.0)

        # Buscar el intervalo que acota la solución de choque débil
        b_prev = b_curr
        found_bracket = False

        while b_curr < max_search_angle:
            try:
                theta_c_calc, _, _ = self._integrate_from_shock(mach_inf, b_curr)
                if theta_c_calc >= theta_c_target:
                    found_bracket = True
                    break
            except Exception:
                pass
            b_prev = b_curr
            b_curr += step

        if not found_bracket:
            raise ValueError(
                f"No se encontró onda de choque unida para cono de {cone_angle_deg}° a Mach {mach_inf}. "
                "Es probable que el ángulo supere el límite de choque desprendido (detached shock)."
            )

        # Encontrar la raíz exacta por el método de Brent en la rama de choque débil
        def objective(b_val: float) -> float:
            c_calc, _, _ = self._integrate_from_shock(mach_inf, b_val)
            return c_calc - theta_c_target

        beta_opt = brentq(objective, b_prev, b_curr, xtol=1e-6)

        # Integración final con el ángulo de choque débil convergido
        cone_angle_rad, t_hist, y_hist = self._integrate_from_shock(mach_inf, beta_opt)

        # Cálculo de presiones
        mn1 = mach_inf * math.sin(beta_opt)
        p2_over_p1 = 1.0 + (2.0 * self.gamma / (self.gamma + 1.0)) * (mn1**2 - 1.0)
        rho2_over_rho1 = ((self.gamma + 1.0) * mn1**2) / (2.0 + (self.gamma - 1.0) * mn1**2)
        p02_over_p01 = p2_over_p1 * (rho2_over_rho1 ** (-self.gamma)) # Pérdida de presión total

        # Presión isentrópica desde el choque hasta el cono
        vr_c = y_hist[0][-1]
        vtheta_c = y_hist[1][-1]
        v_surf_sq = vr_c**2 + vtheta_c**2
        m_surf = math.sqrt(v_surf_sq / (self.gm1_over_2 * (1.0 - v_surf_sq)))

        p01_over_p1 = (1.0 + self.gm1_over_2 * mach_inf**2) ** (self.gamma / (self.gamma - 1.0))
        p02_over_p_surf = (1.0 + self.gm1_over_2 * m_surf**2) ** (self.gamma / (self.gamma - 1.0))

        # Relación de presión en pared respecto a infinito
        p_surf_over_p_inf = (p01_over_p1 * p02_over_p01) / p02_over_p_surf

        return TaylorMaccollResult(
            mach_inf=mach_inf,
            cone_angle_deg=cone_angle_deg,
            shock_angle_deg=math.degrees(beta_opt),
            pressure_ratio_surface=p_surf_over_p_inf,
            mach_surface=m_surf,
            theta_array_deg=np.degrees(t_hist),
            pressure_ratio_array=np.zeros_like(t_hist) # placeholder
        )
