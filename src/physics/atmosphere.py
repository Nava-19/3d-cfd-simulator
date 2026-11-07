"""
Módulo de Atmósfera Estándar Internacional (ISA 1976) y Cinemática de Vuelo.
Calcula propiedades termodinámicas, velocidad del sonido y vectores de velocidad
en función del número de Mach, altitud y ángulo de ataque.
"""

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class FlightCondition:
    """Contenedor inmutable de condiciones aerodinámicas y termodinámicas de vuelo."""
    altitude_m: float
    mach: float
    alpha_deg: float
    temperature_k: float
    pressure_pa: float
    density_kg_m3: float
    dynamic_viscosity_pa_s: float
    kinematic_viscosity_m2_s: float
    speed_of_sound_m_s: float
    velocity_magnitude_m_s: float
    velocity_vector_m_s: tuple[float, float, float]
    dynamic_pressure_pa: float
    stagnation_temperature_k: float


class AtmosphereISA:
    """Implementación de la Atmósfera Estándar Internacional (ISA 1976) hasta 32 km."""

    # Constantes físicas
    G0: float = 9.80665          # Gravedad estándar al nivel del mar [m/s^2]
    R_AIR: float = 287.05287     # Constante específica del aire [J/(kg*K)]
    GAMMA: float = 1.4           # Coeficiente adiabático del aire

    # Nivel del mar (h = 0 m)
    T0: float = 288.15           # [K] (15 °C)
    P0: float = 101325.0         # [Pa]
    RHO0: float = 1.2250         # [kg/m^3]

    # Troposfera (0 a 11 km)
    H_TROPO: float = 11000.0     # [m]
    LAPSE_TROPO: float = -0.0065 # [K/m] (-6.5 K/km)

    # Estratosfera inferior isotérmica (11 a 20 km)
    H_STRATO_1: float = 20000.0  # [m]

    # Estratosfera media con gradiente positivo (20 a 32 km)
    H_STRATO_2: float = 32000.0  # [m]
    LAPSE_STRATO: float = 0.0010 # [K/m] (+1.0 K/km)

    # Parámetros ley de Sutherland para viscosidad dinámica del aire
    MU_REF: float = 1.716e-5     # [Pa*s = kg/(m*s)]
    T_REF: float = 273.15        # [K]
    S_SUTHERLAND: float = 110.4  # [K]

    def __init__(self) -> None:
        # Condiciones en la tropopausa (h = 11000 m)
        self.t_tropo = self.T0 + self.LAPSE_TROPO * self.H_TROPO
        exponent_tropo = -self.G0 / (self.LAPSE_TROPO * self.R_AIR)
        self.p_tropo = self.P0 * (self.t_tropo / self.T0) ** exponent_tropo

        # Condiciones a h = 20000 m (capa isotérmica a t_tropo)
        delta_h_iso = self.H_STRATO_1 - self.H_TROPO
        self.p_strato_1 = self.p_tropo * math.exp(-self.G0 * delta_h_iso / (self.R_AIR * self.t_tropo))

    def get_state_at_altitude(self, altitude_m: float) -> tuple[float, float, float]:
        """
        Calcula (T, p, rho) a una altitud geopotencial dada en metros (0 a 32000 m).
        Retorna:
            temperature_k [K], pressure_pa [Pa], density_kg_m3 [kg/m^3]
        """
        if altitude_m < 0.0 or altitude_m > self.H_STRATO_2:
            raise ValueError(f"Altitud {altitude_m} m fuera del rango válido [0, {self.H_STRATO_2}] m.")

        if altitude_m <= self.H_TROPO:
            # Región Troposfera (0 - 11 km)
            t = self.T0 + self.LAPSE_TROPO * altitude_m
            exponent = -self.G0 / (self.LAPSE_TROPO * self.R_AIR)
            p = self.P0 * (t / self.T0) ** exponent
        elif altitude_m <= self.H_STRATO_1:
            # Región Tropopausa / Estratosfera baja isotérmica (11 - 20 km)
            t = self.t_tropo
            delta_h = altitude_m - self.H_TROPO
            p = self.p_tropo * math.exp(-self.G0 * delta_h / (self.R_AIR * t))
        else:
            # Región Estratosfera media (20 - 32 km)
            delta_h = altitude_m - self.H_STRATO_1
            t = self.t_tropo + self.LAPSE_STRATO * delta_h
            exponent = -self.G0 / (self.LAPSE_STRATO * self.R_AIR)
            p = self.p_strato_1 * (t / self.t_tropo) ** exponent

        # Ecuación de estado del gas ideal
        rho = p / (self.R_AIR * t)
        return t, p, rho

    def get_dynamic_viscosity(self, temperature_k: float) -> float:
        """Calcula la viscosidad dinámica mu [Pa*s] mediante la ley de Sutherland."""
        return (
            self.MU_REF
            * (temperature_k / self.T_REF) ** 1.5
            * (self.T_REF + self.S_SUTHERLAND)
            / (temperature_k + self.S_SUTHERLAND)
        )

    def compute_flight_condition(
        self,
        mach: float,
        altitude_km: float,
        alpha_deg: float = 0.0
    ) -> FlightCondition:
        """
        Genera la condición completa de vuelo requerida por el túnel numérico.
        Args:
            mach: Número de Mach de vuelo (ej. 0.8, 1.5, 3.0)
            altitude_km: Altitud en kilómetros (ej. 0.0, 5.0, 15.0)
            alpha_deg: Ángulo de ataque en grados (ej. 0.0, 2.0, 5.0)
        """
        if mach <= 0.0:
            raise ValueError("El número de Mach debe ser estrictamente positivo.")

        alt_m = altitude_km * 1000.0
        t, p, rho = self.get_state_at_altitude(alt_m)

        mu = self.get_dynamic_viscosity(t)
        nu = mu / rho

        # Velocidad del sonido local
        a = math.sqrt(self.GAMMA * self.R_AIR * t)

        # Módulo de la velocidad
        u_mag = mach * a

        # Vector de velocidad según ángulo de ataque (plano X-Y, guiñada Z = 0)
        alpha_rad = math.radians(alpha_deg)
        ux = u_mag * math.cos(alpha_rad)
        uy = u_mag * math.sin(alpha_rad)
        uz = 0.0

        # Presión dinámica y temperatura de remanso (estancamiento)
        q_inf = 0.5 * rho * (u_mag ** 2)
        t_stag = t * (1.0 + 0.5 * (self.GAMMA - 1.0) * (mach ** 2))

        return FlightCondition(
            altitude_m=alt_m,
            mach=mach,
            alpha_deg=alpha_deg,
            temperature_k=t,
            pressure_pa=p,
            density_kg_m3=rho,
            dynamic_viscosity_pa_s=mu,
            kinematic_viscosity_m2_s=nu,
            speed_of_sound_m_s=a,
            velocity_magnitude_m_s=u_mag,
            velocity_vector_m_s=(ux, uy, uz),
            dynamic_pressure_pa=q_inf,
            stagnation_temperature_k=t_stag
        )
