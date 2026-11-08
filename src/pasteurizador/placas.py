"""Placas chevron: geometría del canal, coeficiente h y caída de presión (v1).

Qué reemplaza: en v0 el U de cada sección era un dato (el de la yema de Gut et
al., 2003). Acá U sale de los dos coeficientes de película y de lo que hay en
el medio, como resistencias en serie:

    1/U = 1/h_caliente + 1/h_fría + e/k_placa + R_ensuciamiento

Cada h sale de una correlación de Nusselt para el canal entre dos placas. La
velocidad en el canal la fija el arreglo: cuántos canales en paralelo tiene
cada paso y cuántos pasos en serie. Menos canales por paso = más velocidad =
más h, pero también más Δp. Es el mismo compromiso que el DN del tubo.

Correlaciones (ángulo chevron medido respecto de la dirección del flujo; a 45°
las dos convenciones que circulan en la literatura coinciden):
- Kumar (1984): la de diseño. Constantes por ángulo y rango de Re, tomadas de
  la tabla que reproducen Kakaç et al. (2012). Cubre Re de 0,1 a 10 000, que
  incluye el laminar de la yema.
- Martin (1996), versión del VDI Heat Atlas: contraste semiteórico. Ajustada
  para Re de 200 a 10 000: con los productos viscosos queda fuera de rango.
- Sharifi et al. (1995): placas planas. Solo para revisar el caso de Gut.

Fluidos ley de potencia (yema, huevo entero): viscosidad generalizada de
Delplace y Leuliet (1995), como la usan Gut et al. (2003), con los parámetros
geométricos de placas paralelas (υ = 2, ξ = 12). El canal chevron no es un par
de placas paralelas: es una aproximación declarada.

Caída de presión (Galiński et al., 2025, ec. 2; Kakaç et al., 2012):
    Δp = 4f·(L/Dh)·G²/(2ρ)·pasos  +  1,4·pasos·Gp²/(2ρ)
con f de Fanning. El primer término es el de los canales; el segundo, el de
los puertos (la entrada y la salida de cada paso).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

from . import propiedades as pr


# -----------------------------------------------------------------------------
# Geometría
# -----------------------------------------------------------------------------
@dataclass(frozen=True)
class Placa:
    area: float        # m², área de transferencia de una placa (ya incluye el factor de área)
    ancho: float       # m, ancho del canal
    b: float           # m, separación entre placas (profundidad del corrugado)
    angulo: float      # °, ángulo chevron
    phi: float         # factor de área: desarrollada / proyectada
    espesor: float     # m
    k: float           # W/(m·K), material de la placa
    D_puerto: float    # m, diámetro de los puertos

    @property
    def Dh(self) -> float:
        """Diámetro hidráulico del canal: 2b/φ (Shah y Sekulić, 2003)."""
        return 2 * self.b / self.phi

    @property
    def De(self) -> float:
        """Diámetro equivalente 2b: el que usan Galiński et al. (2025) para el Re de CIP."""
        return 2 * self.b

    @property
    def L(self) -> float:
        """Largo efectivo de puerto a puerto, como si la placa fuera lisa."""
        return self.area / (self.phi * self.ancho)

    @property
    def A_canal(self) -> float:
        """Sección de paso de un canal, m²."""
        return self.b * self.ancho

    @property
    def R_pared(self) -> float:
        return self.espesor / self.k


@dataclass(frozen=True)
class Arreglo:
    """Cómo recorre un fluido su lado de la sección: pasos en serie × canales en paralelo."""
    pasos: int
    canales: int       # canales por paso

    @property
    def canales_lado(self) -> int:
        return self.pasos * self.canales

    @property
    def placas(self) -> int:
        """Placas térmicas con el mismo arreglo de los dos lados.

        Con P pasos de N canales por lado hay 2·P·N canales en total; las
        placas que separan canales son una menos.
        """
        return 2 * self.canales_lado - 1

    def __str__(self) -> str:
        return f"{self.pasos} × {self.canales}"


# -----------------------------------------------------------------------------
# Fluidos
# -----------------------------------------------------------------------------
@dataclass(frozen=True)
class Fluido:
    nombre: str
    termicas: Callable[[float], pr.Termicas]
    reologia: Callable[[float], pr.Reologia]   # agua: K = μ, n = 1


AGUA = Fluido("agua", pr.agua, pr.agua_reologia)


# -----------------------------------------------------------------------------
# Viscosidad generalizada (Delplace y Leuliet, 1995, vía Gut et al., 2003)
# -----------------------------------------------------------------------------
UPSILON_PLACAS, XI_PLACAS = 2.0, 12.0   # placas paralelas


def viscosidad_generalizada(K: float, n: float, v: float, D: float,
                            upsilon: float = UPSILON_PLACAS, xi: float = XI_PLACAS) -> float:
    """μg = K·ξⁿ⁻¹·(v/D)ⁿ⁻¹·((υn + 1)/((υ + 1)n))ⁿ.  Con n = 1 da μ = K.

    Es la viscosidad que hace valer las correlaciones newtonianas para un
    fluido ley de potencia, evaluada con la tasa de corte en la pared.
    """
    return K * xi ** (n - 1) * (v / D) ** (n - 1) * ((upsilon * n + 1) / ((upsilon + 1) * n)) ** n


def tasa_corte_pared(n: float, v: float, D: float,
                     upsilon: float = UPSILON_PLACAS, xi: float = XI_PLACAS) -> float:
    """γ̇ en la pared, 1/s. Sirve para ver si la reología medida vale ahí."""
    return xi * (v / D) * (upsilon * n + 1) / ((upsilon + 1) * n)


# -----------------------------------------------------------------------------
# Correlaciones
# -----------------------------------------------------------------------------
INF = float("inf")
_KUMAR_ANGULOS = (30.0, 45.0, 50.0, 60.0, 65.0)
# Nu = C·Re^y·Pr^(1/3)·(μ/μw)^0,17 ; cada tramo: (Re hasta, C, y)
_KUMAR_NU = {
    30.0: ((10, 0.718, 0.349), (INF, 0.348, 0.663)),
    45.0: ((10, 0.718, 0.349), (100, 0.400, 0.598), (INF, 0.300, 0.663)),
    50.0: ((20, 0.630, 0.333), (300, 0.291, 0.591), (INF, 0.130, 0.732)),
    60.0: ((20, 0.562, 0.326), (400, 0.306, 0.529), (INF, 0.108, 0.703)),
    65.0: ((20, 0.562, 0.326), (500, 0.331, 0.503), (INF, 0.087, 0.718)),
}
# f (Fanning) = Kp / Re^z ; cada tramo: (Re hasta, Kp, z)
_KUMAR_F = {
    30.0: ((10, 50.0, 1.0), (100, 19.40, 0.589), (INF, 2.990, 0.183)),
    45.0: ((15, 47.0, 1.0), (300, 18.29, 0.652), (INF, 1.441, 0.206)),
    50.0: ((20, 34.0, 1.0), (300, 11.25, 0.631), (INF, 0.772, 0.161)),
    60.0: ((40, 24.0, 1.0), (400, 3.24, 0.457), (INF, 0.760, 0.215)),
    65.0: ((50, 24.0, 1.0), (500, 2.80, 0.451), (INF, 0.639, 0.213)),
}


def _tabla_kumar(tabla: dict, angulo: float, Re: float) -> tuple[float, float]:
    # Fila: el primer ángulo de la tabla que no es menor que el pedido (sin interpolar)
    fila = next((tabla[a] for a in _KUMAR_ANGULOS if angulo <= a), tabla[_KUMAR_ANGULOS[-1]])
    _, C, e = next(t for t in fila if Re <= t[0])
    return C, e


def nu_kumar(Re: float, Pr: float, angulo: float, mu_sobre_mu_pared: float = 1.0) -> float:
    C, y = _tabla_kumar(_KUMAR_NU, angulo, Re)
    return C * Re ** y * Pr ** (1 / 3) * mu_sobre_mu_pared ** 0.17


def f_kumar(Re: float, angulo: float) -> float:
    """Factor de fricción de Fanning."""
    Kp, z = _tabla_kumar(_KUMAR_F, angulo, Re)
    return Kp / Re ** z


def xi_martin(Re: float, angulo: float) -> float:
    """Factor de fricción de Darcy de Martin (1996), constantes del VDI Heat Atlas."""
    phi = math.radians(angulo)
    if Re < 2000:
        xi0, xi1 = 64 / Re, 597 / Re + 3.85
    else:
        xi0, xi1 = (1.8 * math.log10(Re) - 1.5) ** -2, 39 * Re ** -0.289
    inv = (math.cos(phi) / math.sqrt(0.18 * math.tan(phi) + 0.36 * math.sin(phi) + xi0 / math.cos(phi))
           + (1 - math.cos(phi)) / math.sqrt(3.8 * xi1))
    return inv ** -2


def nu_martin(Re: float, Pr: float, angulo: float, mu_sobre_mu_pared: float = 1.0) -> float:
    xi = xi_martin(Re, angulo)
    return (0.122 * Pr ** (1 / 3) * mu_sobre_mu_pared ** (1 / 6)
            * (xi * Re ** 2 * math.sin(2 * math.radians(angulo))) ** 0.374)


def nu_sharifi(Re: float, Pr: float, *_) -> float:
    """Placas planas (Sharifi et al., 1995, vía Gut et al., 2003). Sin corrección de pared."""
    return 1.35 * Re ** 0.36 * Pr ** 0.33


CORRELACIONES = {"kumar", "martin", "sharifi"}
RANGO_RE = {"kumar": (0.1, 1e4), "martin": (200.0, 1e4), "sharifi": (0.0, INF)}


def _nu(correlacion: str, Re: float, Pr: float, angulo: float, ratio: float) -> float:
    if correlacion == "kumar":
        return nu_kumar(Re, Pr, angulo, ratio)
    if correlacion == "martin":
        return nu_martin(Re, Pr, angulo, ratio)
    if correlacion == "sharifi":
        return nu_sharifi(Re, Pr)
    raise ValueError(f"Correlación desconocida: {correlacion}. Opciones: {sorted(CORRELACIONES)}")


def _f_fanning(correlacion: str, Re: float, angulo: float) -> float:
    if correlacion == "martin":
        return xi_martin(Re, angulo) / 4   # Darcy → Fanning
    # Kumar también para Sharifi: el caso de Gut no necesita Δp
    return f_kumar(Re, angulo)


# -----------------------------------------------------------------------------
# Un lado de una sección
# -----------------------------------------------------------------------------
@dataclass
class Canal:
    """Lo que pasa en un lado de una sección: hidráulica y coeficiente h."""
    fluido: str
    arreglo: Arreglo
    v: float            # m/s en el canal
    G: float            # kg/(m²·s), por canal
    mu_g: float         # Pa·s, viscosidad generalizada
    Re: float
    Pr: float
    gamma_pared: float  # 1/s
    Nu: float
    h: float            # W/(m²·K)
    f: float            # Fanning
    dp_canales: float   # Pa
    dp_puertos: float   # Pa

    @property
    def dp(self) -> float:
        return self.dp_canales + self.dp_puertos


def canal(fl: Fluido, m: float, arr: Arreglo, placa: Placa, T: float, T_pared: float,
          correlacion: str = "kumar") -> Canal:
    """Hidráulica y h de un fluido de caudal m [kg/s] con el arreglo dado.

    T: temperatura media del fluido en la sección; T_pared: la de la placa
    (para la corrección μ/μ_pared, que en calentamiento sube h y en
    enfriamiento la baja).
    """
    t = fl.termicas(T)
    reo, reo_p = fl.reologia(T), fl.reologia(T_pared)
    G = m / (arr.canales * placa.A_canal)
    v = G / t.rho
    D = placa.Dh
    mu = viscosidad_generalizada(reo.K, reo.n, v, D)
    mu_p = viscosidad_generalizada(reo_p.K, reo_p.n, v, D)
    Re = G * D / mu
    Pr = t.cp * mu / t.k
    Nu = _nu(correlacion, Re, Pr, placa.angulo, mu / mu_p)
    f = _f_fanning(correlacion, Re, placa.angulo)
    dp_c = 4 * f * (placa.L / D) * G ** 2 / (2 * t.rho) * arr.pasos
    Gp = m / (math.pi * placa.D_puerto ** 2 / 4)
    dp_p = 1.4 * arr.pasos * Gp ** 2 / (2 * t.rho)
    # float(): IAPWS devuelve tipos de numpy; así los reportes muestran números comunes
    return Canal(fl.nombre, arr, float(v), float(G), float(mu), float(Re), float(Pr),
                 float(tasa_corte_pared(reo.n, v, D)), float(Nu), float(Nu * t.k / D),
                 float(f), float(dp_c), float(dp_p))


@dataclass(frozen=True)
class Lado:
    fluido: Fluido
    m: float           # kg/s
    arreglo: Arreglo


@dataclass
class CoeficienteU:
    U: float
    caliente: Canal
    fria: Canal
    R_ensuciamiento: float

    @property
    def U_limpio(self) -> float:
        """El U sin ensuciamiento: el del equipo recién limpio."""
        return 1 / (1 / self.U - self.R_ensuciamiento)

    @property
    def fraccion_resistencias(self) -> dict[str, float]:
        """Qué parte de 1/U pone cada resistencia (suma 1)."""
        r = {"película caliente": 1 / self.caliente.h, "película fría": 1 / self.fria.h,
             "ensuciamiento": self.R_ensuciamiento}
        r["pared"] = 1 / self.U - sum(r.values())
        return {k: v * self.U for k, v in r.items()}


def coeficiente_U(placa: Placa, caliente: Lado, fria: Lado, Tc: float, Tf: float,
                  R_ensuciamiento: float, correlacion: str = "kumar") -> CoeficienteU:
    """U con las medias de cada corriente; la pared, a mitad de camino entre ambas."""
    T_pared = (Tc + Tf) / 2
    c = canal(caliente.fluido, caliente.m, caliente.arreglo, placa, Tc, T_pared, correlacion)
    f = canal(fria.fluido, fria.m, fria.arreglo, placa, Tf, T_pared, correlacion)
    U = 1 / (1 / c.h + 1 / f.h + placa.R_pared + R_ensuciamiento)
    return CoeficienteU(U, c, f, R_ensuciamiento)


def funcion_U(placa: Placa, caliente: Lado, fria: Lado, R_ensuciamiento: float,
              correlacion: str = "kumar") -> Callable[[float, float, float, float], float]:
    """U como función de las cuatro T de la sección, para ε-NTU.

    seccion_phe la llama con (Tc_in, Tc_out, Tf_in, Tf_out) en cada vuelta del
    punto fijo, así las propiedades siguen a las temperaturas.
    """
    def U(Tc_in: float, Tc_out: float, Tf_in: float, Tf_out: float) -> float:
        return coeficiente_U(placa, caliente, fria, (Tc_in + Tc_out) / 2, (Tf_in + Tf_out) / 2,
                             R_ensuciamiento, correlacion).U
    return U


def caudal_CIP(placa: Placa, arr: Arreglo, Re_min: float, T: float) -> float:
    """m³/h de solución de limpieza para llegar a Re_min en los canales de un lado.

    Re como lo definen Galiński et al. (2025), ec. 6: Re = V̇·De·ρ/(N·b·w·η),
    con De = 2b y N canales por paso. Se toma la solución como agua a T.
    """
    t, mu = pr.agua(T), pr.agua_reologia(T).K
    V = Re_min * arr.canales * placa.b * placa.ancho * mu / (placa.De * t.rho)
    return float(V * 3600)
