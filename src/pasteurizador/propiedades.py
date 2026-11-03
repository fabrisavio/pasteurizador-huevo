"""Propiedades de los productos y del agua.

Térmicas (ρ, cp, k): por composición con Choi y Okos (1986), la misma fuente
que usa el modelo del secador, así los dos equipos hablan el mismo idioma.
La idea: un alimento es una mezcla de agua, proteína, grasa, carbohidratos y
cenizas; cada uno tiene su polinomio en T, y la mezcla se arma como una regla
de mezclas (como el cp de una solución ideal).

Reología (K, n): medida, porque Choi–Okos no la da. Una función por producto.

Todas las temperaturas en °C. Salidas en SI: kg/m³, J/(kg·K), W/(m·K), Pa·sⁿ.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np
from iapws import IAPWS97

# -----------------------------------------------------------------------------
# Choi–Okos (1986): coeficientes de a + b·T + c·T², T en °C
# cp en kJ/(kg·K) (se pasa a J/(kg·K) abajo), ρ en kg/m³, k en W/(m·K)
# -----------------------------------------------------------------------------
_CHOI_OKOS = {
    "agua": {
        "cp": (4.1762, -9.0864e-5, 5.4731e-6),
        "rho": (997.18, 3.1439e-3, -3.7574e-3),
        "k": (0.57109, 1.7625e-3, -6.7036e-6),
    },
    "proteina": {
        "cp": (2.0082, 1.2089e-3, -1.3129e-6),
        "rho": (1329.9, -0.5184, 0.0),
        "k": (0.17881, 1.1958e-3, -2.7178e-6),
    },
    "grasa": {
        "cp": (1.9842, 1.4733e-3, -4.8008e-6),
        "rho": (925.59, -0.41757, 0.0),
        "k": (0.18071, -2.7604e-4, -1.7749e-7),
    },
    "carbohidratos": {
        "cp": (1.5488, 1.9625e-3, -5.9399e-6),
        "rho": (1599.1, -0.31046, 0.0),
        "k": (0.20141, 1.3874e-3, -4.3312e-6),
    },
    "cenizas": {
        "cp": (1.0926, 1.8896e-3, -3.6817e-6),
        "rho": (2423.8, -0.28063, 0.0),
        "k": (0.32962, 1.4011e-3, -2.9069e-6),
    },
}

COMPONENTES = tuple(_CHOI_OKOS)


def _poli(coef: tuple[float, float, float], T: float) -> float:
    a, b, c = coef
    return a + b * T + c * T**2


@dataclass
class Termicas:
    rho: float  # kg/m³
    cp: float   # J/(kg·K)
    k: float    # W/(m·K)


def choi_okos(composicion: dict[str, float], T: float) -> Termicas:
    """ρ, cp y k de un alimento a partir de su composición (fracciones másicas).

    cp: promedio másico.  ρ: volúmenes aditivos.  k: modelo en paralelo con
    fracciones volumétricas, que es el que proponen Choi y Okos.
    """
    faltan = set(composicion) - set(COMPONENTES)
    if faltan:
        raise ValueError(f"Componentes desconocidos para Choi–Okos: {faltan}")
    suma = sum(composicion.values())
    if abs(suma - 1.0) > 0.01:
        raise ValueError(f"La composición suma {suma:.4f}; tiene que sumar 1 (±0,01)")

    x = {i: composicion.get(i, 0.0) / suma for i in COMPONENTES}  # normaliza
    # cp: cada kg de mezcla es x_i kg de cada componente → promedio másico
    cp = sum(x[i] * _poli(_CHOI_OKOS[i]["cp"], T) for i in COMPONENTES) * 1000.0
    # ρ: los volúmenes se suman (1/ρ = Σ x_i/ρ_i), como en una mezcla ideal
    rho_i = {i: _poli(_CHOI_OKOS[i]["rho"], T) for i in COMPONENTES}
    rho = 1.0 / sum(x[i] / rho_i[i] for i in COMPONENTES)
    # k: el calor cruza por los componentes en paralelo, cada uno según el
    # volumen que ocupa (como resistencias en paralelo)
    phi = {i: x[i] / rho_i[i] * rho for i in COMPONENTES}  # fracción volumétrica
    k = sum(phi[i] * _poli(_CHOI_OKOS[i]["k"], T) for i in COMPONENTES)
    return Termicas(rho=rho, cp=cp, k=k)


def gut2003_yema(T: float) -> Termicas:
    """Correlaciones medidas para yema, ecs. (1)–(3) de Gut et al. (2003).

    Válidas entre 0,4 y 60,8 °C. Se usan en la validación y como alternativa
    a Choi–Okos para la yema.
    """
    return Termicas(
        rho=1133.19 - 0.06 * T,
        cp=2629.05 + 2.39 * T,
        k=0.390 + 4e-4 * T,
    )


def termicas(modelo: str, composicion: dict[str, float] | None, T: float) -> Termicas:
    if modelo == "choi_okos":
        if composicion is None:
            raise ValueError("Choi–Okos necesita la composición del producto")
        return choi_okos(composicion, T)
    if modelo == "gut2003_yema":
        return gut2003_yema(T)
    raise ValueError(f"Modelo de propiedades térmicas desconocido: {modelo}")


# -----------------------------------------------------------------------------
# Reología: ley de potencia τ = K·γ̇ⁿ  (n = 1 → newtoniano, K = μ)
# -----------------------------------------------------------------------------
@dataclass
class Reologia:
    K: float                 # Pa·sⁿ
    n: float                 # adimensional
    avisos: list[str] = field(default_factory=list)


def _leon_munguia_huevo_entero(T: float) -> Reologia:
    """Huevo entero, León Munguía (tesis UNAM): K = 0,1·(0,0689 − 0,0007·T), n = 0,955.

    Medido entre 10 y 40 °C. Fuera de ese rango es extrapolación lineal. En el
    tubo (60 °C) la viscosidad real puede ser algo mayor: el huevo entero
    empieza a espesar cerca de 56 °C, aunque no coagula hasta ~73 °C
    (manual AEB). Más viscoso en laminar no cambia el largo (el factor depende
    de n, no de K), pero sí sube el Δp.
    """
    avisos = []
    if not 10.0 <= T <= 40.0:
        avisos.append(
            f"HP: reología extrapolada a {T:.1f} °C (medida entre 10 y 40 °C)"
        )
    K = 0.1 * (0.0689 - 0.0007 * T)
    if K <= 0:
        raise ValueError(f"K de huevo entero no positivo a {T:.1f} °C: fuera de toda validez")
    return Reologia(K=K, n=0.955, avisos=avisos)


# Gut et al. (2003), tablas 4 y 5 (promedios de tres repeticiones)
_GUT_T = np.array([0.4, 8.9, 20.3, 28.4, 37.8, 48.2, 55.2, 60.8])
_GUT_K = np.array([2.26, 1.24, 0.59, 0.35, 0.20, 0.16, 0.14, 0.10])
_GUT_N = np.array([0.840, 0.843, 0.861, 0.854, 0.859, 0.861, 0.878, 0.874])


def _gut2003_yema(T: float) -> Reologia:
    """Yema: interpolación en las tablas de Gut et al. (2003).

    ln K se interpola linealmente (K cae de forma casi exponencial con T).
    Se usan las tablas y no la ec. (4), que sobre 48 °C subestima K en 23–38 %.
    Fuera de 0,4–60,8 °C se extrapola con el último tramo y se avisa.
    Válido para γ̇ entre 70 y 512 s⁻¹.
    """
    avisos = []
    lnK = np.log(_GUT_K)
    if T < _GUT_T[0] or T > _GUT_T[-1]:
        avisos.append(
            f"YP: reología extrapolada a {T:.1f} °C (tablas de Gut entre 0,4 y 60,8 °C)"
        )
        i0, i1 = (0, 1) if T < _GUT_T[0] else (-2, -1)
        pend_lnK = (lnK[i1] - lnK[i0]) / (_GUT_T[i1] - _GUT_T[i0])
        pend_n = (_GUT_N[i1] - _GUT_N[i0]) / (_GUT_T[i1] - _GUT_T[i0])
        K = float(np.exp(lnK[i0] + pend_lnK * (T - _GUT_T[i0])))
        n = float(_GUT_N[i0] + pend_n * (T - _GUT_T[i0]))
    else:
        K = float(np.exp(np.interp(T, _GUT_T, lnK)))
        n = float(np.interp(T, _GUT_T, _GUT_N))
    return Reologia(K=K, n=n, avisos=avisos)


# Scalzo et al. (1970) vía manual AEB, Tabla 3: clara estabilizada pH 7.
# Promedio de Ohio y California (las dos series completas), en mPa·s.
_SCALZO_T = np.array([5.0, 30.0, 40.0, 50.0])
_SCALZO_MU = np.array([(4.8 + 5.5) / 2, (2.7 + 2.8) / 2, (1.9 + 2.3) / 2, (1.7 + 1.9) / 2])
_MU_CP_30C = 9.3e-3  # Pa·s: Thapon y Bourgeois, interpolado a 23 % de sólidos


def _clara_uf_thapon_scalzo(T: float) -> Reologia:
    """Clara concentrada (23 % sólidos), estabilizada: newtoniana.

    μ(T) = 9,3 mPa·s × μ_Scalzo(T)/μ_Scalzo(30 °C). El nivel sale de Thapon
    (concentración) y la forma con T de Scalzo (clara al 12 %): supuesto declarado.
    Sobre 50 °C se extrapola ln μ con el último tramo.

    Es clara estabilizada con ácido láctico (pH ~7) y sal de aluminio, la misma
    condición de Scalzo. Esa estabilización es la que permite pasteurizarla a
    60 °C sin que coagule (Cunningham y Lineweaver, 1965).
    """
    avisos = []
    ln_mu = np.log(_SCALZO_MU)
    if T > _SCALZO_T[-1] or T < _SCALZO_T[0]:
        avisos.append(
            f"CP: forma μ(T) extrapolada a {T:.1f} °C (Scalzo entre 5 y 50 °C)"
        )
        # Fuera de la tabla: se prolonga la recta del tramo más cercano en ln μ
        i0, i1 = (0, 1) if T < _SCALZO_T[0] else (-2, -1)
        pend = (ln_mu[i1] - ln_mu[i0]) / (_SCALZO_T[i1] - _SCALZO_T[i0])
        mu_rel = np.exp(ln_mu[i0] + pend * (T - _SCALZO_T[i0]))
    else:
        mu_rel = np.exp(np.interp(T, _SCALZO_T, ln_mu))
    mu_30 = np.exp(np.interp(30.0, _SCALZO_T, ln_mu))
    return Reologia(K=float(_MU_CP_30C * mu_rel / mu_30), n=1.0, avisos=avisos)


_REOLOGIAS = {
    "leon_munguia_huevo_entero": _leon_munguia_huevo_entero,
    "gut2003_yema": _gut2003_yema,
    "clara_uf_thapon_scalzo": _clara_uf_thapon_scalzo,
}


def reologia(modelo: str, T: float) -> Reologia:
    """K y n del modelo pedido a T. El nombre del modelo viene del YAML."""
    # _REOLOGIAS[modelo] devuelve la función; el (T) la llama.
    # try/except: si el nombre no existe, en vez de un error críptico se
    # explica cuáles hay.
    try:
        return _REOLOGIAS[modelo](T)
    except KeyError:
        raise ValueError(
            f"Modelo de reología desconocido: {modelo}. Opciones: {list(_REOLOGIAS)}"
        ) from None


# -----------------------------------------------------------------------------
# Agua de servicio: IAPWS-IF97
# -----------------------------------------------------------------------------
P_AGUA_MPA = 0.3  # presión de servicio supuesta; en líquido casi no influye


# @lru_cache: la función recuerda los resultados ya calculados. IAPWS es lento
# y los puntos fijos piden el agua a la misma T cientos de veces; con la T
# redondeada a 0,01 °C, la segunda vez la respuesta sale de memoria.
@lru_cache(maxsize=4096)
def _agua_cache(T_redondeada: float) -> Termicas:
    w = IAPWS97(T=T_redondeada + 273.15, P=P_AGUA_MPA)
    return Termicas(rho=w.rho, cp=w.cp * 1000.0, k=w.k)


def agua(T: float) -> Termicas:
    """ρ, cp y k del agua líquida (IAPWS-IF97), T en °C."""
    if not 0.0 < T < 130.0:
        raise ValueError(f"Agua fuera de rango líquido razonable: {T:.1f} °C")
    return _agua_cache(round(T, 2))
