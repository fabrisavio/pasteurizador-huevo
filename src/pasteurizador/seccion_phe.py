"""Una sección del intercambiador de placas, en contracorriente pura (ε-NTU).

El mismo bloque sirve para regeneración (producto contra producto),
calentamiento (agua contra producto) y enfriamiento (en la validación).

Dos modos:
- diseñar: se conocen las cuatro temperaturas → área necesaria.
- verificar: se conoce el área y las dos entradas → salidas.

v0: U constante por sección, dado desde la config. Las propiedades se evalúan
a la T media de cada corriente, que es lo que recomiendan Gut et al. (2003).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Callable

from .propiedades import Termicas

TOL_T = 1e-6       # °C, convergencia de los puntos fijos
MAX_ITER = 200


@dataclass
class Corriente:
    nombre: str
    m: float                              # kg/s
    props: Callable[[float], Termicas]    # T [°C] → propiedades

    def C(self, T1: float, T2: float) -> float:
        """Capacidad calorífica m·cp en W/K, con cp a la T media."""
        return self.m * self.props((T1 + T2) / 2).cp


@dataclass
class Seccion:
    nombre: str
    A: float      # m²
    U: float      # W/(m²·K)
    Q: float      # W
    NTU: float
    eps: float
    Cr: float
    Tc_in: float
    Tc_out: float
    Tf_in: float
    Tf_out: float
    m_c: float
    m_f: float

    @property
    def approach_extremo_caliente(self) -> float:
        return self.Tc_in - self.Tf_out

    @property
    def approach_extremo_frio(self) -> float:
        return self.Tc_out - self.Tf_in


# -----------------------------------------------------------------------------
# ε-NTU en contracorriente
# -----------------------------------------------------------------------------
def eps_contracorriente(NTU: float, Cr: float) -> float:
    if NTU < 0 or not 0 <= Cr <= 1:
        raise ValueError(f"NTU={NTU}, Cr={Cr} fuera de dominio")
    if abs(1 - Cr) < 1e-9:
        return NTU / (1 + NTU)
    e = math.exp(-NTU * (1 - Cr))
    return (1 - e) / (1 - Cr * e)


def ntu_contracorriente(eps: float, Cr: float) -> float:
    """Inversa exacta de eps_contracorriente."""
    if not 0 <= eps < 1 or not 0 <= Cr <= 1:
        raise ValueError(f"ε={eps}, Cr={Cr} fuera de dominio")
    if abs(1 - Cr) < 1e-9:
        return eps / (1 - eps)
    if eps * Cr >= 1:
        raise ValueError(f"ε·Cr = {eps*Cr:.3f} ≥ 1: no hay área que alcance")
    return math.log((1 - eps * Cr) / (1 - eps)) / (1 - Cr)


# -----------------------------------------------------------------------------
# Balances
# -----------------------------------------------------------------------------
def temperatura_por_balance(c: Corriente, T_conocida: float, Q: float) -> float:
    """T de la otra punta de una corriente que intercambia Q [W].

    Q > 0: la corriente se calienta.  Q < 0: se enfría.  cp a la T media,
    así que se itera.
    """
    T = T_conocida + Q / c.C(T_conocida, T_conocida)
    for _ in range(MAX_ITER):
        T_nueva = T_conocida + Q / c.C(T_conocida, T)
        if abs(T_nueva - T) < TOL_T:
            return T_nueva
        T = T_nueva
    raise RuntimeError(f"Balance de {c.nombre} no converge")


# -----------------------------------------------------------------------------
# Diseño y verificación
# -----------------------------------------------------------------------------
def disenar_seccion(nombre: str, caliente: Corriente, fria: Corriente,
                    Tc_in: float, Tc_out: float, Tf_in: float, Tf_out: float,
                    U: float) -> Seccion:
    """Área necesaria para cumplir las cuatro temperaturas dadas."""
    if Tc_in <= Tf_out or Tc_out <= Tf_in:
        raise ValueError(
            f"{nombre}: cruce de temperaturas (caliente {Tc_in:.2f}→{Tc_out:.2f}, "
            f"fría {Tf_in:.2f}→{Tf_out:.2f}); no hay contracorriente que lo haga"
        )
    C_c = caliente.C(Tc_in, Tc_out)
    C_f = fria.C(Tf_in, Tf_out)
    Q_c = C_c * (Tc_in - Tc_out)
    Q_f = C_f * (Tf_out - Tf_in)
    if abs(Q_c - Q_f) > 5e-3 * max(Q_c, Q_f):
        raise ValueError(f"{nombre}: balance no cierra (Q_c={Q_c:.1f} W, Q_f={Q_f:.1f} W)")
    Q = (Q_c + Q_f) / 2
    C_min, C_max = min(C_c, C_f), max(C_c, C_f)
    Cr = C_min / C_max
    eps = Q / (C_min * (Tc_in - Tf_in))
    NTU = ntu_contracorriente(eps, Cr)
    A = NTU * C_min / U
    return Seccion(nombre, A, U, Q, NTU, eps, Cr, Tc_in, Tc_out, Tf_in, Tf_out,
                   caliente.m, fria.m)


def verificar_seccion(nombre: str, caliente: Corriente, fria: Corriente,
                      Tc_in: float, Tf_in: float, A: float, U: float) -> Seccion:
    """Temperaturas de salida de una sección de área A conocida."""
    if Tc_in <= Tf_in:
        raise ValueError(f"{nombre}: la corriente caliente entra más fría que la fría")
    Tc_out, Tf_out = Tc_in, Tf_in  # arranque: sin intercambio
    for _ in range(MAX_ITER):
        C_c = caliente.C(Tc_in, Tc_out)
        C_f = fria.C(Tf_in, Tf_out)
        C_min, C_max = min(C_c, C_f), max(C_c, C_f)
        Cr = C_min / C_max
        NTU = U * A / C_min
        eps = eps_contracorriente(NTU, Cr)
        Q = eps * C_min * (Tc_in - Tf_in)
        Tc_nueva = Tc_in - Q / C_c
        Tf_nueva = Tf_in + Q / C_f
        if abs(Tc_nueva - Tc_out) < TOL_T and abs(Tf_nueva - Tf_out) < TOL_T:
            return Seccion(nombre, A, U, Q, NTU, eps, Cr, Tc_in, Tc_nueva,
                           Tf_in, Tf_nueva, caliente.m, fria.m)
        Tc_out, Tf_out = Tc_nueva, Tf_nueva
    raise RuntimeError(f"{nombre}: la verificación no converge")
