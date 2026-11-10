"""Pérdida de calor del tubo de retención, según cómo se aísle.

Por qué importa: el sensor legal está a la salida del tubo. Todo lo que el
tubo pierde hay que dárselo antes, en el calentamiento, así que el agua
caliente tiene que estar más caliente. Con la clara eso se come el margen
contra su tope de 62 °C. Y la energía perdida se paga todo el año.

Formas de aislar:
- tubo: cada metro de cañería con su coquilla. La superficie exterior es la
  de un cilindro de ~100 m de largo.
- tubo con las curvas desnudas (solo en horquillas): coquilla en los tramos
  rectos; las curvas de 180° y sus uniones quedan al aire, que es más simple
  de desarmar e inspeccionar. Ese pedazo pierde como un caño pelado.
- bloque (solo en espiral): la espiral entera dentro de una camisa aislante,
  como un termo. La superficie exterior es la de un cilindro de ~1 m × 1,5 m.
La normativa pide que la aislación del tubo de retención se pueda sacar para
inspeccionarlo [cfia2026htst]: coquilla o camisa, tienen que ser desmontables.

Modelo: el tubo es un intercambiador cuyo "otro fluido" es el ambiente, con
capacidad infinita (Cr = 0). Con UA total [W/K] y la capacidad del producto
C = m·cp, la T decae exponencial hacia T_ambiente a lo largo del tubo:

    T_sale − T_amb = (T_entra − T_amb) · exp(−UA/C)

Como el control fija la T de SALIDA del tubo, se despeja cuánto tiene que
entrar más caliente:  ΔT = (T_sale − T_amb) · (exp(UA/C) − 1).

Resistencias (por qué casi todo es el aislante): la película interna del
producto y la pared de acero suman ~1 % de la resistencia total y se
desprecian. La convección exterior h_ext tampoco pesa: con 50 mm de lana, la
resistencia del aislante es 10 a 20 veces la del aire.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class Aislacion:
    forma: str = "bloque"     # bloque | tubo
    espesor: float = 0.05     # m
    k: float = 0.040          # W/(m·K), lana mineral
    h_ext: float = 10.0       # W/(m²·K), convección natural + radiación
    T_amb: float = 20.0       # °C, sala de proceso
    holgura: float = 0.05     # m, luz entre la espiral y la camisa (bloque)
    # Horquillas con coquilla: ¿se aíslan también las curvas de 180°?
    curvas_aisladas: bool = True
    h_desnudo: float = 10.0   # W/(m²·K), caño de acero pelado al aire
    union: float = 0.05       # m sin coquilla a cada lado de una curva (abrazadera)

    def __post_init__(self):
        if self.forma not in ("bloque", "tubo"):
            raise ValueError(f"Forma de aislación desconocida: {self.forma} (bloque | tubo)")
        if self.espesor < 0 or self.k <= 0 or self.h_ext <= 0 or self.h_desnudo <= 0:
            raise ValueError("Aislación: espesor ≥ 0, k > 0, h_ext > 0 y h_desnudo > 0")


@dataclass
class ProductoTermico:
    """Lo que hace falta de cada producto para calcular su pérdida."""
    clave: str
    m: float            # kg/s
    cp: float           # J/(kg·K) a la T del tubo
    T_sale: float       # °C, salida del tubo (la fija el control)
    horas_anio: float


@dataclass
class OpcionAislacion:
    forma: str
    area_ext: float             # m², superficie exterior (revestimiento)
    volumen_aislante: float     # m³ de aislante
    UA: float                   # W/K, conductancia total tubo → ambiente
    dT: dict[str, float]        # °C que pierde cada producto en el tubo
    Q: dict[str, float]         # W que pierde cada producto
    E_MWh_anio: float           # calor perdido al año, que repone el agua caliente


def ua_tubo(L: float, D_ext: float, a: Aislacion) -> tuple[float, float, float]:
    """UA, área exterior y volumen de aislante con coquilla a lo largo del tubo.

    Cilindro: la resistencia de una capa es ln(D_afuera/D_adentro)/(2πk) por
    metro; la del aire exterior, 1/(h·π·D). Se suman en serie.
    """
    D_ais = D_ext + 2 * a.espesor
    R_aislante = math.log(D_ais / D_ext) / (2 * math.pi * a.k) if a.espesor > 0 else 0.0
    R_aire = 1 / (a.h_ext * math.pi * D_ais)
    UA = L / (R_aislante + R_aire)
    area = math.pi * D_ais * L
    volumen = L * math.pi / 4 * (D_ais**2 - D_ext**2)
    return UA, area, volumen


def ua_desnudo(L: float, D_ext: float, h: float) -> tuple[float, float]:
    """UA y área de un caño pelado: solo la convección exterior, h·π·D·L.

    La película interna y la pared se desprecian: pierde un poco más de lo
    real (lado seguro). En laminar la película interna de la yema puede
    sumar ~15 % de resistencia; con la clara y el huevo entero, menos.
    """
    area = math.pi * D_ext * L
    return h * area, area


def ua_bloque(D_espiral: float, alto: float, D_ext: float, a: Aislacion) -> tuple[float, float, float]:
    """UA, área exterior y volumen de aislante con la espiral dentro de una camisa.

    La camisa es un cilindro cerrado (lateral + dos tapas) alrededor de la
    espiral. Se supone que el aire de adentro está a la T del producto: es
    la cota superior de la pérdida (en realidad queda algo más frío).
    Como el espesor es chico frente al diámetro (~5 cm contra ~1 m), la
    pared se trata como plana con el área media.
    """
    D_in = D_espiral + D_ext + 2 * a.holgura
    H_in = alto + 2 * a.holgura
    D_out, H_out = D_in + 2 * a.espesor, H_in + 2 * a.espesor

    def area(D, H):
        return math.pi * D * H + 2 * math.pi * D**2 / 4

    A_in, A_out = area(D_in, H_in), area(D_out, H_out)
    A_media = (A_in + A_out) / 2
    R_aislante = a.espesor / (a.k * A_media)  # K/W, pared plana
    R_aire = 1 / (a.h_ext * A_out)
    return 1 / (R_aislante + R_aire), A_out, A_media * a.espesor


def perdida(UA: float, p: ProductoTermico, T_amb: float) -> float:
    """°C que el producto tiene que entrar más caliente para salir a T_sale."""
    C = p.m * p.cp
    return (p.T_sale - T_amb) * (math.exp(UA / C) - 1)


def evaluar(forma: str, UA: float, area: float, volumen: float,
            productos: list[ProductoTermico], T_amb: float) -> OpcionAislacion:
    dT = {p.clave: perdida(UA, p, T_amb) for p in productos}
    Q = {p.clave: p.m * p.cp * dT[p.clave] for p in productos}
    # cada producto pierde su Q solo durante sus horas de campaña
    E = sum(Q[p.clave] * p.horas_anio for p in productos) / 1e6  # W·h → MWh
    return OpcionAislacion(forma, area, volumen, UA, dT, Q, E)


def opciones(L: float, D_ext: float, D_espiral: float | None, alto: float | None,
             productos: list[ProductoTermico], a: Aislacion,
             L_desnudo: float | None = None) -> dict[str, OpcionAislacion]:
    """Las formas de aislar que tienen sentido para este tubo, para comparar.

    Siempre la coquilla entera. En espiral, además la camisa en bloque. En
    horquillas (L_desnudo = lo que queda al aire si no se aíslan las curvas),
    además la coquilla solo en los tramos rectos.
    """
    res = {"tubo": evaluar("tubo", *ua_tubo(L, D_ext, a), productos, a.T_amb)}
    if D_espiral and alto:
        res["bloque"] = evaluar("bloque", *ua_bloque(D_espiral, alto, D_ext, a), productos, a.T_amb)
    if L_desnudo:
        L_desnudo = min(L_desnudo, L)
        UA_t, area_t, vol_t = ua_tubo(L - L_desnudo, D_ext, a)
        UA_d, area_d = ua_desnudo(L_desnudo, D_ext, a.h_desnudo)
        res["tubo_curvas_desnudas"] = evaluar("tubo_curvas_desnudas", UA_t + UA_d,
                                              area_t + area_d, vol_t, productos, a.T_amb)
    return res


def clave_elegida(a: Aislacion, ops: dict[str, OpcionAislacion]) -> str:
    """Cuál de las opciones calculadas es la de diseño, según la config."""
    if a.forma == "tubo" and not a.curvas_aisladas and "tubo_curvas_desnudas" in ops:
        return "tubo_curvas_desnudas"
    return a.forma
