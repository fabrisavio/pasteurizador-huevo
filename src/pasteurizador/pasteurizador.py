"""El pasteurizador completo, en dos pasadas.

1. Diseño por producto: qué área necesitaría cada uno por separado.
2. Equipo único: en cada sección manda el producto que más área pide.
3. Verificación: cómo funciona cada producto en ese equipo real.

Esquema por producto (sin enfriamiento: sale al secador):

    crudo ─► [regeneración, lado frío] ─► [calentamiento] ─► tubo de retención ─┐
                                                                                │
    al secador ◄── [regeneración, lado caliente] ◄─────────────────────────────┘
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from scipy.optimize import brentq

from . import propiedades as pr
from .config import Config, Producto
from .holding import ProductoTubo, Tubo, elegir_tubo
from .seccion_phe import (Corriente, Seccion, disenar_seccion,
                          temperatura_por_balance, verificar_seccion)

MARGEN_REGEN_C = 0.5  # si el producto entra a menos de esto de T_salida, no se regenera


def _corriente_producto(p: Producto) -> Corriente:
    return Corriente(p.clave, p.m, p.termicas)


# -----------------------------------------------------------------------------
# 1. Diseño de un producto
# -----------------------------------------------------------------------------
@dataclass
class DisenoProducto:
    clave: str
    regen: Seccion | None
    calent: Seccion
    m_agua: float          # kg/s de agua caliente
    T_regen_salida: float  # °C, producto crudo a la salida de la regeneración


def disenar_producto(p: Producto, cfg: Config) -> DisenoProducto:
    prod = _corriente_producto(p)
    T_calent_salida = p.T_holding + cfg.perdida_holding  # el sensor legal está a la salida del tubo

    # Regeneración: el pasteurizado (caliente) baja de T_holding a T_salida
    # y le entrega ese calor al crudo (frío), que sube por balance.
    if p.T_salida > p.T_entrada + MARGEN_REGEN_C:
        Q_regen = prod.C(p.T_holding, p.T_salida) * (p.T_holding - p.T_salida)
        T_r = temperatura_por_balance(prod, p.T_entrada, Q_regen)
        regen = disenar_seccion("regeneración", prod, prod, p.T_holding, p.T_salida,
                                p.T_entrada, T_r, cfg.U_regen)
    else:
        regen, T_r = None, p.T_entrada

    # Calentamiento: el agua entra al tope permitido (protege la pared) y baja DT
    Q_cal = prod.C(T_r, T_calent_salida) * (T_calent_salida - T_r)
    T_w_in, T_w_out = p.T_max_agua, p.T_max_agua - cfg.DT_agua
    m_w = Q_cal / (pr.agua((T_w_in + T_w_out) / 2).cp * cfg.DT_agua)
    agua = Corriente("agua caliente", m_w, pr.agua)
    # cerrar el balance del agua con su cp exacto
    T_w_out = temperatura_por_balance(agua, T_w_in, -Q_cal)
    calent = disenar_seccion("calentamiento", agua, prod, T_w_in, T_w_out,
                             T_r, T_calent_salida, cfg.U_calent)
    return DisenoProducto(p.clave, regen, calent, m_w, T_r)


# -----------------------------------------------------------------------------
# 2. Equipo único
# -----------------------------------------------------------------------------
@dataclass
class Equipo:
    placas_regen: int
    placas_calent: int
    A_regen: float         # m² instalados
    A_calent: float
    gobierna_regen: str
    gobierna_calent: str
    m_agua: float          # kg/s, caudal de la bomba de agua caliente
    tubo: Tubo
    disenos: dict[str, DisenoProducto]


def _placas(A: float, A_placa: float) -> int:
    """Placas térmicas necesarias (cada una aporta A_placa de área)."""
    return max(1, math.ceil(A / A_placa - 1e-9))


def productos_tubo(cfg: Config) -> list[ProductoTubo]:
    lista = []
    for p in cfg.productos.values():
        T_media = p.T_holding + cfg.perdida_holding / 2
        reo = p.reologia(T_media)
        lista.append(ProductoTubo(
            clave=p.clave, m=p.m, rho=p.termicas(T_media).rho, K=reo.K, n=reo.n,
            T_holding=p.T_holding, t_holding=p.t_holding, cinetica=p.cinetica,
        ))
    return lista


def disenar_equipo(cfg: Config) -> Equipo:
    disenos = {k: disenar_producto(p, cfg) for k, p in cfg.productos.items()}

    A_reg = {k: (d.regen.A if d.regen else 0.0) for k, d in disenos.items()}
    A_cal = {k: d.calent.A for k, d in disenos.items()}
    gob_reg = max(A_reg, key=A_reg.get)
    gob_cal = max(A_cal, key=A_cal.get)
    n_reg = _placas(A_reg[gob_reg], cfg.area_placa)
    n_cal = _placas(A_cal[gob_cal], cfg.area_placa)

    tubo = elegir_tubo(productos_tubo(cfg), cfg.DN_diam, cfg.criterios,
                       cfg.dp_max_holding, cfg.v_CIP, cfg.geometria_tubo, cfg.DN_ext,
                       cfg.factores_sensibilidad)
    return Equipo(
        placas_regen=n_reg, placas_calent=n_cal,
        A_regen=n_reg * cfg.area_placa, A_calent=n_cal * cfg.area_placa,
        gobierna_regen=gob_reg, gobierna_calent=gob_cal,
        m_agua=max(d.m_agua for d in disenos.values()),
        tubo=tubo, disenos=disenos,
    )


# -----------------------------------------------------------------------------
# 3. Verificación de cada producto en el equipo real
# -----------------------------------------------------------------------------
@dataclass
class Verificacion:
    clave: str
    T_entrada: float
    regen: Seccion
    calent: Seccion | None
    T_agua_necesaria: float  # °C a la entrada del calentamiento
    T_salida: float          # °C, producto hacia el secador
    Q_calent_kW: float
    Q_regen_kW: float
    E_calor_MWh_anio: float
    K_salida: float          # reología a la T de salida (para el atomizador)
    n_salida: float
    alertas: list[str] = field(default_factory=list)  # no cumple algo: hay que actuar
    avisos: list[str] = field(default_factory=list)   # informativo

    @property
    def ok(self) -> bool:
        return self.calent is not None


def verificar_producto(p: Producto, eq: Equipo, cfg: Config,
                       T_entrada: float | None = None) -> Verificacion:
    T_in = p.T_entrada if T_entrada is None else T_entrada
    prod = _corriente_producto(p)
    T_calent_salida = p.T_holding + cfg.perdida_holding
    alertas: list[str] = []
    avisos: list[str] = []

    # Regeneración: la T de holding la fija el control, así que el lazo se abre:
    # entra caliente a T_holding y frío a T_entrada, con el área instalada.
    regen = verificar_seccion("regeneración", prod, prod, p.T_holding, T_in,
                              eq.A_regen, cfg.U_regen)
    T_r, T_salida = regen.Tf_out, regen.Tc_out

    # Calentamiento: con el caudal de agua fijo, ¿a qué T tiene que entrar el agua
    # para que el producto llegue a T_holding + pérdida?
    agua = Corriente("agua caliente", eq.m_agua, pr.agua)

    def falta(T_w: float) -> float:
        s = verificar_seccion("calentamiento", agua, prod, T_w, T_r, eq.A_calent, cfg.U_calent)
        return s.Tf_out - T_calent_salida

    T_w_max = 99.0
    if falta(T_w_max) < 0:
        alertas.append(f"{p.clave} (entrada {T_in:.1f} °C): ni con agua a {T_w_max:.0f} °C llega a {T_calent_salida:.1f} °C")
        calent, T_w = None, float("nan")
    else:
        T_w = brentq(falta, T_calent_salida + 1e-6, T_w_max, xtol=1e-6)
        calent = verificar_seccion("calentamiento", agua, prod, T_w, T_r,
                                   eq.A_calent, cfg.U_calent)
        if T_w > p.T_max_agua + 1e-6:
            alertas.append(
                f"{p.clave} (entrada {T_in:.1f} °C): necesita agua a {T_w:.1f} °C, sobre el tope de {p.T_max_agua:.0f} °C"
            )

    if T_salida > p.T_salida + 0.05:
        alertas.append(f"{p.clave} (entrada {T_in:.1f} °C): sale a {T_salida:.1f} °C, sobre los {p.T_salida:.0f} °C pedidos")

    reo = p.reologia(T_salida)
    K_ref = p.reologia(p.T_salida).K
    if reo.K > 1.2 * K_ref:
        avisos.append(
            f"{p.clave}: a {T_salida:.1f} °C su K es {reo.K/K_ref:.1f} veces el de "
            f"{p.T_salida:.0f} °C; revisar la atomización"
        )
    Q_cal = calent.Q / 1000 if calent else float("nan")
    return Verificacion(
        clave=p.clave, T_entrada=T_in, regen=regen, calent=calent,
        T_agua_necesaria=T_w, T_salida=T_salida,
        Q_calent_kW=Q_cal, Q_regen_kW=regen.Q / 1000,
        E_calor_MWh_anio=Q_cal * p.horas_anio / 1000,
        K_salida=reo.K, n_salida=reo.n, alertas=alertas, avisos=avisos,
    )


# -----------------------------------------------------------------------------
# Todo junto
# -----------------------------------------------------------------------------
@dataclass
class Resultado:
    cfg: Config
    equipo: Equipo
    verificaciones: dict[str, list[Verificacion]]   # por producto, una por T de entrada
    alertas: list[str]   # algo no cumple
    avisos: list[str]    # supuestos, extrapolaciones y datos para otros equipos

    @property
    def cumple(self) -> bool:
        return not self.alertas

    def en_diseno(self, clave: str) -> Verificacion:
        """La verificación del producto a su T de entrada de diseño."""
        T = self.cfg.productos[clave].T_entrada
        return next(v for v in self.verificaciones[clave] if v.T_entrada == T)


def correr(cfg: Config) -> Resultado:
    avisos = [
        "v0: U fijo por sección tomado del caso de yema de Gut et al. (2003); "
        "los tres productos usan el mismo U hasta que v1 meta correlaciones",
    ]
    eq = disenar_equipo(cfg)
    avisos += eq.tubo.avisos
    alertas: list[str] = []

    verifs: dict[str, list[Verificacion]] = {}
    for k, p in cfg.productos.items():
        temps = sorted({p.T_entrada, *p.T_entrada_rango})
        verifs[k] = [verificar_producto(p, eq, cfg, T) for T in temps]
        for v in verifs[k]:
            alertas += v.alertas
            if v.T_entrada == p.T_entrada:  # lo informativo, solo en el punto de diseño
                avisos += v.avisos

    # avisos de propiedades: reología en el tubo y a la salida
    for p in cfg.productos.values():
        avisos += p.reologia(p.T_holding + cfg.perdida_holding / 2).avisos
        r = eq.tubo.elegido.productos[p.clave]
        if r.logs_min < cfg.logs_objetivo:
            avisos.append(
                f"{p.clave}: la partícula más rápida logra {r.logs_min:.1f} log en el tubo "
                f"(objetivo {cfg.logs_objetivo:.0f}); el diseño se apoya en el par T–t "
                "normativo y la letalidad de las placas entra en v2"
            )
    return Resultado(cfg, eq, verifs, list(dict.fromkeys(alertas)), list(dict.fromkeys(avisos)))
