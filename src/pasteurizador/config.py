"""Carga la configuración YAML y la pasa a SI (salvo T, que queda en °C).

El YAML está en unidades prácticas (kg/h, min, bar) porque así se leen los
datos de planta. Acá se convierten una sola vez, al cargar: de ahí en adelante
todo el código trabaja en SI y no hay que pensar en unidades.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

from . import propiedades as pr
from .aislacion import Aislacion
from .cinetica import Cinetica
from .holding import Criterios, Geometria


@dataclass
class Producto:
    clave: str
    nombre: str
    m: float                      # kg/s
    T_entrada: float              # °C, valor de diseño
    T_entrada_rango: tuple[float, float]  # °C, se verifica en las dos puntas
    T_holding: float              # °C
    t_holding: float              # s
    T_salida: float               # °C, objetivo hacia el secador
    T_max_agua: float             # °C, tope del agua caliente (protege la pared)
    composicion: dict[str, float]
    modelo_termicas: str
    modelo_reologia: str
    cinetica: Cinetica
    horas_anio: float
    pH: float | None = None        # informativo, sale en el reporte
    T_espesa: float | None = None  # °C donde empieza a subir la viscosidad (aviso)

    # Métodos: funciones que viven dentro de la ficha. p.termicas(60) da ρ, cp, k
    # del producto p a 60 °C sin tener que pasarle su composición cada vez.
    def termicas(self, T: float) -> pr.Termicas:
        return pr.termicas(self.modelo_termicas, self.composicion, T)

    def reologia(self, T: float) -> pr.Reologia:
        return pr.reologia(self.modelo_reologia, T)


@dataclass
class Config:
    productos: dict[str, Producto]
    area_placa: float             # m²
    U_regen: float                # W/(m²·K)
    U_calent: float
    DN_diam: dict[int, float]     # DN → diámetro interno [m]
    DN_ext: dict[int, float]      # DN → diámetro externo [m]
    geometria_tubo: Geometria
    # Pérdida de T en el tubo: se calcula con la aislación; si no hay sección
    # de aislación en el YAML, se usa el valor fijo perdida_T_C.
    aislacion: Aislacion | None
    perdida_fija: float | None    # °C
    margen_control: float         # °C sobre la T legal a la que opera el control
    dp_max_holding: float         # Pa
    factor_dp: float              # el Δp calculado × esto tiene que entrar en el máximo
    v_CIP: float                  # m/s
    DT_agua: float                # °C, salto del agua caliente en el diseño
    criterios: Criterios
    logs_objetivo: float

    def T_operacion(self, p: Producto) -> float:
        """°C a la salida del tubo en operación: la T legal más el margen del control."""
        return p.T_holding + self.margen_control


class ConfigInvalida(ValueError):
    """El archivo existe pero no es una configuración de diseño."""


def _producto(clave: str, d: dict) -> Producto:
    # d.get("x", defecto): si la clave no está en el YAML, usa el defecto
    rango = d.get("T_entrada_rango_C", [d["T_entrada_C"], d["T_entrada_C"]])
    return Producto(
        clave=clave,
        nombre=d.get("nombre", clave),
        m=d["caudal_kg_h"] / 3600.0,
        T_entrada=float(d["T_entrada_C"]),
        T_entrada_rango=(float(rango[0]), float(rango[1])),
        T_holding=float(d["T_holding_C"]),
        t_holding=d["t_holding_min"] * 60.0,
        T_salida=float(d["T_salida_C"]),
        T_max_agua=float(d["T_max_agua_C"]),
        composicion={k: float(v) for k, v in d["composicion"].items()},
        modelo_termicas=d.get("termicas", "choi_okos"),
        modelo_reologia=d["reologia"],
        cinetica=Cinetica(**d["cinetica"]),  # ** reparte el diccionario en argumentos
        horas_anio=float(d["horas_anio"]),
        pH=d.get("pH"),
        T_espesa=d.get("T_espesa_C"),
    )


def _aislacion(d: dict | None) -> Aislacion | None:
    if not d:
        return None
    return Aislacion(
        forma=d.get("forma", "bloque"),
        espesor=d["espesor_m"],
        k=d["k_W_mK"],
        h_ext=d.get("h_ext_W_m2K", 10.0),
        T_amb=d.get("T_ambiente_C", 20.0),
        holgura=d.get("holgura_m", 0.05),
    )


def cargar(ruta: str | Path) -> Config:
    with open(ruta, encoding="utf-8") as f:
        c = yaml.safe_load(f)  # el YAML queda como diccionarios y listas anidados
    if not isinstance(c, dict) or "productos" not in c:
        raise ConfigInvalida(
            f"{Path(ruta).name} no tiene la sección 'productos': no es una configuración "
            "de diseño (¿es el caso de validación? se corre con --validar)"
        )
    eq, hold, crit = c["equipo"], c["equipo"]["holding"], c["criterios"]
    geo = hold.get("geometria", {})
    if "aislacion" not in hold and "perdida_T_C" not in hold:
        raise ConfigInvalida("equipo.holding necesita 'aislacion' (se calcula la pérdida) "
                             "o 'perdida_T_C' (valor fijo)")
    return Config(
        productos={k: _producto(k, v) for k, v in c["productos"].items()},
        area_placa=eq["placa"]["area_m2"],
        U_regen=eq["U_fijo_W_m2K"]["regeneracion"],
        U_calent=eq["U_fijo_W_m2K"]["calentamiento"],
        DN_diam={int(k): float(v) for k, v in hold["DN_diam_int_m"].items()},
        DN_ext={int(k): float(v) for k, v in hold.get("DN_diam_ext_m", {}).items()},
        geometria_tubo=Geometria(
            tipo=geo.get("tipo", "recto"),
            D_espiral=geo.get("diam_espiral_m"),
            separacion=geo.get("separacion_vueltas_m", 0.01),
            t_min_sobre_t_medio=geo.get("t_min_sobre_t_medio", 0.586),
            Dean_min=geo.get("Dean_min", 11.6),
        ),
        aislacion=_aislacion(hold.get("aislacion")),
        perdida_fija=hold.get("perdida_T_C"),
        margen_control=c.get("control", {}).get("margen_sobre_T_legal_C", 0.0),
        dp_max_holding=hold["dp_max_bar"] * 1e5,
        factor_dp=hold.get("factor_seguridad_dp", 1.0),
        v_CIP=hold["v_CIP_min_m_s"],
        DT_agua=c["servicios"]["agua_caliente"]["DT_C"],
        criterios=Criterios(
            Re_critico=crit["Re_critico"],
            Re_turbulento=crit["Re_turbulento"],
            factor_turbulento_vmax=crit["factor_turbulento_vmax"],
        ),
        logs_objetivo=crit["logs_objetivo"],
    )
