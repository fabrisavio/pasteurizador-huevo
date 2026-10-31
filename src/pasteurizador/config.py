"""Carga la configuración YAML y la pasa a SI (salvo T, que queda en °C)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

from . import propiedades as pr
from .cinetica import Cinetica
from .holding import Criterios, Geometria


@dataclass
class Producto:
    clave: str
    nombre: str
    m: float                      # kg/s
    T_entrada: float              # °C, valor de diseño
    T_entrada_rango: tuple[float, float]
    T_holding: float              # °C
    t_holding: float              # s
    T_salida: float               # °C
    T_max_agua: float             # °C
    composicion: dict[str, float]
    modelo_termicas: str
    modelo_reologia: str
    cinetica: Cinetica
    horas_anio: float

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
    factores_sensibilidad: list[float]
    perdida_holding: float        # °C
    dp_max_holding: float         # Pa
    v_CIP: float                  # m/s
    DT_agua: float                # °C
    criterios: Criterios
    logs_objetivo: float


def _producto(clave: str, d: dict) -> Producto:
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
        cinetica=Cinetica(**d["cinetica"]),
        horas_anio=float(d["horas_anio"]),
    )


def cargar(ruta: str | Path) -> Config:
    with open(ruta, encoding="utf-8") as f:
        c = yaml.safe_load(f)
    eq, hold, crit = c["equipo"], c["equipo"]["holding"], c["criterios"]
    return Config(
        productos={k: _producto(k, v) for k, v in c["productos"].items()},
        area_placa=eq["placa"]["area_m2"],
        U_regen=eq["U_fijo_W_m2K"]["regeneracion"],
        U_calent=eq["U_fijo_W_m2K"]["calentamiento"],
        DN_diam={int(k): float(v) for k, v in hold["DN_diam_int_m"].items()},
        DN_ext={int(k): float(v) for k, v in hold.get("DN_diam_ext_m", {}).items()},
        geometria_tubo=Geometria(
            tipo=hold.get("geometria", {}).get("tipo", "recto"),
            D_espiral=hold.get("geometria", {}).get("diam_espiral_m"),
            separacion=hold.get("geometria", {}).get("separacion_vueltas_m", 0.01),
        ),
        factores_sensibilidad=[float(x) for x in hold.get("factores_sensibilidad", [])],
        perdida_holding=hold["perdida_T_C"],
        dp_max_holding=hold["dp_max_bar"] * 1e5,
        v_CIP=hold["v_CIP_min_m_s"],
        DT_agua=c["servicios"]["agua_caliente"]["DT_C"],
        criterios=Criterios(
            Re_critico=crit["Re_critico"],
            Re_turbulento=crit["Re_turbulento"],
            factor_turbulento_vmax=crit["factor_turbulento_vmax"],
        ),
        logs_objetivo=crit["logs_objetivo"],
    )
