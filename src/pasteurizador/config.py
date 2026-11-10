"""Carga la configuración YAML y la pasa a SI (salvo T, que queda en °C).

El YAML está en unidades prácticas (kg/h, min, bar) porque así se leen los
datos de planta. Acá se convierten una sola vez, al cargar: de ahí en adelante
todo el código trabaja en SI y no hay que pensar en unidades.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml

from . import propiedades as pr
from .aislacion import Aislacion
from .cinetica import Cinetica
from .dimensiones import Bastidor
from .holding import Criterios, Geometria
from .placas import CORRELACIONES, Placa


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
    # Rango de tasa de corte donde se midió la reología (aviso si en las placas se sale)
    gamma_rango: tuple[float, float] | None = None

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
    # --- placas (v1) ---
    placa: Placa | None = None
    modelo_U: str = "fijo"        # fijo (v0) | correlacion (v1)
    correlacion: str = "kumar"
    correlacion_contraste: str | None = None
    Rf_regen: float = 0.0         # m²·K/W, producto, cada lado de la regeneración
    Rf_calent: float = 0.0        # m²·K/W, producto en el calentamiento
    Rf_agua: float = 0.0          # m²·K/W, agua caliente
    pasos_candidatos: tuple[int, ...] = (1,)
    dp_max_placas: float = 1e5    # Pa, por sección y por corriente
    Re_CIP_placas: float | None = None
    T_CIP: float = 70.0           # °C
    bastidor: Bastidor = field(default_factory=lambda: Bastidor())  # para el plano
    DN_fijo: int | None = None    # None: el DN más chico que cumple el Δp

    @property
    def con_correlacion(self) -> bool:
        return self.modelo_U == "correlacion"

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
        gamma_rango=tuple(d["reologia_gamma_1_s"]) if "reologia_gamma_1_s" in d else None,
    )


def _placa(d: dict) -> Placa:
    return Placa(
        area=d["area_m2"],
        ancho=d["ancho_m"],
        b=d["separacion_b_m"],
        angulo=d["angulo_chevron_grados"],
        phi=d["factor_area"],
        espesor=d["espesor_m"],
        k=d["k_W_mK"],
        D_puerto=d.get("diam_puerto_m", 0.05),
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
        curvas_aisladas=bool(d.get("curvas_aisladas", True)),
        h_desnudo=d.get("h_desnudo_W_m2K", d.get("h_ext_W_m2K", 10.0)),
        union=d.get("union_sin_aislar_m", 0.05),
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
    modelo_U = eq.get("modelo_U", "fijo")
    if modelo_U not in ("fijo", "correlacion"):
        raise ConfigInvalida(f"equipo.modelo_U = {modelo_U}: tiene que ser 'fijo' o 'correlacion'")
    for clave in ("correlacion", "correlacion_contraste"):
        if eq.get(clave) not in CORRELACIONES | {None}:
            raise ConfigInvalida(f"equipo.{clave} = {eq.get(clave)}: opciones {sorted(CORRELACIONES)}")
    ens = eq.get("ensuciamiento_m2K_W", {})
    arr = eq.get("arreglo", {})
    cip = eq.get("cip", {})
    if "aislacion" not in hold and "perdida_T_C" not in hold:
        raise ConfigInvalida("equipo.holding necesita 'aislacion' (se calcula la pérdida) "
                             "o 'perdida_T_C' (valor fijo)")
    aislacion = _aislacion(hold.get("aislacion"))
    # Con coquilla, cada tramo es más gordo: eso define cuánto separar las horquillas
    extra = 2 * aislacion.espesor if aislacion and aislacion.forma == "tubo" else 0.0
    bas = eq.get("bastidor", {})
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
            largo_tramo=geo.get("largo_tramo_max_m", 6.0),
            radio_curva=geo.get("radio_curva_m", 0.06),
            luz=geo.get("luz_entre_tubos_m", 0.02),
            extra_envolvente=extra,
            K1_curva=geo.get("K1_curva", 1000.0),
            Kinf_curva=geo.get("Kinf_curva", 0.30),
            pendiente=geo.get("pendiente_min", 0.02),
        ),
        aislacion=aislacion,
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
        placa=_placa(eq["placa"]),
        modelo_U=modelo_U,
        correlacion=eq.get("correlacion", "kumar"),
        correlacion_contraste=eq.get("correlacion_contraste"),
        Rf_regen=ens.get("producto_regeneracion", 0.0),
        Rf_calent=ens.get("producto_calentamiento", 0.0),
        Rf_agua=ens.get("agua_caliente", 0.0),
        pasos_candidatos=tuple(arr.get("pasos_candidatos", [1])),
        dp_max_placas=arr.get("dp_max_seccion_bar", 1.0) * 1e5,
        Re_CIP_placas=cip.get("Re_min_placas"),
        T_CIP=cip.get("T_C", 70.0),
        bastidor=Bastidor(
            alto=bas.get("alto_m", 0.92),
            ancho=bas.get("ancho_m", 0.32),
            largo_rango=tuple(bas.get("largo_rango_m", [0.515, 1.43])),
            conector=bas.get("conector_m", 0.10),
            extremos_por_seccion=bas.get("placas_extremo_por_seccion", 2),
        ),
        DN_fijo=hold.get("DN_fijo"),
    )
