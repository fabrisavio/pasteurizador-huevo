"""Tubo de retención (holding tube).

Es un reactor de flujo pistón imperfecto: el criterio legal es que TODA
partícula pase al menos t_holding a T_holding (9 CFR 590.570), así que el
largo lo fija la partícula más rápida, no la media.

Fluidos ley de potencia: el Reynolds se generaliza con Metzner–Reed (1955), que
define la viscosidad aparente en la pared para que sigan valiendo f = 16/Re y
la transición en ~2100. Con n = 1 todo vuelve al caso newtoniano.

Dos geometrías:
- recto: tramos rectos; las curvas de unión se desprecian.
- espiral: tubo enrollado en hélice. En cada vuelta, la fuerza centrífuga
  arma dos remolinos (vórtices de Dean) que suben la fricción y retrasan la
  turbulencia. También angostan la distribución de tiempos de residencia, pero
  cuánto se mide con un ensayo de trazador: por eso el diseño sigue usando el
  factor de tubo recto y la mejora posible se informa como sensibilidad.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .cinetica import Cinetica


@dataclass
class ProductoTubo:
    """Lo que el tubo necesita saber de cada producto."""
    clave: str
    m: float          # kg/s
    rho: float        # kg/m³ a la T media del tubo
    K: float          # Pa·sⁿ
    n: float
    T_holding: float  # °C, T legal (salida del tubo)
    t_holding: float  # s
    cinetica: Cinetica


@dataclass
class Criterios:
    Re_critico: float = 2100.0
    Re_turbulento: float = 4000.0
    factor_turbulento_vmax: float = 1.25


@dataclass(frozen=True)
class Geometria:
    tipo: str = "recto"            # recto | espiral
    D_espiral: float | None = None  # m, diámetro de la hélice (al eje del tubo)
    separacion: float = 0.01        # m, luz entre vueltas

    def __post_init__(self):
        if self.tipo not in ("recto", "espiral"):
            raise ValueError(f"Geometría de tubo desconocida: {self.tipo}")
        if self.tipo == "espiral" and not self.D_espiral:
            raise ValueError("La espiral necesita su diámetro (diam_espiral_m)")

    @property
    def es_espiral(self) -> bool:
        return self.tipo == "espiral"


RECTO = Geometria()


@dataclass
class Flujo:
    v_media: float      # m/s
    mu_ef: float        # Pa·s, viscosidad aparente de Metzner–Reed
    Re: float           # Reynolds generalizado
    regimen: str        # laminar | transicion | turbulento
    factor: float       # v_max / v_media usado para el largo
    f_fanning: float    # ya corregido por curvatura si es espiral
    gamma_pared: float  # 1/s, tasa de corte en la pared (laminar)
    Dean: float = 0.0   # Re·√(D/Dc); 0 en tubo recto
    correccion_f: float = 1.0  # f_espiral / f_recto


# -----------------------------------------------------------------------------
# Correlaciones
# -----------------------------------------------------------------------------
def viscosidad_metzner_reed(K: float, n: float, v: float, D: float) -> float:
    """μ_ef = K·((3n+1)/(4n))ⁿ·(8v/D)^(n−1).  Con n = 1 da μ = K."""
    return K * ((3 * n + 1) / (4 * n)) ** n * (8 * v / D) ** (n - 1)


def correccion_white_laminar(Dean: float) -> float:
    """f_espiral/f_recto en laminar, White (1929). Vale 1 para Dean ≤ 11,6."""
    if Dean <= 11.6:
        return 1.0
    return 1.0 / (1.0 - (1.0 - (11.6 / Dean) ** 0.45) ** (1 / 0.45))


def correccion_ito_turbulento(Re: float, D: float, D_espiral: float) -> float:
    """f_espiral/f_recto en turbulento, Ito (1959): [Re·(D/Dc)²]^(1/20)."""
    x = Re * (D / D_espiral) ** 2
    return x ** 0.05 if x > 6 else 1.0


def Re_critico_espiral(D: float, D_espiral: float) -> float:
    """Ito (1959): la curvatura estabiliza el laminar hasta 2·10⁴·(D/Dc)^0,32."""
    return 2e4 * (D / D_espiral) ** 0.32


# -----------------------------------------------------------------------------
# Flujo de un producto en un diámetro
# -----------------------------------------------------------------------------
def flujo_en_tubo(p: ProductoTubo, D: float, crit: Criterios,
                  geo: Geometria = RECTO) -> Flujo:
    A = math.pi * D**2 / 4
    v = p.m / (p.rho * A)
    mu_ef = viscosidad_metzner_reed(p.K, p.n, v, D)
    Re = p.rho * v * D / mu_ef
    factor_laminar = (3 * p.n + 1) / (p.n + 1)  # perfil de ley de potencia; 2 si n = 1

    Re_c, Re_t = crit.Re_critico, crit.Re_turbulento
    if geo.es_espiral:
        Re_c = max(Re_c, Re_critico_espiral(D, geo.D_espiral))
        Re_t = Re_c * crit.Re_turbulento / crit.Re_critico  # misma banda de transición, corrida

    if Re < Re_c:
        regimen, factor = "laminar", factor_laminar
    elif Re < Re_t:
        # En transición el perfil no está desarrollado como turbulento:
        # se diseña con el factor laminar, que es el lado seguro.
        regimen, factor = "transicion", factor_laminar
    else:
        regimen, factor = "turbulento", crit.factor_turbulento_vmax

    # Fricción: laminar exacta; transición y turbulento con Blasius (conservador en transición)
    f = 16.0 / Re if regimen == "laminar" else 0.079 / Re**0.25
    Dean, corr = 0.0, 1.0
    if geo.es_espiral:
        Dean = Re * math.sqrt(D / geo.D_espiral)
        corr = (correccion_white_laminar(Dean) if regimen == "laminar"
                else correccion_ito_turbulento(Re, D, geo.D_espiral))
    gamma = (3 * p.n + 1) / (4 * p.n) * 8 * v / D
    return Flujo(v, mu_ef, Re, regimen, factor, f * corr, gamma, Dean, corr)


def caida_presion(p: ProductoTubo, flujo: Flujo, D: float, L: float) -> float:
    """Δp = 4f·(L/D)·ρv²/2, en Pa."""
    return 4 * flujo.f_fanning * (L / D) * p.rho * flujo.v_media**2 / 2


# -----------------------------------------------------------------------------
# Un diámetro con todos los productos
# -----------------------------------------------------------------------------
@dataclass
class ResultadoProductoTubo:
    clave: str
    flujo: Flujo
    L_requerido: float   # m: lo que este producto solo necesitaría
    dp: float            # Pa, con el largo instalado
    t_rapida_min: float  # residencia de la partícula más rápida
    t_media_min: float   # residencia media
    logs_min: float      # logs de la partícula más rápida
    logs_medio: float


@dataclass
class OpcionDN:
    DN: int
    D: float
    L: float             # m, el que gobierna (máximo entre productos)
    volumen_L: float
    gobierna: str
    dp_max: float        # Pa, el peor producto
    Q_CIP_m3_h: float
    productos: dict[str, ResultadoProductoTubo]
    D_ext: float | None = None
    vueltas: float | None = None   # solo espiral
    alto: float | None = None      # m, solo espiral
    factible: bool = True
    motivo: str = ""


def _forma_espiral(L: float, D_ext: float | None, geo: Geometria):
    if not geo.es_espiral:
        return None, None
    vueltas = L / (math.pi * geo.D_espiral)
    alto = vueltas * (D_ext + geo.separacion) if D_ext else None
    return vueltas, alto


def evaluar_DN(DN: int, D: float, productos: list[ProductoTubo], crit: Criterios,
               v_CIP: float, geo: Geometria = RECTO, D_ext: float | None = None) -> OpcionDN:
    flujos = {p.clave: flujo_en_tubo(p, D, crit, geo) for p in productos}
    L_req = {p.clave: flujos[p.clave].factor * flujos[p.clave].v_media * p.t_holding
             for p in productos}
    gobierna = max(L_req, key=L_req.get)
    L = L_req[gobierna]  # un solo tubo físico: manda el que más largo pide

    res = {}
    for p in productos:
        fl = flujos[p.clave]
        t_rap = L / (fl.factor * fl.v_media) / 60
        t_med = L / fl.v_media / 60
        res[p.clave] = ResultadoProductoTubo(
            clave=p.clave,
            flujo=fl,
            L_requerido=L_req[p.clave],
            dp=caida_presion(p, fl, D, L),
            t_rapida_min=t_rap,
            t_media_min=t_med,
            logs_min=p.cinetica.logs(t_rap, p.T_holding),
            logs_medio=p.cinetica.logs(t_med, p.T_holding),
        )
    A = math.pi * D**2 / 4
    vueltas, alto = _forma_espiral(L, D_ext, geo)
    return OpcionDN(
        DN=DN, D=D, L=L, volumen_L=A * L * 1000, gobierna=gobierna,
        dp_max=max(r.dp for r in res.values()),
        Q_CIP_m3_h=v_CIP * A * 3600, productos=res,
        D_ext=D_ext, vueltas=vueltas, alto=alto,
    )


# -----------------------------------------------------------------------------
# Sensibilidad: ¿qué pasa si el ensayo mide una partícula rápida más lenta?
# -----------------------------------------------------------------------------
@dataclass
class PuntoSensibilidad:
    factor: float        # v_max / v_media supuesto (igual para todos los productos)
    L: float
    volumen_L: float
    dp_max: float        # Pa
    vueltas: float | None
    alto: float | None


def sensibilidad_factor(op: OpcionDN, productos: list[ProductoTubo], factores: list[float],
                        geo: Geometria) -> list[PuntoSensibilidad]:
    """Largo y Δp del DN elegido si el factor de la partícula rápida fuera otro.

    El Δp de cada producto escala con el largo (f no depende de L), así que no
    hace falta recalcular el flujo.
    """
    t = {p.clave: p.t_holding for p in productos}
    puntos = []
    for fac in factores:
        L = max(fac * r.flujo.v_media * t[k] for k, r in op.productos.items())
        dp = max(r.dp * L / op.L for r in op.productos.values())
        vueltas, alto = _forma_espiral(L, op.D_ext, geo)
        puntos.append(PuntoSensibilidad(fac, L, math.pi * op.D**2 / 4 * L * 1000, dp, vueltas, alto))
    return puntos


# -----------------------------------------------------------------------------
# Elección del DN
# -----------------------------------------------------------------------------
@dataclass
class Tubo:
    elegido: OpcionDN
    opciones: list[OpcionDN]
    geometria: Geometria = RECTO
    sensibilidad: list[PuntoSensibilidad] = field(default_factory=list)
    avisos: list[str] = field(default_factory=list)


def elegir_tubo(productos: list[ProductoTubo], DN_diam: dict[int, float],
                crit: Criterios, dp_max_Pa: float, v_CIP: float,
                geo: Geometria = RECTO, DN_ext: dict[int, float] | None = None,
                factores_sensibilidad: list[float] | None = None) -> Tubo:
    """Evalúa todos los DN y elige el más chico que respeta el Δp máximo.

    Por qué el más chico: en laminar el volumen no depende del DN, así que el
    DN solo reparte ese volumen entre largo y Δp (Δp ∝ 1/D⁶). El Δp máximo pone
    la cota inferior; achicar el DN baja el caudal de CIP y sube la velocidad.
    """
    DN_ext = DN_ext or {}
    opciones = [evaluar_DN(DN, D, productos, crit, v_CIP, geo, DN_ext.get(DN))
                for DN, D in sorted(DN_diam.items())]
    for op in opciones:
        if op.dp_max > dp_max_Pa:
            op.factible = False
            op.motivo = f"Δp {op.dp_max/1e5:.2f} bar > {dp_max_Pa/1e5:.2f} bar"

    avisos = []
    factibles = [op for op in opciones if op.factible]
    if factibles:
        elegido = factibles[0]
    else:
        elegido = opciones[-1]
        avisos.append("Ningún DN cumple el Δp máximo; se toma el más grande disponible")

    margen = 1 - elegido.dp_max / dp_max_Pa
    if elegido.factible and margen < 0.10:
        avisos.append(
            f"Tubo: DN{elegido.DN} cumple el Δp máximo con un margen de {margen:.0%}; "
            "una variación chica en la reología lo deja afuera"
        )
    for clave, r in elegido.productos.items():
        if r.flujo.regimen == "transicion":
            avisos.append(
                f"{clave}: tubo en transición (Re = {r.flujo.Re:.0f}); se usó el factor "
                f"laminar {r.flujo.factor:.2f} para la partícula más rápida"
            )

    if geo.es_espiral:
        relacion = geo.D_espiral / elegido.D
        if not 15 <= relacion <= 860:
            avisos.append(
                f"Tubo: Dc/D = {relacion:.0f}, fuera del rango de las correlaciones de Ito (15–860)"
            )
        avisos.append(
            "Tubo en espiral: corrección de fricción por curvatura (White 1929 en laminar, "
            "Ito 1959 en turbulento) pensada para newtonianos, aplicada con el Re de "
            "Metzner–Reed: aproximación"
        )
        avisos.append(
            "Tubo en espiral: el largo se diseña con el factor de partícula rápida de tubo "
            "recto (lado seguro). El factor real se mide con un ensayo de trazador "
            "(Kaufman et al., 1968); la tabla de sensibilidad muestra lo que podría ahorrar"
        )

    sens = sensibilidad_factor(elegido, productos, factores_sensibilidad or [], geo)
    return Tubo(elegido=elegido, opciones=opciones, geometria=geo,
                sensibilidad=sens, avisos=avisos)
