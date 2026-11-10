"""Tubo de retención (holding tube).

Qué resuelve: el largo y el diámetro del tubo donde el producto se queda
t_holding a T_holding.

Por qué el largo lo fija la partícula más rápida y no la media: el criterio
legal es que TODA partícula pase al menos t_holding a T_holding
(9 CFR 590.570). En flujo laminar el perfil de velocidad es una parábola: el
centro va al doble de la media (newtoniano). Si el tubo se dimensionara con la
velocidad media, el centro saldría con la mitad del tiempo. Es la diferencia
entre un PFR ideal y uno real.

Fluidos ley de potencia (yema, huevo entero): el Reynolds se generaliza con
Metzner–Reed (1955). Define una viscosidad aparente evaluada en la pared para
que sigan valiendo f = 16/Re y la transición en ~2100. Con n = 1 todo vuelve
al caso newtoniano.

Tres geometrías:
- recto: tramos rectos; las curvas de unión se desprecian.
- horquillas: tramos rectos unidos por curvas de 180°, como un serpentín
  (la forma habitual de un tubo de retención). Para la partícula más rápida
  se toma el factor de tubo recto: las curvas mezclan algo el centro con la
  pared, pero no hay un valor validado para tramos de metros entre curvas, y
  no contarlo es el lado seguro. Las curvas sí suman Δp (Hooper, 1981).
- espiral: tubo enrollado en hélice. En cada vuelta la fuerza centrífuga
  empuja el fluido rápido del centro hacia afuera y arma dos remolinos
  (vórtices de Dean). Eso sube la fricción, retrasa la turbulencia y mezcla
  el centro con la pared: la partícula más rápida ya no viaja siempre por el
  centro y tarda más. Para el diseño se toma t_min/t_medio = 0,586, el valor
  más bajo (el más conservador) de los análisis de flujo laminar en hélice
  de Ruthven (1971) y Trivedi y Vasudeva (1974), según la revisión de
  Gutierrez et al. (2010), que además midieron 0,70 o más en un tubo
  helicoidal de retención.

Construcciones de Python que aparecen acá (y en el resto del paquete):
- @dataclass: una "ficha" con campos con nombre, como una fila de tabla. Python
  arma solo el constructor: ProductoTubo(clave="YP", m=0.3, ...).
- frozen=True: la ficha no se puede modificar después de creada.
- `x: float | None = None`: el campo puede quedar vacío (None) si no aplica.
- `-> Flujo`: anota qué devuelve la función. Es solo documentación.
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
    K: float          # Pa·sⁿ, índice de consistencia (ley de potencia)
    n: float          # índice de flujo: 1 newtoniano, < 1 pseudoplástico
    T_holding: float  # °C, T legal (salida del tubo)
    t_holding: float  # s
    cinetica: Cinetica


@dataclass
class Criterios:
    Re_critico: float = 2100.0          # fin del laminar en tubo recto
    Re_turbulento: float = 4000.0       # comienzo del turbulento desarrollado
    factor_turbulento_vmax: float = 1.25  # v_max/v_media en turbulento


@dataclass(frozen=True)
class Geometria:
    tipo: str = "recto"             # recto | espiral | horquillas
    D_espiral: float | None = None  # m, diámetro de la hélice (al eje del tubo)
    separacion: float = 0.01        # m, luz entre vueltas
    # Partícula más rápida en la hélice: t_min/t_medio. 0,586 → v_max/v_media = 1,71
    t_min_sobre_t_medio: float = 0.586
    # Por debajo de este Dean el flujo secundario es despreciable (White, 1929):
    # el tubo se comporta como recto y se usa el factor de tubo recto.
    Dean_min: float = 11.6
    # --- horquillas: tramos rectos unidos por curvas de 180° ---
    largo_tramo: float = 6.0        # m, largo máximo de cada tramo (los tubos vienen de ~6 m)
    radio_curva: float = 0.06       # m, al eje (curva DIN 11852 de DN40: R = 60 mm)
    luz: float = 0.02               # m, luz mínima entre tubos vecinos, con su aislación
    # m que la aislación le suma al diámetro exterior de cada tramo (2 × espesor
    # si se aísla con coquilla): define cuánto hay que separar los tramos.
    extra_envolvente: float = 0.0
    # Hooper (1981), curva de 180° de radio largo (R/D ≈ 1,5): K = K1/Re + K∞·(1 + 1/D[pulg])
    K1_curva: float = 1000.0
    Kinf_curva: float = 0.30
    # --- todas: subida continua mínima en el sentido del flujo, curvas
    # incluidas, desde el punto más bajo hasta la válvula de desvío (CFIA) ---
    pendiente: float = 0.02

    def __post_init__(self):
        # Se ejecuta al crear la ficha: valida antes de que el error aparezca lejos.
        if self.tipo not in ("recto", "espiral", "horquillas"):
            raise ValueError(f"Geometría de tubo desconocida: {self.tipo} (recto | espiral | horquillas)")
        if self.tipo == "espiral" and not self.D_espiral:
            raise ValueError("La espiral necesita su diámetro (diam_espiral_m)")
        if self.tipo == "horquillas" and (self.largo_tramo <= 0 or self.radio_curva <= 0):
            raise ValueError("Horquillas: largo_tramo_m y radio_curva_m tienen que ser positivos")
        if not 0.5 <= self.t_min_sobre_t_medio <= 1.0:
            # 0,5 es el laminar newtoniano recto; 1 es flujo pistón ideal
            raise ValueError("t_min/t_medio tiene que estar entre 0,5 y 1")

    @property
    def es_espiral(self) -> bool:
        # @property: se usa como un campo (geo.es_espiral), sin paréntesis
        return self.tipo == "espiral"

    @property
    def es_horquillas(self) -> bool:
        return self.tipo == "horquillas"

    @property
    def factor_espiral(self) -> float:
        """v_max/v_media de diseño en la hélice."""
        return 1.0 / self.t_min_sobre_t_medio


RECTO = Geometria()  # geometría por defecto, para no repetirla en cada llamada


@dataclass
class Horquillas:
    """Cómo queda armado un tubo en horquillas, y cuánto lugar ocupa."""
    tramos: int
    largo_tramo: float     # m, cada tramo recto (todos iguales)
    curvas: int            # curvas de 180° (una menos que los tramos)
    largo_curva: float     # m, desarrollo de cada curva: π·R más el niple si hace falta
    niple: float           # m, tramo recto dentro de la curva para separar los planos
    entre_planos: float    # m, distancia entre ejes de los dos planos
    pendiente: float       # m/m, la que resulta (≥ la mínima pedida)
    desnivel: float        # m, de la entrada a la salida del tubo
    D_envolvente: float    # m, tubo + aislación
    largo: float           # m, envolvente del tubo con su aislación (en planta)
    ancho: float           # m
    alto: float            # m, del eje del tramo más bajo a la cara superior del más alto


def disponer_horquillas(L: float, D_ext: float, geo: Geometria) -> Horquillas:
    """Arma L metros de tubo en horquillas con subida continua.

    Disposición: dos planos verticales, uno adelante y otro atrás. Los tramos
    alternan de plano y cada curva de 180° pasa de uno al otro, así el tubo
    da vueltas alrededor de un rectángulo largo y angosto, subiendo siempre
    (una hélice aplastada). Es la forma de menor altura posible: si todo el
    tubo sube al menos la pendiente mínima, el desnivel total no puede ser
    menor que pendiente × L, se lo arme como se lo arme.

    Dos tramos seguidos del mismo plano quedan separados en altura por dos
    subidas; si eso no alcanza para que entren los tubos con su aislación,
    se sube la pendiente. Entre planos, si el diámetro de la curva (2R) es
    menor que el tubo aislado más la luz, la curva se arma con dos codos de
    90° y un niple recto en el medio.
    """
    D_env = D_ext + geo.extra_envolvente
    R = geo.radio_curva
    entre = max(2 * R, D_env + geo.luz)
    niple = entre - 2 * R
    Lc = math.pi * R + niple
    # L = N·Lt + (N−1)·Lc con Lt ≤ largo_tramo → el N más chico que alcanza
    N = max(1, math.ceil((L + Lc) / (geo.largo_tramo + Lc) - 1e-9))
    Lt = (L - (N - 1) * Lc) / N
    s = geo.pendiente
    if N > 2:
        s = max(s, (D_env + geo.luz) / (2 * (Lt + Lc)))
    desnivel = s * L
    if N > 1:
        largo, ancho = Lt + 2 * (R + D_env / 2), entre + D_env
    else:
        largo, ancho = Lt, D_env
    return Horquillas(N, Lt, N - 1, Lc, niple, entre, s, desnivel, D_env,
                      largo, ancho, desnivel + D_env)


def K_curva(Re: float, D: float, geo: Geometria) -> float:
    """Coeficiente de pérdida de una curva de 180°, método 2-K de Hooper (1981).

    K = K1/Re + K∞·(1 + 1/D), con D en pulgadas. El primer término pesa en
    laminar (la yema, Re ~ 170, da K ≈ 6) y el segundo en turbulento.
    """
    return geo.K1_curva / Re + geo.Kinf_curva * (1 + 0.0254 / D)


def pendiente_espiral(D_ext: float, geo: Geometria) -> float:
    """Lo que sube la hélice por metro de tubo: un paso por cada vuelta de π·Dc."""
    return (D_ext + geo.separacion) / (math.pi * geo.D_espiral)


@dataclass
class Flujo:
    v_media: float      # m/s
    mu_ef: float        # Pa·s, viscosidad aparente de Metzner–Reed
    Re: float           # Reynolds generalizado
    regimen: str        # laminar | transicion | turbulento
    factor: float       # v_max/v_media usado para el largo (diseño)
    factor_recto: float  # el que correspondería en tubo recto, para comparar
    f_fanning: float    # ya corregido por curvatura si es espiral
    gamma_pared: float  # 1/s, tasa de corte en la pared (laminar)
    Dean: float = 0.0   # Re·√(D/Dc); 0 en tubo recto
    correccion_f: float = 1.0  # f_espiral / f_recto


# -----------------------------------------------------------------------------
# Correlaciones
# -----------------------------------------------------------------------------
def viscosidad_metzner_reed(K: float, n: float, v: float, D: float) -> float:
    """μ_ef = K·((3n+1)/(4n))ⁿ·(8v/D)^(n−1).  Con n = 1 da μ = K.

    8v/D es la tasa de corte "nominal" en la pared. Un pseudoplástico (n < 1)
    se vuelve más fluido cuanto más rápido lo cortás: por eso el exponente n−1
    es negativo y μ_ef baja al subir la velocidad.
    """
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
    """Régimen, factor de partícula rápida y fricción de un producto en un tubo de diámetro D."""
    A = math.pi * D**2 / 4
    v = p.m / (p.rho * A)                       # caudal másico → velocidad media
    mu_ef = viscosidad_metzner_reed(p.K, p.n, v, D)
    Re = p.rho * v * D / mu_ef
    # Perfil laminar de ley de potencia: v_max/v_media = (3n+1)/(n+1).
    # n = 1 → 2 (la parábola de siempre); n → 0 → 1 (flujo pistón).
    factor_laminar = (3 * p.n + 1) / (p.n + 1)

    # Límites de régimen. En espiral la curvatura estabiliza el laminar,
    # así que toda la banda de transición se corre hacia Re más altos.
    Re_c, Re_t = crit.Re_critico, crit.Re_turbulento
    if geo.es_espiral:
        Re_c = max(Re_c, Re_critico_espiral(D, geo.D_espiral))
        Re_t = Re_c * crit.Re_turbulento / crit.Re_critico

    if Re < Re_c:
        regimen, factor_recto = "laminar", factor_laminar
    elif Re < Re_t:
        # En transición el perfil no está desarrollado como turbulento:
        # se diseña con el factor laminar, que es el lado seguro.
        regimen, factor_recto = "transicion", factor_laminar
    else:
        regimen, factor_recto = "turbulento", crit.factor_turbulento_vmax

    # Fricción de Fanning: laminar exacta (16/Re); transición y turbulento con
    # Blasius (conservador en transición). Δp = 4f·(L/D)·ρv²/2.
    f = 16.0 / Re if regimen == "laminar" else 0.079 / Re**0.25

    factor, Dean, corr = factor_recto, 0.0, 1.0
    if geo.es_espiral:
        Dean = Re * math.sqrt(D / geo.D_espiral)  # Re corregido por la curvatura
        corr = (correccion_white_laminar(Dean) if regimen == "laminar"
                else correccion_ito_turbulento(Re, D, geo.D_espiral))
        if regimen != "turbulento" and Dean > geo.Dean_min:
            # Los vórtices de Dean sacan a la partícula del centro: factor menor.
            # min(): si el perfil ya es más plano que eso (n bajo), queda el recto.
            factor = min(factor_laminar, geo.factor_espiral)
    gamma = (3 * p.n + 1) / (4 * p.n) * 8 * v / D  # corte real en la pared (Rabinowitsch)
    return Flujo(v, mu_ef, Re, regimen, factor, factor_recto, f * corr, gamma, Dean, corr)


def caida_presion(p: ProductoTubo, flujo: Flujo, D: float, L: float) -> float:
    """Δp = 4f·(L/D)·ρv²/2, en Pa (Fanning; el 4 lo pasa a Darcy)."""
    return 4 * flujo.f_fanning * (L / D) * p.rho * flujo.v_media**2 / 2


# -----------------------------------------------------------------------------
# Un diámetro con todos los productos
# -----------------------------------------------------------------------------
@dataclass
class ResultadoProductoTubo:
    clave: str
    flujo: Flujo
    L_requerido: float   # m: lo que este producto solo necesitaría
    dp: float            # Pa, con el largo instalado (con curvas, si hay)
    t_rapida_min: float  # residencia de la partícula más rápida
    t_media_min: float   # residencia media
    logs_min: float      # logs de la partícula más rápida
    logs_medio: float
    dp_curvas: float = 0.0  # Pa, la parte del Δp que ponen las curvas


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
    horquillas: Horquillas | None = None  # solo horquillas
    factible: bool = True
    motivo: str = ""


def _forma_espiral(L: float, D_ext: float | None, geo: Geometria):
    """Vueltas y altura de la hélice. Cada vuelta mide π·Dc de largo."""
    if not geo.es_espiral:
        return None, None
    vueltas = L / (math.pi * geo.D_espiral)
    # cada vuelta ocupa en altura un diámetro exterior más la luz entre vueltas
    alto = vueltas * (D_ext + geo.separacion) if D_ext else None
    return vueltas, alto


def evaluar_DN(DN: int, D: float, productos: list[ProductoTubo], crit: Criterios,
               v_CIP: float, geo: Geometria = RECTO, D_ext: float | None = None) -> OpcionDN:
    """Largo, volumen, Δp y tiempos de todos los productos en un DN dado."""
    # {clave: valor for p in lista} arma un diccionario recorriendo la lista
    flujos = {p.clave: flujo_en_tubo(p, D, crit, geo) for p in productos}
    # Largo que pide cada producto: L = v_max · t = factor · v_media · t
    L_req = {p.clave: flujos[p.clave].factor * flujos[p.clave].v_media * p.t_holding
             for p in productos}
    gobierna = max(L_req, key=L_req.get)  # la clave con el mayor valor
    L = L_req[gobierna]  # un solo tubo físico: manda el que más largo pide
    # En horquillas, L incluye el desarrollo de las curvas: el producto sigue
    # a temperatura ahí adentro y cuenta para el tiempo.
    horq = disponer_horquillas(L, D_ext or D, geo) if geo.es_horquillas else None

    res = {}
    for p in productos:
        fl = flujos[p.clave]
        t_rap = L / (fl.factor * fl.v_media) / 60
        t_med = L / fl.v_media / 60
        # Curvas: K de Hooper por la energía cinética. La fricción del largo
        # total ya incluye el de las curvas: se cuenta dos veces (lado seguro).
        dp_c = (horq.curvas * K_curva(fl.Re, D, geo) * p.rho * fl.v_media**2 / 2
                if horq else 0.0)
        res[p.clave] = ResultadoProductoTubo(
            clave=p.clave,
            flujo=fl,
            L_requerido=L_req[p.clave],
            dp=caida_presion(p, fl, D, L) + dp_c,
            t_rapida_min=t_rap,
            t_media_min=t_med,
            logs_min=p.cinetica.logs(t_rap, p.T_holding),
            logs_medio=p.cinetica.logs(t_med, p.T_holding),
            dp_curvas=dp_c,
        )
    A = math.pi * D**2 / 4
    vueltas, alto = _forma_espiral(L, D_ext, geo)
    return OpcionDN(
        DN=DN, D=D, L=L, volumen_L=A * L * 1000, gobierna=gobierna,
        dp_max=max(r.dp for r in res.values()),
        Q_CIP_m3_h=v_CIP * A * 3600, productos=res,
        D_ext=D_ext, vueltas=vueltas, alto=alto, horquillas=horq,
    )


# -----------------------------------------------------------------------------
# Elección del DN
# -----------------------------------------------------------------------------
@dataclass
class Tubo:
    elegido: OpcionDN
    opciones: list[OpcionDN]
    geometria: Geometria = RECTO
    # El mismo DN en tramos rectos, para mostrar cuánto ahorra enrollarlo
    equivalente_recto: OpcionDN | None = None
    avisos: list[str] = field(default_factory=list)  # lista vacía nueva por cada Tubo


def elegir_tubo(productos: list[ProductoTubo], DN_diam: dict[int, float],
                crit: Criterios, dp_max_Pa: float, v_CIP: float,
                geo: Geometria = RECTO, DN_ext: dict[int, float] | None = None,
                factor_dp: float = 1.0, DN_fijo: int | None = None) -> Tubo:
    """Evalúa todos los DN y elige el más chico que respeta el Δp máximo.

    Por qué el más chico: en laminar el volumen no depende del DN, así que el
    DN solo reparte ese volumen entre largo y Δp (Δp ∝ 1/D⁶ en newtoniano:
    D⁴ por la velocidad y D² más por el largo). El Δp máximo pone la cota
    inferior; un DN más grande solo agranda el caudal de CIP necesario.

    factor_dp: el Δp calculado se multiplica por esto antes de compararlo con
    el máximo. En laminar Δp es proporcional a K, y K de la yema es incierto:
    sin factor, un DN que entra raspando queda afuera con un lote más espeso.

    DN_fijo: se usa ese DN en vez del más chico (por ejemplo, para acortar un
    tubo en horquillas, que sube 2 % de su largo); se avisa si no cumple el Δp.
    """
    DN_ext = DN_ext or {}  # si no vino, diccionario vacío
    opciones = [evaluar_DN(DN, D, productos, crit, v_CIP, geo, DN_ext.get(DN))
                for DN, D in sorted(DN_diam.items())]
    for op in opciones:
        if op.dp_max * factor_dp > dp_max_Pa:
            op.factible = False
            op.motivo = (f"Δp {op.dp_max/1e5:.2f} bar × {factor_dp} = "
                         f"{op.dp_max*factor_dp/1e5:.2f} bar > {dp_max_Pa/1e5:.2f} bar")

    avisos = []
    factibles = [op for op in opciones if op.factible]
    if DN_fijo is not None:
        if DN_fijo not in DN_diam:
            raise ValueError(f"DN_fijo = {DN_fijo} no está entre los DN disponibles {sorted(DN_diam)}")
        elegido = next(op for op in opciones if op.DN == DN_fijo)
        if not elegido.factible:
            avisos.append(f"Tubo: DN{DN_fijo} fijado en la config no cumple el Δp máximo ({elegido.motivo})")
        elif factibles and factibles[0].DN != DN_fijo:
            avisos.append(f"Tubo: DN{DN_fijo} fijado en la config; el más chico que cumple el Δp "
                          f"sería DN{factibles[0].DN}")
    elif factibles:
        elegido = factibles[0]  # están ordenados de menor a mayor DN
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

    recto = None
    if geo.es_horquillas and elegido.horquillas:
        h = elegido.horquillas
        avisos.append(
            "Tubo en horquillas: la partícula más rápida se toma como en tubo recto "
            "(las curvas mezclan algo y no se cuenta: lado seguro). Δp de curvas por "
            "Hooper (1981), pensado para newtonianos, con el Re de Metzner–Reed: aproximación"
        )
        if h.niple > 1e-6:
            avisos.append(
                f"Tubo: las curvas de 180° se arman con dos codos de 90° (R = {geo.radio_curva*1000:.0f} mm) "
                f"y un niple de {h.niple*1000:.0f} mm, para que los tramos aislados de los dos "
                f"planos no se toquen ({h.entre_planos*1000:.0f} mm entre ejes)"
            )
        if h.pendiente > geo.pendiente + 1e-9:
            avisos.append(
                f"Tubo: con tramos de {h.largo_tramo:.2f} m la pendiente sube a {h.pendiente:.1%} "
                f"para que entren los tubos aislados uno sobre otro (mínimo {geo.pendiente:.0%})"
            )
    if geo.es_espiral and elegido.D_ext:
        s_esp = pendiente_espiral(elegido.D_ext, geo)
        if s_esp < geo.pendiente:
            avisos.append(
                f"Tubo: la espiral sube {s_esp:.1%} por metro, menos que el {geo.pendiente:.0%} "
                "continuo que pide CFIA para el tubo de retención"
            )
    if geo.es_espiral:
        relacion = geo.D_espiral / elegido.D
        if not 15 <= relacion <= 860:
            avisos.append(
                f"Tubo: Dc/D = {relacion:.0f}, fuera del rango de las correlaciones de Ito (15–860)"
            )
        for clave, r in elegido.productos.items():
            if r.flujo.regimen != "turbulento" and r.flujo.Dean <= geo.Dean_min:
                avisos.append(
                    f"{clave}: Dean = {r.flujo.Dean:.1f} ≤ {geo.Dean_min}: sin flujo secundario "
                    "apreciable, se usa el factor de tubo recto"
                )
        avisos.append(
            "Tubo en espiral: corrección de fricción por curvatura (White 1929 en laminar, "
            "Ito 1959 en turbulento) pensada para newtonianos, aplicada con el Re de "
            "Metzner–Reed: aproximación"
        )
        recto = evaluar_DN(elegido.DN, elegido.D, productos, crit, v_CIP, RECTO, elegido.D_ext)

    return Tubo(elegido=elegido, opciones=opciones, geometria=geo,
                equivalente_recto=recto, avisos=avisos)
