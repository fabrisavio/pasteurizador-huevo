"""El pasteurizador completo, en dos pasadas.

1. Diseño por producto: qué área necesitaría cada uno por separado.
2. Equipo único: en cada sección manda el producto que más área pide.
3. Verificación: cómo funciona cada producto en ese equipo real.

Esquema por producto (sin enfriamiento: sale al secador):

    crudo ─► [regeneración, lado frío] ─► [calentamiento] ─► tubo de retención ─┐
                                              ▲ agua caliente                   │
    al secador ◄── [regeneración, lado caliente] ◄─────────────────────────────┘

Cómo se controla en planta, y cómo lo imita el modelo:
- El caudal de producto es fijo: lo impone la bomba de tiempo (timing pump).
  Es la garantía legal del tiempo de retención, así que no se toca para
  controlar la temperatura.
- El caudal de agua caliente también es fijo (bomba del circuito de agua).
- Lo que se regula es la T del agua caliente (vapor al circuito de agua),
  para que el producto llegue a la T de holding. Si a la salida del tubo
  no llega a la T legal, la válvula de desvío (FDV) lo devuelve al tanque de
  crudo [gut2003yolk], [cunningham1995pasteurization].
- El control no apunta a la T legal justa sino un poco arriba (margen de
  control): oscila unas décimas, y la FDV vuelve a dejar pasar recién un
  cuarto de grado por encima de donde desvió [cfia2026htst]. Sin margen,
  desviaría a cada rato.
- En la verificación (paso 3) el modelo hace lo mismo que ese controlador en
  régimen: con los caudales fijos, busca la T de agua que lleva el producto
  justo a T_operación + pérdida en el tubo.

Temperaturas en el tubo, de arriba a abajo:
    T_legal (norma) + margen de control = T_operación (sale del tubo)
    T_operación + pérdida en el tubo     = sale del calentamiento
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from scipy.optimize import brentq

from . import aislacion as ai
from . import propiedades as pr
from .config import Config, Producto
from .holding import ProductoTubo, Tubo, elegir_tubo
from .seccion_phe import (Corriente, Seccion, disenar_seccion,
                          temperatura_por_balance, verificar_seccion)

# Guardia, no dato de diseño. En la regeneración, en la punta fría se enfrentan
# el pasteurizado que sale hacia el secador (T_salida) y el crudo que entra
# (T_entrada). Si esa diferencia tiende a cero, el área tiende a infinito; si es
# negativa, hay cruce de temperaturas y no existe intercambiador que lo haga.
# Debajo de este valor no se regenera. En la config base no actúa: las
# diferencias son 35 °C (HP), 38 °C (YP) y 10 °C (CP).
ACERCAMIENTO_MIN_REGEN_C = 0.5


def _corriente_producto(p: Producto) -> Corriente:
    # p.termicas es una función (T → propiedades); la Corriente la guarda y la
    # llama cada vez que necesita cp a otra temperatura.
    return Corriente(p.clave, p.m, p.termicas)


# -----------------------------------------------------------------------------
# 1. Diseño de un producto
# -----------------------------------------------------------------------------
@dataclass
class DisenoProducto:
    clave: str
    regen: Seccion | None  # None si el producto no se regenera
    calent: Seccion
    m_agua: float          # kg/s de agua caliente
    T_regen_salida: float  # °C, producto crudo a la salida de la regeneración


def disenar_producto(p: Producto, cfg: Config, perdida_tubo: float) -> DisenoProducto:
    """Las áreas de regeneración y calentamiento que pediría este producto solo."""
    prod = _corriente_producto(p)
    T_op = cfg.T_operacion(p)  # a esto sale del tubo
    # El sensor legal está a la salida del tubo: para que ahí haya T_op, el
    # calentamiento tiene que entregar T_op + lo que se pierde en el tubo.
    T_calent_salida = T_op + perdida_tubo

    # Regeneración: el pasteurizado (caliente) baja de T_op a T_salida
    # y le entrega ese calor al crudo (frío), que sube por balance.
    # Es el mismo producto de los dos lados: con cp parecidos, el crudo sube
    # casi lo mismo que baja el pasteurizado.
    if p.T_salida - p.T_entrada > ACERCAMIENTO_MIN_REGEN_C:
        Q_regen = prod.C(T_op, p.T_salida) * (T_op - p.T_salida)
        T_r = temperatura_por_balance(prod, p.T_entrada, Q_regen)
        regen = disenar_seccion("regeneración", prod, prod, T_op, p.T_salida,
                                p.T_entrada, T_r, cfg.U_regen)
    else:
        regen, T_r = None, p.T_entrada

    # Calentamiento: el agua entra al tope permitido (protege la pared) y baja DT.
    # Diseñar con el agua al tope da el área mínima; en operación, el
    # controlador la baja a lo que haga falta (paso 3).
    Q_cal = prod.C(T_r, T_calent_salida) * (T_calent_salida - T_r)
    T_w_in, T_w_out = p.T_max_agua, p.T_max_agua - cfg.DT_agua
    m_w = Q_cal / (pr.agua((T_w_in + T_w_out) / 2).cp * cfg.DT_agua)
    agua = Corriente("agua caliente", m_w, pr.agua)
    # cerrar el balance del agua con su cp exacto (Q negativo: el agua se enfría)
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
    perdida_tubo: dict[str, float]   # °C que pierde cada producto en el tubo
    # Las formas de aislar evaluadas (tubo / bloque); vacío si la pérdida es fija
    aislacion: dict[str, ai.OpcionAislacion] = field(default_factory=dict)
    # Placas de calentamiento que haría falta con cada forma de aislar
    placas_calent_por_aislacion: dict[str, int] = field(default_factory=dict)
    avisos: list[str] = field(default_factory=list)


def _placas(A: float, A_placa: float) -> int:
    """Placas térmicas necesarias (cada una aporta A_placa de área)."""
    # ceil redondea para arriba: 11,2 placas → 12. El −1e-9 evita que un
    # 12,0000000001 por redondeo numérico se convierta en 13.
    return max(1, math.ceil(A / A_placa - 1e-9))


def T_media_tubo(p: Producto, cfg: Config, perdida_tubo: float) -> float:
    """Entra a T_op + pérdida, sale a T_op: la media está a mitad de camino."""
    return cfg.T_operacion(p) + perdida_tubo / 2


def productos_tubo(cfg: Config, perdidas: dict[str, float]) -> list[ProductoTubo]:
    """Cada producto con sus propiedades a la T media del tubo."""
    lista = []
    for p in cfg.productos.values():
        T_media = T_media_tubo(p, cfg, perdidas[p.clave])
        reo = p.reologia(T_media)
        lista.append(ProductoTubo(
            clave=p.clave, m=p.m, rho=p.termicas(T_media).rho, K=reo.K, n=reo.n,
            # los log se cuentan a la T legal, no a la de operación: lado seguro
            T_holding=p.T_holding, t_holding=p.t_holding, cinetica=p.cinetica,
        ))
    return lista


def _perdidas_por_aislacion(cfg: Config, tubo: Tubo, perdidas: dict[str, float]):
    """Pérdida de cada producto con las dos formas de aislar el tubo elegido."""
    t, geo = tubo.elegido, tubo.geometria
    termicos = [ai.ProductoTermico(
        clave=p.clave, m=p.m, cp=p.termicas(T_media_tubo(p, cfg, perdidas[p.clave])).cp,
        T_sale=cfg.T_operacion(p), horas_anio=p.horas_anio,
    ) for p in cfg.productos.values()]
    D_ext = t.D_ext or t.D  # sin dato de diámetro exterior, el interior (aproximación)
    return ai.opciones(t.L, D_ext, geo.D_espiral if geo.es_espiral else None, t.alto,
                       termicos, cfg.aislacion)


def disenar_equipo(cfg: Config) -> Equipo:
    """Un solo equipo para los tres: en cada sección, el área más grande que se pidió."""
    avisos: list[str] = []
    # La pérdida en el tubo depende del largo, y el largo (vía propiedades a la
    # T media) depende un poquito de la pérdida. Se resuelve por punto fijo:
    # se arranca con 0,5 °C, se elige el tubo, se calcula la pérdida, y se
    # repite hasta que no cambie (converge en 2 o 3 vueltas).
    semilla = cfg.perdida_fija if cfg.perdida_fija is not None else 0.5
    perdidas = {k: semilla for k in cfg.productos}
    opciones_ais: dict[str, ai.OpcionAislacion] = {}
    for _ in range(20):
        tubo = elegir_tubo(productos_tubo(cfg, perdidas), cfg.DN_diam, cfg.criterios,
                           cfg.dp_max_holding, cfg.v_CIP, cfg.geometria_tubo, cfg.DN_ext,
                           cfg.factor_dp)
        if cfg.aislacion is None:
            break  # pérdida fija: no hay nada que iterar
        opciones_ais = _perdidas_por_aislacion(cfg, tubo, perdidas)
        forma = cfg.aislacion.forma
        if forma not in opciones_ais:  # bloque pedido para un tubo recto
            avisos.append("Aislación en bloque solo tiene sentido con el tubo en espiral; "
                          "se usa coquilla a lo largo del tubo")
            forma = "tubo"
        nuevas = opciones_ais[forma].dT
        cambio = max(abs(nuevas[k] - perdidas[k]) for k in perdidas)
        perdidas = nuevas
        if cambio < 1e-4:
            break
    else:
        # el else de un for corre si el for terminó sin break
        raise RuntimeError("La pérdida del tubo no converge")
    if cfg.aislacion is not None and cfg.aislacion.forma == "bloque" and "bloque" in opciones_ais:
        avisos.append(
            "Aislación en bloque: el aire dentro de la camisa se supone a la T del producto "
            "(cota superior de la pérdida). La camisa tiene que ser desmontable para "
            "inspeccionar el tubo (CFIA)"
        )

    disenos = {k: disenar_producto(p, cfg, perdidas[k]) for k, p in cfg.productos.items()}

    # Para comparar: con la otra forma de aislar, el producto tendría que
    # salir más caliente del calentamiento, y cerca del tope del agua cada
    # décima de grado cuesta área (el acercamiento en la punta caliente se
    # achica y el NTU se dispara).
    placas_por_aislacion = {}
    for forma, op in opciones_ais.items():
        A = max(disenar_producto(p, cfg, op.dT[k]).calent.A for k, p in cfg.productos.items())
        placas_por_aislacion[forma] = _placas(A, cfg.area_placa)

    A_reg = {k: (d.regen.A if d.regen else 0.0) for k, d in disenos.items()}
    A_cal = {k: d.calent.A for k, d in disenos.items()}
    gob_reg = max(A_reg, key=A_reg.get)
    gob_cal = max(A_cal, key=A_cal.get)
    n_reg = _placas(A_reg[gob_reg], cfg.area_placa)
    n_cal = _placas(A_cal[gob_cal], cfg.area_placa)

    return Equipo(
        placas_regen=n_reg, placas_calent=n_cal,
        A_regen=n_reg * cfg.area_placa, A_calent=n_cal * cfg.area_placa,
        gobierna_regen=gob_reg, gobierna_calent=gob_cal,
        # la bomba se elige para el que más agua pide; los demás usan el mismo caudal
        m_agua=max(d.m_agua for d in disenos.values()),
        tubo=tubo, disenos=disenos, perdida_tubo=perdidas,
        aislacion=opciones_ais, placas_calent_por_aislacion=placas_por_aislacion,
        avisos=list(dict.fromkeys(avisos)),
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

    @property
    def T_regen_salida(self) -> float:
        """°C, crudo a la salida de la regeneración (entra al calentamiento)."""
        return self.regen.Tf_out

    @property
    def T_calent_salida(self) -> float:
        """°C, producto a la salida del calentamiento (entra al tubo)."""
        return self.calent.Tf_out if self.calent else float("nan")


def verificar_producto(p: Producto, eq: Equipo, cfg: Config,
                       T_entrada: float | None = None) -> Verificacion:
    """Cómo opera el producto p en el equipo eq, entrando a T_entrada."""
    T_in = p.T_entrada if T_entrada is None else T_entrada
    prod = _corriente_producto(p)
    T_op = cfg.T_operacion(p)
    T_calent_salida = T_op + eq.perdida_tubo[p.clave]
    alertas: list[str] = []
    avisos: list[str] = []

    # Regeneración: la T de salida del tubo la fija el control, así que el lazo
    # se abre: entra caliente a T_op y frío a T_entrada, con el área instalada.
    # Como el área es la del producto que más pide, los demás regeneran de más.
    regen = verificar_seccion("regeneración", prod, prod, T_op, T_in,
                              eq.A_regen, cfg.U_regen)
    T_r, T_salida = regen.Tf_out, regen.Tc_out

    # Calentamiento: con el caudal de agua fijo, ¿a qué T tiene que entrar el agua
    # para que el producto llegue a T_op + pérdida? Es lo que hace el
    # controlador de temperatura en régimen.
    agua = Corriente("agua caliente", eq.m_agua, pr.agua)

    def falta(T_w: float) -> float:
        # Función interna: cuánto le falta (o sobra) al producto con el agua a T_w.
        # Cero = justo. brentq busca ese cero achicando un intervalo que lo
        # encierra, como una bisección pero más rápida.
        s = verificar_seccion("calentamiento", agua, prod, T_w, T_r, eq.A_calent, cfg.U_calent)
        return s.Tf_out - T_calent_salida

    T_w_max = 99.0  # techo físico del agua sin presurizar
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

    # Más frío que lo pedido se acepta, pero más frío = más viscoso: se avisa
    # si la consistencia sube más de 20 % respecto de la de diseño del secador.
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
        E_calor_MWh_anio=Q_cal * p.horas_anio / 1000,  # kW·h/año → MWh/año
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
        return not self.alertas  # lista vacía = no hay alertas = cumple

    def en_diseno(self, clave: str) -> Verificacion:
        """La verificación del producto a su T de entrada de diseño."""
        T = self.cfg.productos[clave].T_entrada
        # next(...) devuelve el primero de la lista que cumple la condición
        return next(v for v in self.verificaciones[clave] if v.T_entrada == T)


def correr(cfg: Config) -> Resultado:
    """Diseña el equipo y verifica cada producto en todo su rango de entrada."""
    avisos = [
        "v0: U fijo por sección tomado del caso de yema de Gut et al. (2003); "
        "los tres productos usan el mismo U hasta que v1 meta correlaciones",
    ]
    eq = disenar_equipo(cfg)
    avisos += eq.tubo.avisos + eq.avisos
    alertas: list[str] = []

    verifs: dict[str, list[Verificacion]] = {}
    for k, p in cfg.productos.items():
        # {a, *lista} arma un conjunto: elimina repetidos (p. ej. CP: 35 y 35)
        temps = sorted({p.T_entrada, *p.T_entrada_rango})
        verifs[k] = [verificar_producto(p, eq, cfg, T) for T in temps]
        for v in verifs[k]:
            alertas += v.alertas
            if v.T_entrada == p.T_entrada:  # lo informativo, solo en el punto de diseño
                avisos += v.avisos

    # avisos de propiedades: reología en el tubo y a la salida
    for p in cfg.productos.values():
        avisos += p.reologia(T_media_tubo(p, cfg, eq.perdida_tubo[p.clave])).avisos
        # El producto más caliente está a la entrada del tubo. Si pasa la T donde
        # empieza a espesar, las reologías medidas (que no llegan ahí) no lo ven.
        T_max = cfg.T_operacion(p) + eq.perdida_tubo[p.clave]
        if p.T_espesa is not None and T_max > p.T_espesa:
            avisos.append(
                f"{p.clave}: en el tubo llega a {T_max:.1f} °C y empieza a espesar desde "
                f"{p.T_espesa:.1f} °C; la reología del modelo no lo ve: el Δp real del tubo "
                "puede ser mayor"
            )
        r = eq.tubo.elegido.productos[p.clave]
        if r.logs_min < cfg.logs_objetivo:
            avisos.append(
                f"{p.clave}: la partícula más rápida logra {r.logs_min:.1f} log en el tubo "
                f"(objetivo {cfg.logs_objetivo:.0f}); el diseño se apoya en el par T–t "
                "normativo y la letalidad de las placas entra en v2"
            )
    # dict.fromkeys elimina repetidos sin cambiar el orden
    return Resultado(cfg, eq, verifs, list(dict.fromkeys(alertas)), list(dict.fromkeys(avisos)))
