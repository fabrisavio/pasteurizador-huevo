"""Validación contra el pasteurizador de yema de Gut et al. (2003).

Por qué este caso: es un pasteurizador de placas real para yema, publicado
con todos sus datos (caudales, áreas, U) y con las temperaturas que calcula
su modelo en cada punto. Si el código, con esos mismos datos, llega a las
mismas temperaturas, los bloques (secciones, balances, propiedades,
encadenamiento) están bien armados.

Los puntos de temperatura usan la numeración del paper (T1…T6), para que se
puedan comparar directo con su Tabla 7:

    T1 yema cruda (5 °C)
     │
     ▼
    [regeneración, lado frío] ─T2─► [calentamiento] ─T3─► tubo de retención ─T4─┐
                                        ▲ agua a 66 °C                         │
    T6 ◄─ [enfriamiento] ◄─T5─ [regeneración, lado caliente] ◄─────────────────┘
               ▲ agua a 2 °C

A diferencia de nuestra planta, el caso de Gut tiene enfriamiento: la yema
pasteurizada va a almacenamiento en frío, no a un secador.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml
from scipy.optimize import brentq

from . import placas as pl
from . import propiedades as pr
from .seccion_phe import Corriente, temperatura_por_balance, verificar_seccion

# Qué es cada punto, en palabras. Mismo orden que la Tabla 7 de Gut.
PUNTOS = {
    "T2": "Yema cruda a la salida de la regeneración (entra al calentamiento)",
    "T3": "Yema a la salida del calentamiento (entra al tubo de retención)",
    "T4": "Yema a la salida del tubo de retención (punto legal de control)",
    "T5": "Yema pasteurizada a la salida de la regeneración (entra al enfriamiento)",
    "T6": "Yema a la salida del enfriamiento (va a almacenamiento)",
}


@dataclass
class ResultadoValidacion:
    calculado: dict[str, float]
    esperado: dict[str, float]
    tolerancia: float
    T1: float  # °C, entrada de la yema (dato, no se calcula)

    @property
    def diferencias(self) -> dict[str, float]:
        return {k: self.calculado[k] - self.esperado[k] for k in self.esperado}

    @property
    def ok(self) -> bool:
        return all(abs(d) <= self.tolerancia for d in self.diferencias.values())

    def tabla(self) -> str:
        """Tabla de texto de ancho fijo, para leer en la terminal."""
        lineas = [f"Dato: la yema entra a T1 = {self.T1:.1f} °C.\n",
                  f"{'Punto':<6}{'Gut 2003':>10}{'Modelo':>9}{'Dif.':>7}   Qué es",
                  f"{'':<6}{'(°C)':>10}{'(°C)':>9}{'(°C)':>7}"]
        for k, d in self.diferencias.items():
            lineas.append(f"{k:<6}{self.esperado[k]:>10.2f}{self.calculado[k]:>9.2f}"
                          f"{d:>+7.2f}   {PUNTOS.get(k, '')}")
        estado = "PASA" if self.ok else "NO PASA"
        lineas.append(f"\nTolerancia ±{self.tolerancia} °C: {estado}")
        return "\n".join(lineas)


def validar(ruta: str | Path) -> ResultadoValidacion:
    with open(ruta, encoding="utf-8") as f:
        c = yaml.safe_load(f)

    # lambda: una función de una línea. Acá: T → propiedades de la yema según Gut.
    termica = lambda T: pr.termicas(c["producto"]["termicas"], None, T)  # noqa: E731
    T1 = c["producto"]["T_entrada_C"]
    m = c["producto"]["caudal_L_h"] / 3.6e6 * termica(T1).rho  # L/h → m³/s → kg/s, a T de entrada
    yema = Corriente("yema", m, termica)
    caliente = Corriente("agua caliente", c["agua_caliente"]["caudal_kg_s"], pr.agua)
    fria = Corriente("agua fría", c["agua_fria"]["caudal_kg_s"], pr.agua)

    placa = c["placa"]
    A_placa = placa["largo_m"] * placa["ancho_m"] * placa["factor_area"]
    # Con NC canales hay NC − 1 placas que transfieren calor entre dos canales
    A = {s: (nc - 1) * A_placa for s, nc in c["canales"].items()}
    U = c["U_W_m2K"]
    perdida = c["holding"]["perdida_T_C"]
    Tw_cal = c["agua_caliente"]["T_entrada_C"]
    Tw_fri = c["agua_fria"]["T_entrada_C"]

    def vuelta(T4: float):
        """Una vuelta al lazo suponiendo T4: regeneración y calentamiento."""
        reg = verificar_seccion("regeneración", yema, yema, T4, T1, A["regeneracion"], U["regeneracion"])
        cal = verificar_seccion("calentamiento", caliente, yema, Tw_cal, reg.Tf_out,
                                A["calentamiento"], U["calentamiento"])
        return reg, cal

    # El lazo se muerde la cola: T4 calienta al crudo en la regeneración, y de
    # cuánto se calienta el crudo depende T3 y por lo tanto T4. Se busca el T4
    # que es consistente consigo mismo: T4 = T3(T4) − pérdida en el tubo.
    T4 = brentq(lambda T: vuelta(T)[1].Tf_out - perdida - T, T1 + 1, Tw_cal - 1e-3)
    reg, cal = vuelta(T4)
    enf = verificar_seccion("enfriamiento", yema, fria, reg.Tc_out, Tw_fri,
                            A["enfriamiento"], U["enfriamiento"])
    calculado = {"T2": reg.Tf_out, "T3": cal.Tf_out, "T4": T4, "T5": reg.Tc_out, "T6": enf.Tc_out}
    return ResultadoValidacion(calculado, c["esperado_C"], c["tolerancia_C"], T1)


# -----------------------------------------------------------------------------
# Chequeo del U de Gut (v1)
# -----------------------------------------------------------------------------
@dataclass
class ChequeoU:
    seccion: str
    U_publicado: float
    U_calculado: float
    h_caliente: float
    h_fria: float

    @property
    def relacion(self) -> float:
        return self.U_publicado / self.U_calculado


def chequear_U_gut(ruta: str | Path) -> list[ChequeoU]:
    """U de cada sección con la correlación y los datos que publica Gut.

    No es una validación que tenga que pasar: muestra si el U que v0 tomó de
    Gut (Tabla 7) se puede reconstruir con lo que dice el propio paper.
    """
    with open(ruta, encoding="utf-8") as f:
        c = yaml.safe_load(f)
    ch = c["chequeo_U"]
    placa = pl.Placa(area=c["placa"]["largo_m"] * c["placa"]["ancho_m"], ancho=c["placa"]["ancho_m"],
                     b=ch["placa"]["b_m"], angulo=45.0, phi=c["placa"]["factor_area"],
                     espesor=ch["placa"]["espesor_m"], k=ch["placa"]["k_W_mK"], D_puerto=0.05)
    yema = pl.Fluido("yema", lambda T: pr.termicas(c["producto"]["termicas"], None, T),
                     lambda T: pr.reologia("gut2003_yema", T))
    T1 = c["producto"]["T_entrada_C"]
    m = c["producto"]["caudal_L_h"] / 3.6e6 * yema.termicas(T1).rho
    T = ch["T_constant_U_C"]
    mw_c, mw_f = c["agua_caliente"]["caudal_kg_s"], c["agua_fria"]["caudal_kg_s"]
    Tw_c, Tw_f = c["agua_caliente"]["T_entrada_C"], c["agua_fria"]["T_entrada_C"]

    # Salida de cada agua por balance con la yema
    y = Corriente("yema", m, yema.termicas)
    Q_cal = y.C(T["T2"], T["T3"]) * (T["T3"] - T["T2"])
    Q_enf = y.C(T["T5"], T["T6"]) * (T["T5"] - T["T6"])
    Tw_c_sale = temperatura_por_balance(Corriente("agua", mw_c, pr.agua), Tw_c, -Q_cal)
    Tw_f_sale = temperatura_por_balance(Corriente("agua", mw_f, pr.agua), Tw_f, Q_enf)

    def arr(sec, lado):
        return pl.Arreglo(*ch["arreglos"][sec][lado])

    casos = {
        "regeneracion": (pl.Lado(yema, m, arr("regeneracion", "caliente")), (T["T4"] + T["T5"]) / 2,
                         pl.Lado(yema, m, arr("regeneracion", "fria")), (T1 + T["T2"]) / 2),
        "calentamiento": (pl.Lado(pl.AGUA, mw_c, arr("calentamiento", "caliente")), (Tw_c + Tw_c_sale) / 2,
                          pl.Lado(yema, m, arr("calentamiento", "fria")), (T["T2"] + T["T3"]) / 2),
        "enfriamiento": (pl.Lado(yema, m, arr("enfriamiento", "caliente")), (T["T5"] + T["T6"]) / 2,
                         pl.Lado(pl.AGUA, mw_f, arr("enfriamiento", "fria")), (Tw_f + Tw_f_sale) / 2),
    }
    res = []
    for sec, (cal, Tc, fri, Tf) in casos.items():
        det = pl.coeficiente_U(placa, cal, fri, Tc, Tf, ch["ensuciamiento_m2K_W"][sec], "sharifi")
        res.append(ChequeoU(sec, ch["U_publicado_W_m2K"][sec], det.U, det.caliente.h, det.fria.h))
    return res


def tabla_chequeo_U(filas: list[ChequeoU]) -> str:
    lineas = [f"{'Sección':<15}{'U Gut':>8}{'U calculado':>13}{'Gut/calc.':>11}{'h cal.':>9}{'h fría':>9}",
              f"{'':<15}{'(W/m²K)':>8}{'(W/m²K)':>13}{'':>11}{'(W/m²K)':>9}{'(W/m²K)':>9}"]
    for f in filas:
        lineas.append(f"{f.seccion:<15}{f.U_publicado:>8.0f}{f.U_calculado:>13.0f}"
                      f"{f.relacion:>11.1f}{f.h_caliente:>9.0f}{f.h_fria:>9.0f}")
    return "\n".join(lineas)
