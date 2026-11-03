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

from . import propiedades as pr
from .seccion_phe import Corriente, verificar_seccion

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
