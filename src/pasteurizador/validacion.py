"""Validación contra el pasteurizador de yema de Gut et al. (2003).

Encadena regeneración → calentamiento → tubo → regeneración → enfriamiento con
los mismos bloques que usa el diseño, y compara T2…T6 con la Tabla 7.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml
from scipy.optimize import brentq

from . import propiedades as pr
from .seccion_phe import Corriente, verificar_seccion


@dataclass
class ResultadoValidacion:
    calculado: dict[str, float]
    esperado: dict[str, float]
    tolerancia: float

    @property
    def diferencias(self) -> dict[str, float]:
        return {k: self.calculado[k] - self.esperado[k] for k in self.esperado}

    @property
    def ok(self) -> bool:
        return all(abs(d) <= self.tolerancia for d in self.diferencias.values())

    def tabla(self) -> str:
        lineas = ["| T | Gut 2003 (°C) | Modelo (°C) | Diferencia (°C) |", "|---|---|---|---|"]
        for k, d in self.diferencias.items():
            lineas.append(f"| {k} | {self.esperado[k]:.2f} | {self.calculado[k]:.2f} | {d:+.2f} |")
        estado = "PASA" if self.ok else "NO PASA"
        lineas.append(f"\nTolerancia ±{self.tolerancia} °C: **{estado}**")
        return "\n".join(lineas)


def validar(ruta: str | Path) -> ResultadoValidacion:
    with open(ruta, encoding="utf-8") as f:
        c = yaml.safe_load(f)

    termica = lambda T: pr.termicas(c["producto"]["termicas"], None, T)  # noqa: E731
    T1 = c["producto"]["T_entrada_C"]
    m = c["producto"]["caudal_L_h"] / 3.6e6 * termica(T1).rho  # caudal medido a T de entrada
    yema = Corriente("yema", m, termica)
    caliente = Corriente("agua caliente", c["agua_caliente"]["caudal_kg_s"], pr.agua)
    fria = Corriente("agua fría", c["agua_fria"]["caudal_kg_s"], pr.agua)

    placa = c["placa"]
    A_placa = placa["largo_m"] * placa["ancho_m"] * placa["factor_area"]
    A = {s: (nc - 1) * A_placa for s, nc in c["canales"].items()}
    U = c["U_W_m2K"]
    perdida = c["holding"]["perdida_T_C"]
    Tw_cal = c["agua_caliente"]["T_entrada_C"]
    Tw_fri = c["agua_fria"]["T_entrada_C"]

    def vuelta(T4: float):
        reg = verificar_seccion("regeneración", yema, yema, T4, T1, A["regeneracion"], U["regeneracion"])
        cal = verificar_seccion("calentamiento", caliente, yema, Tw_cal, reg.Tf_out,
                                A["calentamiento"], U["calentamiento"])
        return reg, cal

    # La T de salida del tubo (T4) depende de sí misma a través de la regeneración:
    # se busca el punto fijo T4 = T3(T4) − pérdida.
    T4 = brentq(lambda T: vuelta(T)[1].Tf_out - perdida - T, T1 + 1, Tw_cal - 1e-3)
    reg, cal = vuelta(T4)
    enf = verificar_seccion("enfriamiento", yema, fria, reg.Tc_out, Tw_fri,
                            A["enfriamiento"], U["enfriamiento"])
    calculado = {"T2": reg.Tf_out, "T3": cal.Tf_out, "T4": T4, "T5": reg.Tc_out, "T6": enf.Tc_out}
    return ResultadoValidacion(calculado, c["esperado_C"], c["tolerancia_C"])
