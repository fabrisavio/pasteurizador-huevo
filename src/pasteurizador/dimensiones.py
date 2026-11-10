"""Dimensiones para el plano: paquete de placas, bastidor y tubo de retención.

Qué resuelve: las cotas que necesita el plano de planta. El modelo sabe
cuántas placas hay y cuánto mide el tubo; el bastidor depende del proveedor,
así que se toma de un equipo comercial con la misma placa (canal libre de
3 mm, conexiones DN50) y se marca como referencia, no como dato.

El paquete de placas apretado mide (placas) × (paso), con paso = canal + chapa:
cada placa ocupa el espesor de su chapa más el canal que deja hasta la
siguiente. Entre regeneración y calentamiento va un bastidor intermedio
(conector) que lleva las conexiones de las dos secciones.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Bastidor:
    """Bastidor de referencia (Alfa Laval M6) y supuestos del paquete."""
    alto: float = 0.92                  # m
    ancho: float = 0.32                 # m
    largo_rango: tuple[float, float] = (0.515, 1.43)  # m, según el largo de bastidor
    conector: float = 0.10              # m, bastidor intermedio entre secciones (supuesto)
    extremos_por_seccion: int = 2       # placas de principio y fin de cada sección (supuesto)


@dataclass
class Paquete:
    placas_termicas: int
    placas_total: int      # con las de extremo
    paso: float            # m por placa apretada
    largo: float           # m, placas apretadas + conector
    # Largo de bastidor estimado: el paquete, otro tanto libre para correr la
    # placa de presión y sacar placas, y las dos placas gruesas del bastidor.
    largo_bastidor: float


def paquete(placas_regen: int, placas_calent: int, b: float, espesor: float,
            bas: Bastidor, chapas_bastidor: float = 0.15) -> Paquete:
    """Largo del paquete y del bastidor que hace falta para abrirlo.

    chapas_bastidor: placa fija + placa de presión + columna, del orden de
    5 cm cada una (supuesto). El espacio para abrir se toma igual al paquete,
    que es lo que hay que correr la placa de presión para sacar cualquier
    placa sin desarmar las demás.
    """
    termicas = placas_regen + placas_calent
    total = termicas + 2 * bas.extremos_por_seccion
    paso = b + espesor
    largo = total * paso + bas.conector
    return Paquete(termicas, total, paso, largo, 2 * largo + chapas_bastidor)
