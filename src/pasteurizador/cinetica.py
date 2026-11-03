"""Cinética de muerte térmica de primer orden (modelo D–z).

Verifica, no dimensiona: el diseño se basa en el par T–t de la norma (safe
harbor de FSIS o RBN). Esto solo informa cuántos logaritmos se logran.

Intuición: es una reacción de primer orden donde el "reactivo" es la
población de Salmonella. D es el tiempo para que sobreviva 1 de cada 10
(una reducción decimal, "1 log"). z es la sensibilidad a la temperatura:
subir z grados divide D por 10. Es la misma idea que la energía de activación
de Arrhenius, escrita en el lenguaje de la microbiología.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Cinetica:
    D_ref_min: float  # min: tiempo para reducir 1 log a T_ref
    T_ref_C: float
    z_C: float        # °C que hay que subir para dividir D por 10

    def D(self, T: float) -> float:
        """D(T) = D_ref · 10^((T_ref − T)/z), en minutos."""
        return self.D_ref_min * 10 ** ((self.T_ref_C - T) / self.z_C)

    def logs(self, t_min: float, T: float) -> float:
        """Reducciones decimales logradas a T constante durante t minutos.

        Ejemplo: D = 1,25 min y t = 3,5 min → 2,8 log: sobrevive 1 de cada
        10^2,8 ≈ 630.
        """
        return t_min / self.D(T)
