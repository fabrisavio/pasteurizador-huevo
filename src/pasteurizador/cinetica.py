"""Cinética de muerte térmica de primer orden (modelo D–z).

Verifica, no dimensiona: el diseño se basa en el par T–t de la norma (safe
harbor de FSIS o RBN). Esto solo informa cuántos logaritmos se logran.
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
        """Reducciones decimales logradas a T constante durante t minutos."""
        return t_min / self.D(T)
