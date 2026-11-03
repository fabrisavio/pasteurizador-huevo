"""Tests de la pérdida de calor del tubo: fórmulas contra cuentas a mano y
sentido físico (más aislante, menos pérdida; menos superficie, menos pérdida).
"""

import math

import pytest

from pasteurizador import aislacion as ai

A = ai.Aislacion(forma="tubo", espesor=0.05, k=0.04, h_ext=10.0, T_amb=20.0)


def producto(m=0.3, cp=3500.0, T=60.5, horas=1000.0, clave="X"):
    return ai.ProductoTermico(clave, m, cp, T, horas)


def test_ua_tubo_contra_cuenta_a_mano():
    # 1 m de tubo de 41 mm de diámetro exterior con 50 mm de lana
    UA, area, vol = ai.ua_tubo(1.0, 0.041, A)
    R = math.log(0.141 / 0.041) / (2 * math.pi * 0.04) + 1 / (10 * math.pi * 0.141)
    assert UA == pytest.approx(1 / R)
    assert area == pytest.approx(math.pi * 0.141)
    assert vol == pytest.approx(math.pi / 4 * (0.141**2 - 0.041**2))


def test_perdida_chica_es_el_balance_lineal():
    # con UA ≪ m·cp la exponencial se vuelve UA·ΔT/(m·cp)
    p = producto()
    UA = 1.0
    assert ai.perdida(UA, p, 20.0) == pytest.approx(UA * (60.5 - 20.0) / (0.3 * 3500), rel=1e-3)


def test_mas_aislante_menos_perdida():
    p = producto()
    dT = [ai.perdida(ai.ua_tubo(90.0, 0.041, ai.Aislacion("tubo", e, 0.04))[0], p, 20.0)
          for e in (0.0, 0.025, 0.05, 0.1)]
    assert dT == sorted(dT, reverse=True)


def test_bloque_pierde_menos_que_coquilla_en_una_espiral_compacta():
    ops = ai.opciones(94.0, 0.041, 1.0, 1.53, [producto()], A)
    assert ops["bloque"].dT["X"] < ops["tubo"].dT["X"]
    assert ops["bloque"].area_ext < ops["tubo"].area_ext


def test_tubo_recto_solo_tiene_coquilla():
    ops = ai.opciones(94.0, 0.041, None, None, [producto()], A)
    assert set(ops) == {"tubo"}


def test_energia_anual_es_potencia_por_horas():
    p = producto(horas=2000.0)
    op = ai.opciones(94.0, 0.041, 1.0, 1.53, [p], A)["tubo"]
    assert op.E_MWh_anio == pytest.approx(op.Q["X"] * 2000.0 / 1e6)


def test_forma_desconocida_es_error():
    with pytest.raises(ValueError):
        ai.Aislacion(forma="manta")
