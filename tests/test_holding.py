import math

import pytest

from pasteurizador.cinetica import Cinetica
from pasteurizador.holding import (Criterios, Geometria, ProductoTubo, Re_critico_espiral,
                                   correccion_white_laminar, elegir_tubo, evaluar_DN,
                                   flujo_en_tubo, sensibilidad_factor, viscosidad_metzner_reed)

CIN = Cinetica(1.0, 60.0, 5.0)
CRIT = Criterios()


def producto(n=1.0, K=0.05, m=0.3, clave="X"):
    return ProductoTubo(clave, m=m, rho=1030.0, K=K, n=n, T_holding=60.0,
                        t_holding=210.0, cinetica=CIN)


def test_metzner_reed_con_n_1_es_la_viscosidad():
    assert viscosidad_metzner_reed(K=0.01, n=1.0, v=0.5, D=0.04) == pytest.approx(0.01)


def test_reynolds_newtoniano_de_siempre():
    p = producto(n=1.0, K=0.01)
    f = flujo_en_tubo(p, 0.038, CRIT)
    v = p.m / (p.rho * math.pi * 0.038**2 / 4)
    assert f.Re == pytest.approx(p.rho * v * 0.038 / 0.01)


def test_factor_laminar_2_si_newtoniano():
    assert flujo_en_tubo(producto(n=1.0, K=0.1), 0.05, CRIT).factor == pytest.approx(2.0)


def test_factor_laminar_ley_de_potencia():
    n = 0.85
    f = flujo_en_tubo(producto(n=n, K=0.1), 0.05, CRIT)
    assert f.regimen == "laminar"
    assert f.factor == pytest.approx((3 * n + 1) / (n + 1))


def test_transicion_usa_factor_laminar():
    # K chico para caer entre 2100 y 4000
    f = flujo_en_tubo(producto(n=1.0, K=3.0e-3), 0.032, CRIT)
    assert f.regimen == "transicion"
    assert f.factor == pytest.approx(2.0)


def test_volumen_laminar_no_depende_del_DN():
    p = producto(n=0.85, K=0.1)
    vols = [evaluar_DN(DN, D, [p], CRIT, 1.5).volumen_L for DN, D in [(40, 0.038), (50, 0.05)]]
    assert vols[0] == pytest.approx(vols[1], rel=1e-9)


def test_dp_laminar_newtoniano_va_con_D_a_la_menos_6():
    p = producto(n=1.0, K=0.1)
    a, b = evaluar_DN(40, 0.038, [p], CRIT, 1.5), evaluar_DN(50, 0.05, [p], CRIT, 1.5)
    assert a.dp_max / b.dp_max == pytest.approx((0.05 / 0.038) ** 6, rel=1e-6)


def test_partícula_mas_rapida_del_que_gobierna_dura_t_holding():
    lento, rapido = producto(n=0.85, K=0.1, clave="viscoso"), producto(n=1.0, K=0.002, clave="fluido")
    op = evaluar_DN(40, 0.038, [lento, rapido], CRIT, 1.5)
    gob = op.productos[op.gobierna]
    assert gob.t_rapida_min == pytest.approx(3.5)
    assert all(r.t_rapida_min >= 3.5 - 1e-9 for r in op.productos.values())


def test_elige_el_DN_mas_chico_que_cumple_dp():
    p = producto(n=0.85, K=0.15)
    DNs = {25: 0.026, 32: 0.032, 40: 0.038, 50: 0.05}
    tubo = elegir_tubo([p], DNs, CRIT, dp_max_Pa=0.5e5, v_CIP=1.5)
    factibles = [op.DN for op in tubo.opciones if op.factible]
    assert tubo.elegido.DN == min(factibles)
    assert tubo.elegido.dp_max <= 0.5e5


# --- tubo enrollado en espiral ------------------------------------------------
ESPIRAL = Geometria(tipo="espiral", D_espiral=1.0)


def test_white_vale_1_por_debajo_de_Dean_11_6_y_crece_despues():
    assert correccion_white_laminar(10.0) == 1.0
    assert 1.0 < correccion_white_laminar(30.0) < correccion_white_laminar(300.0)


def test_la_espiral_sube_la_friccion_pero_no_cambia_el_largo():
    p = producto(n=0.85, K=0.1)
    recto = evaluar_DN(40, 0.038, [p], CRIT, 1.5)
    enrollado = evaluar_DN(40, 0.038, [p], CRIT, 1.5, ESPIRAL, D_ext=0.041)
    assert enrollado.L == pytest.approx(recto.L)          # el diseño sigue con el factor recto
    assert enrollado.dp_max > recto.dp_max
    assert enrollado.vueltas * math.pi * 1.0 == pytest.approx(enrollado.L)


def test_la_curvatura_retrasa_la_turbulencia():
    assert Re_critico_espiral(0.038, 1.0) > 2100
    # un Re que en recto es transición, en espiral sigue laminar
    p = producto(n=1.0, K=3.0e-3)
    assert flujo_en_tubo(p, 0.032, CRIT).regimen == "transicion"
    assert flujo_en_tubo(p, 0.032, CRIT, ESPIRAL).regimen == "laminar"


def test_sensibilidad_escala_lineal_con_el_factor():
    p = producto(n=1.0, K=0.1)  # laminar newtoniano: factor de diseño 2
    op = evaluar_DN(40, 0.038, [p], CRIT, 1.5, ESPIRAL, D_ext=0.041)
    pts = sensibilidad_factor(op, [p], [2.0, 1.5, 1.0], ESPIRAL)
    assert pts[0].L == pytest.approx(op.L)
    assert pts[1].L == pytest.approx(0.75 * op.L)
    assert pts[2].dp_max == pytest.approx(0.5 * op.dp_max)


def test_espiral_sin_diametro_es_error():
    with pytest.raises(ValueError):
        Geometria(tipo="espiral")
