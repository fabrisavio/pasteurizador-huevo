"""Tests del tubo de retención: casos límite con solución conocida
(newtoniano, Δp ∝ 1/D⁶, volumen laminar independiente del DN) y el
comportamiento de la espiral.
"""

import math

import pytest

from pasteurizador.cinetica import Cinetica
from pasteurizador.holding import (Criterios, Geometria, K_curva, ProductoTubo, Re_critico_espiral,
                                   correccion_white_laminar, disponer_horquillas, elegir_tubo,
                                   evaluar_DN, flujo_en_tubo, pendiente_espiral,
                                   viscosidad_metzner_reed)

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


def test_la_espiral_sube_la_friccion_y_acorta_el_largo():
    p = producto(n=0.85, K=0.1)
    recto = evaluar_DN(40, 0.038, [p], CRIT, 1.5)
    enrollado = evaluar_DN(40, 0.038, [p], CRIT, 1.5, ESPIRAL, D_ext=0.041)
    # el largo baja en la proporción de los factores: 1,71 contra (3n+1)/(n+1)
    assert enrollado.L / recto.L == pytest.approx((1 / 0.586) / ((3 * 0.85 + 1) / 1.85))
    assert enrollado.productos["X"].flujo.correccion_f > 1     # más fricción por vuelta
    assert enrollado.vueltas * math.pi * 1.0 == pytest.approx(enrollado.L)


def test_factor_espiral_es_1_sobre_0586_y_nunca_mayor_que_el_recto():
    f_newt = flujo_en_tubo(producto(n=1.0, K=0.1), 0.038, CRIT, ESPIRAL)
    assert f_newt.Dean > 11.6
    assert f_newt.factor == pytest.approx(1 / 0.586)
    assert f_newt.factor_recto == pytest.approx(2.0)
    # pseudoplástico muy marcado: su perfil recto ya es más plano que 1,71
    f_plano = flujo_en_tubo(producto(n=0.4, K=0.1), 0.038, CRIT, ESPIRAL)
    assert f_plano.factor == pytest.approx((3 * 0.4 + 1) / 1.4)


def test_sin_vortices_de_Dean_usa_el_factor_recto():
    # muy viscoso: Re y Dean chiquitos, la curvatura no hace nada
    f = flujo_en_tubo(producto(n=1.0, K=50.0), 0.038, CRIT, ESPIRAL)
    assert f.Dean <= 11.6
    assert f.factor == pytest.approx(2.0)


def test_la_curvatura_retrasa_la_turbulencia():
    assert Re_critico_espiral(0.038, 1.0) > 2100
    # un Re que en recto es transición, en espiral sigue laminar
    p = producto(n=1.0, K=3.0e-3)
    assert flujo_en_tubo(p, 0.032, CRIT).regimen == "transicion"
    assert flujo_en_tubo(p, 0.032, CRIT, ESPIRAL).regimen == "laminar"


def test_espiral_informa_el_equivalente_recto():
    p = producto(n=1.0, K=0.1)
    tubo = elegir_tubo([p], {40: 0.038}, CRIT, dp_max_Pa=1e6, v_CIP=1.5, geo=ESPIRAL,
                       DN_ext={40: 0.041})
    assert tubo.equivalente_recto.L > tubo.elegido.L
    assert elegir_tubo([p], {40: 0.038}, CRIT, 1e6, 1.5).equivalente_recto is None


def test_espiral_sin_diametro_es_error():
    with pytest.raises(ValueError):
        Geometria(tipo="espiral")


# --- horquillas ---------------------------------------------------------------
HORQ = Geometria(tipo="horquillas", largo_tramo=6.0, radio_curva=0.06, luz=0.02,
                 extra_envolvente=0.10, pendiente=0.02)


def test_horquillas_reparten_el_largo_entre_tramos_y_curvas():
    h = disponer_horquillas(114.0, 0.041, HORQ)
    assert h.largo_tramo <= 6.0 + 1e-12
    assert h.curvas == h.tramos - 1
    assert h.tramos * h.largo_tramo + h.curvas * h.largo_curva == pytest.approx(114.0)
    # con el N anterior no alcanzaría: es el mínimo de tramos
    assert (h.tramos - 1) * 6.0 + (h.tramos - 2) * h.largo_curva < 114.0


def test_horquillas_separan_los_planos_para_que_entre_la_coquilla():
    h = disponer_horquillas(114.0, 0.041, HORQ)
    assert h.entre_planos == pytest.approx(0.041 + 0.10 + 0.02)
    assert h.niple == pytest.approx(h.entre_planos - 2 * 0.06)
    assert h.largo_curva == pytest.approx(math.pi * 0.06 + h.niple)


def test_desnivel_es_la_pendiente_por_el_largo():
    h = disponer_horquillas(114.0, 0.041, HORQ)
    assert h.pendiente == pytest.approx(0.02)
    assert h.desnivel == pytest.approx(0.02 * 114.0)


def test_tramos_cortos_obligan_a_subir_mas_empinado():
    geo = Geometria(tipo="horquillas", largo_tramo=2.0, radio_curva=0.06, luz=0.02,
                    extra_envolvente=0.10, pendiente=0.02)
    h = disponer_horquillas(114.0, 0.041, geo)
    # dos tramos del mismo plano tienen que quedar separados un tubo aislado más la luz
    assert 2 * h.pendiente * (h.largo_tramo + h.largo_curva) == pytest.approx(0.041 + 0.10 + 0.02)
    assert h.pendiente > 0.02


def test_K_de_curva_de_hooper():
    D = 0.0254  # 1 pulgada: 1 + 1/D = 2
    assert K_curva(1e12, D, HORQ) == pytest.approx(0.30 * 2)
    assert K_curva(1000.0, D, HORQ) == pytest.approx(1.0 + 0.6)


def test_las_curvas_suman_dp_y_no_cambian_el_largo():
    p = producto(n=0.85, K=0.1)
    recto = evaluar_DN(40, 0.038, [p], CRIT, 1.5, Geometria(), 0.041)
    horq = evaluar_DN(40, 0.038, [p], CRIT, 1.5, HORQ, 0.041)
    assert horq.L == pytest.approx(recto.L)            # misma partícula rápida que el recto
    x = horq.productos["X"]
    assert x.dp_curvas > 0
    assert x.dp == pytest.approx(recto.productos["X"].dp + x.dp_curvas)


def test_pendiente_de_la_espiral():
    geo = Geometria(tipo="espiral", D_espiral=1.0, separacion=0.01)
    assert pendiente_espiral(0.041, geo) == pytest.approx(0.051 / math.pi)


def test_DN_fijo_usa_ese_DN_y_avisa():
    p = producto(n=0.85, K=0.1)
    DNs = {40: 0.038, 50: 0.050}
    t = elegir_tubo([p], DNs, CRIT, 1e5, 1.5, DN_fijo=50)
    assert t.elegido.DN == 50
    assert any("DN50 fijado" in a for a in t.avisos)


def test_radio_de_curva_por_DN():
    geo = Geometria(tipo="horquillas", radio_curva=0.06, radios_curva=((40, 0.06), (50, 0.07)))
    p = producto(n=0.85, K=0.1)
    assert evaluar_DN(50, 0.050, [p], CRIT, 1.5, geo, 0.053).horquillas.radio == pytest.approx(0.07)
    assert evaluar_DN(40, 0.038, [p], CRIT, 1.5, geo, 0.041).horquillas.radio == pytest.approx(0.06)
    assert geo.radio_para(65) == pytest.approx(0.06)   # sin dato: el de la ficha
