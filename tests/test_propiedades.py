"""Tests de propiedades: contra valores conocidos (agua pura) y límites físicos.

Cómo funciona pytest: busca archivos test_*.py, corre cada función test_*, y
cada `assert` es una afirmación que tiene que ser verdadera. Si una falla,
pytest dice cuál, en qué línea y con qué valores. pytest.approx compara con
tolerancia, porque en números con coma "igual" casi nunca es exacto.
"""

import pytest

from pasteurizador import propiedades as pr


def test_choi_okos_agua_pura_a_20C():
    t = pr.choi_okos({"agua": 1.0}, 20.0)
    # el polinomio de Choi–Okos para agua da 995,7: ~0,25 % por debajo del real (998,2)
    assert t.rho == pytest.approx(998.2, rel=0.003)
    assert t.cp == pytest.approx(4182, rel=0.005)
    assert t.k == pytest.approx(0.60, abs=0.01)


def test_choi_okos_coincide_con_iapws_para_agua():
    for T in (10.0, 40.0, 60.0):
        co, iapws = pr.choi_okos({"agua": 1.0}, T), pr.agua(T)
        assert co.cp == pytest.approx(iapws.cp, rel=0.005)
        assert co.rho == pytest.approx(iapws.rho, rel=0.003)


def test_composicion_que_no_suma_1_falla():
    with pytest.raises(ValueError):
        pr.choi_okos({"agua": 0.5, "proteina": 0.2}, 20.0)


def test_mas_solidos_baja_cp():
    diluido = pr.choi_okos({"agua": 0.9, "proteina": 0.1}, 20.0)
    concentrado = pr.choi_okos({"agua": 0.7, "proteina": 0.3}, 20.0)
    assert concentrado.cp < diluido.cp
    assert concentrado.rho > diluido.rho


def test_yema_gut_reproduce_tabla_en_sus_puntos():
    # Tablas 4 y 5 de Gut et al. (2003)
    for T, K, n in [(8.9, 1.24, 0.843), (37.8, 0.20, 0.859), (60.8, 0.10, 0.874)]:
        r = pr.reologia("gut2003_yema", T)
        assert r.K == pytest.approx(K, rel=1e-9)
        assert r.n == pytest.approx(n, rel=1e-9)
        assert not r.avisos


def test_yema_gut_avisa_si_extrapola():
    assert pr.reologia("gut2003_yema", 61.4).avisos


def test_clara_vale_9_3_mPas_a_30C_y_es_newtoniana():
    r = pr.reologia("clara_uf_thapon_scalzo", 30.0)
    assert r.K == pytest.approx(9.3e-3, rel=1e-9)
    assert r.n == 1.0


def test_viscosidades_bajan_con_la_temperatura():
    for modelo in ("leon_munguia_huevo_entero", "gut2003_yema", "clara_uf_thapon_scalzo"):
        assert pr.reologia(modelo, 50.0).K < pr.reologia(modelo, 20.0).K


def test_modelo_agua_para_comparar():
    # el producto "como agua": térmicas IAPWS y viscosidad newtoniana del agua
    t = pr.termicas("agua", None, 40.0)
    assert t.cp == pytest.approx(pr.agua(40.0).cp)
    r = pr.reologia("agua", 40.0)
    assert r.n == 1.0 and r.K == pytest.approx(pr.agua_reologia(40.0).K)
