import pytest

from pasteurizador import propiedades as pr
from pasteurizador.seccion_phe import (Corriente, disenar_seccion, eps_contracorriente,
                                       ntu_contracorriente, temperatura_por_balance,
                                       verificar_seccion)


def test_eps_tiende_a_1_con_area_infinita():
    for Cr in (0.0, 0.5, 1.0):
        assert eps_contracorriente(200.0, Cr) == pytest.approx(1.0, abs=1e-2)


def test_Cr_1_da_NTU_igual_eps_sobre_1_menos_eps():
    assert ntu_contracorriente(0.6, 1.0) == pytest.approx(0.6 / 0.4)


@pytest.mark.parametrize("Cr", [0.0, 0.2, 0.7, 0.999, 1.0])
def test_ntu_es_la_inversa_de_eps(Cr):
    for NTU in (0.3, 1.0, 3.0):
        assert ntu_contracorriente(eps_contracorriente(NTU, Cr), Cr) == pytest.approx(NTU, rel=1e-6)


def test_verificar_lo_disenado_devuelve_las_mismas_temperaturas():
    prod = Corriente("producto", 0.25, lambda T: pr.choi_okos({"agua": 0.77, "proteina": 0.23}, T))
    agua = Corriente("agua", 1.0, pr.agua)
    T_out = 60.5
    Q = prod.C(40.0, T_out) * (T_out - 40.0)
    Tw_out = temperatura_por_balance(agua, 63.0, -Q)
    d = disenar_seccion("cal", agua, prod, 63.0, Tw_out, 40.0, T_out, U=850.0)
    v = verificar_seccion("cal", agua, prod, 63.0, 40.0, d.A, U=850.0)
    assert v.Tf_out == pytest.approx(T_out, abs=1e-4)
    assert v.Tc_out == pytest.approx(Tw_out, abs=1e-4)


def test_balance_de_energia_cierra():
    prod = Corriente("producto", 0.3, lambda T: pr.choi_okos({"agua": 0.75, "grasa": 0.25}, T))
    agua = Corriente("agua", 0.8, pr.agua)
    s = verificar_seccion("cal", agua, prod, 65.0, 20.0, A=2.0, U=800.0)
    Q_caliente = agua.C(s.Tc_in, s.Tc_out) * (s.Tc_in - s.Tc_out)
    Q_fria = prod.C(s.Tf_in, s.Tf_out) * (s.Tf_out - s.Tf_in)
    assert Q_caliente == pytest.approx(Q_fria, rel=1e-4)


def test_cruce_de_temperaturas_se_rechaza():
    c = Corriente("c", 1.0, pr.agua)
    with pytest.raises(ValueError):
        disenar_seccion("x", c, c, 50.0, 30.0, 20.0, 55.0, U=800.0)
