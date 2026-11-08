"""Tests de las placas (v1): geometría, correlaciones contra valores de
referencia y sentido físico del U y del Δp.

Los valores de referencia de Kumar (fricción) y Martin salen de la biblioteca
fluids 1.3.1 de Caleb Bell, que implementa las mismas tablas; sirven para
detectar un error de tipeo en las constantes.
"""

import math

import pytest

from pasteurizador import placas as pl
from pasteurizador import propiedades as pr
from pasteurizador.seccion_phe import Corriente, verificar_seccion

PLACA = pl.Placa(area=0.15, ancho=0.20, b=0.003, angulo=45.0, phi=1.17,
                 espesor=0.0005, k=16.3, D_puerto=0.05)
CLARA = pl.Fluido("CP", lambda T: pr.choi_okos({"agua": 0.77, "proteina": 0.23}, T),
                  lambda T: pr.reologia("clara_uf_thapon_scalzo", T))


# --- geometría ---------------------------------------------------------------
def test_geometria_del_canal():
    assert PLACA.Dh == pytest.approx(2 * 0.003 / 1.17)
    assert PLACA.De == pytest.approx(0.006)
    assert PLACA.L == pytest.approx(0.15 / (1.17 * 0.20))


def test_placas_de_un_arreglo():
    # 3 pasos de 8 canales por lado: 48 canales en total, 47 placas que los separan
    assert pl.Arreglo(3, 8).placas == 47
    assert pl.Arreglo(1, 1).placas == 1


# --- viscosidad generalizada ------------------------------------------------
def test_viscosidad_generalizada_newtoniana_es_K():
    assert pl.viscosidad_generalizada(0.01, 1.0, 0.3, 0.005) == pytest.approx(0.01)
    assert pl.tasa_corte_pared(1.0, 0.3, 0.005) == pytest.approx(12 * 0.3 / 0.005)


def test_pseudoplastico_mas_rapido_es_menos_viscoso():
    lento = pl.viscosidad_generalizada(0.2, 0.86, 0.05, 0.005)
    rapido = pl.viscosidad_generalizada(0.2, 0.86, 0.5, 0.005)
    assert rapido < lento


# --- correlaciones -------------------------------------------------------------
@pytest.mark.parametrize("Re, f_darcy", [(5, 37.6), (100, 3.633069423730658),
                                         (3000, 1.107737845948731)])
def test_friccion_de_kumar_contra_fluids(Re, f_darcy):
    assert pl.f_kumar(Re, 45) == pytest.approx(f_darcy / 4, rel=1e-9)


def test_nusselt_de_kumar_en_sus_tramos():
    # Re < 10: C = 0,718, y = 0,349 ; 10–100: 0,400 / 0,598 ; > 100: 0,300 / 0,663 (45°)
    assert pl.nu_kumar(5, 500, 45) == pytest.approx(0.718 * 5**0.349 * 500**(1 / 3))
    assert pl.nu_kumar(50, 100, 45) == pytest.approx(0.400 * 50**0.598 * 100**(1 / 3))
    assert pl.nu_kumar(2000, 3, 45) == pytest.approx(0.300 * 2000**0.663 * 3**(1 / 3))
    # ht usa Pr^0,33 en vez de 1/3: coinciden dentro del 3 %
    assert pl.nu_kumar(5, 500, 45) == pytest.approx(9.78870131970703, rel=0.03)


@pytest.mark.parametrize("Re, xi, Nu", [(500, 1.0668845111588283, 18.82546548662704),
                                        (5000, 0.8341738282706309, 96.11247602778299)])
def test_martin_contra_fluids_y_ht(Re, xi, Nu):
    assert pl.xi_martin(Re, 45) == pytest.approx(xi, rel=1e-9)
    assert pl.nu_martin(Re, 3.0, 45) == pytest.approx(Nu, rel=1e-9)


def test_correccion_de_pared_sube_h_al_calentar():
    # calentando, la pared está más caliente y el fluido allí es menos viscoso
    assert pl.nu_kumar(50, 100, 45, mu_sobre_mu_pared=1.3) > pl.nu_kumar(50, 100, 45)


# --- un lado y el U ----------------------------------------------------------
def test_mas_canales_por_paso_baja_h_y_dp():
    a = pl.canal(CLARA, 0.23, pl.Arreglo(1, 4), PLACA, 45.0, 45.0)
    b = pl.canal(CLARA, 0.23, pl.Arreglo(1, 8), PLACA, 45.0, 45.0)
    assert b.v == pytest.approx(a.v / 2)
    assert b.h < a.h and b.dp < a.dp


def test_dp_de_canales_proporcional_a_los_pasos():
    uno = pl.canal(CLARA, 0.23, pl.Arreglo(1, 4), PLACA, 45.0, 45.0)
    tres = pl.canal(CLARA, 0.23, pl.Arreglo(3, 4), PLACA, 45.0, 45.0)
    assert tres.dp_canales == pytest.approx(3 * uno.dp_canales)


def test_U_es_resistencias_en_serie():
    lado = pl.Lado(CLARA, 0.23, pl.Arreglo(1, 6))
    det = pl.coeficiente_U(PLACA, lado, lado, 50.0, 40.0, 6e-4)
    esperado = 1 / (1 / det.caliente.h + 1 / det.fria.h + PLACA.R_pared + 6e-4)
    assert det.U == pytest.approx(esperado)
    assert det.U_limpio > det.U
    assert sum(det.fraccion_resistencias.values()) == pytest.approx(1.0)


def test_caudal_CIP_da_el_Re_pedido():
    arr = pl.Arreglo(3, 8)
    Q = pl.caudal_CIP(PLACA, arr, 2150, 70.0) / 3600   # m³/s
    t, mu = pr.agua(70.0), pr.agua_reologia(70.0).K
    Re = Q * PLACA.De * t.rho / (arr.canales * PLACA.b * PLACA.ancho * mu)   # Galiński, ec. 6
    assert Re == pytest.approx(2150)


def test_seccion_con_U_funcion_constante_igual_que_numero():
    prod = Corriente("p", 0.25, lambda T: pr.choi_okos({"agua": 0.77, "proteina": 0.23}, T))
    a = verificar_seccion("r", prod, prod, 60.0, 10.0, 2.0, 600.0)
    b = verificar_seccion("r", prod, prod, 60.0, 10.0, 2.0, lambda *T: 600.0)
    assert b.Tf_out == pytest.approx(a.Tf_out, abs=1e-9)


def test_agua_viscosidad_iapws():
    assert pr.agua_reologia(20.0).K == pytest.approx(1.0016e-3, rel=0.01)
    assert pr.agua_reologia(20.0).n == 1.0
    assert math.isclose(pr.agua_reologia(70.0).K, 4.04e-4, rel_tol=0.02)
