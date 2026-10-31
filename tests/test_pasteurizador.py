from pathlib import Path

import pytest

from pasteurizador import config
from pasteurizador.pasteurizador import correr
from pasteurizador.validacion import validar

RAIZ = Path(__file__).parents[1]


@pytest.fixture(scope="module")
def resultado():
    return correr(config.cargar(RAIZ / "config" / "base.yaml"))


def test_validacion_gut2003_dentro_de_tolerancia():
    res = validar(RAIZ / "config" / "validacion_gut2003.yaml")
    assert res.ok, res.tabla()


def test_el_equipo_cubre_a_cada_producto(resultado):
    eq = resultado.equipo
    for d in eq.disenos.values():
        if d.regen:
            assert eq.A_regen >= d.regen.A
        assert eq.A_calent >= d.calent.A


def test_la_regeneracion_la_gobierna_la_clara(resultado):
    # entra a 35 °C y necesita ε ≈ 0,6 contra ~0,3 de huevo y yema
    assert resultado.equipo.gobierna_regen == "CP"


def test_todos_llegan_a_T_holding_sin_pasarse_del_tope_de_agua(resultado):
    cfg = resultado.cfg
    for k, vs in resultado.verificaciones.items():
        p = cfg.productos[k]
        for v in vs:
            assert v.calent is not None
            assert v.calent.Tf_out == pytest.approx(p.T_holding + cfg.perdida_holding, abs=1e-3)
            assert v.T_agua_necesaria <= p.T_max_agua + 1e-6


def test_nadie_sale_mas_caliente_que_lo_pedido(resultado):
    for k, vs in resultado.verificaciones.items():
        for v in vs:
            assert v.T_salida <= resultado.cfg.productos[k].T_salida + 0.05


def test_ninguna_particula_pasa_menos_que_t_holding(resultado):
    for k, r in resultado.equipo.tubo.elegido.productos.items():
        t_min = resultado.cfg.productos[k].t_holding / 60
        assert r.t_rapida_min >= t_min - 1e-9


def test_corrida_base_cumple(resultado):
    assert resultado.cumple, resultado.alertas


def test_base_enrollada_descarta_DN32(resultado):
    # en recto DN32 entraba con 2 % de margen; enrollado, la fricción extra lo saca
    tubo = resultado.equipo.tubo
    assert tubo.geometria.es_espiral
    assert not next(op for op in tubo.opciones if op.DN == 32).factible
    assert tubo.elegido.DN == 40
