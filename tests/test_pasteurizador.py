"""Tests del modelo completo con la config base y la validación contra Gut.

El fixture `resultado` corre el diseño una sola vez y se lo pasa a cada test
que lo pide como argumento (scope="module": una vez por archivo).
"""

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
            T_objetivo = cfg.T_operacion(p) + resultado.equipo.perdida_tubo[k]
            assert v.calent.Tf_out == pytest.approx(T_objetivo, abs=1e-3)
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
    # DN32 entra raspando en Δp; con el factor de seguridad sobre K queda afuera
    tubo = resultado.equipo.tubo
    assert tubo.geometria.es_espiral
    assert not next(op for op in tubo.opciones if op.DN == 32).factible
    assert tubo.elegido.DN == 40


def test_temperaturas_del_recorrido_en_orden(resultado):
    # entra < sale de regeneración < entra al tubo, y al secador más frío que el tubo
    for k in resultado.cfg.productos:
        v = resultado.en_diseno(k)
        assert v.T_entrada < v.T_regen_salida < v.T_calent_salida
        assert v.T_salida < resultado.cfg.productos[k].T_holding


def test_config_de_validacion_no_se_confunde_con_una_de_diseno():
    with pytest.raises(config.ConfigInvalida):
        config.cargar(RAIZ / "config" / "validacion_gut2003.yaml")


def test_el_bloque_ahorra_placas_frente_a_la_coquilla(resultado):
    eq = resultado.equipo
    assert eq.aislacion["bloque"].E_MWh_anio < eq.aislacion["tubo"].E_MWh_anio
    assert eq.placas_calent_por_aislacion["bloque"] <= eq.placas_calent_por_aislacion["tubo"]
    assert eq.placas_calent == eq.placas_calent_por_aislacion["bloque"]  # base: bloque


def test_el_tubo_sale_a_la_T_legal_mas_el_margen(resultado):
    cfg = resultado.cfg
    for p in cfg.productos.values():
        assert cfg.T_operacion(p) == pytest.approx(p.T_holding + cfg.margen_control)
        assert cfg.T_operacion(p) >= p.T_holding


def test_sin_seccion_de_aislacion_usa_la_perdida_fija(tmp_path):
    # tmp_path: carpeta temporal que pytest crea y borra solo
    import yaml
    with open(RAIZ / "config" / "base.yaml", encoding="utf-8") as f:
        c = yaml.safe_load(f)
    del c["equipo"]["holding"]["aislacion"]
    c["equipo"]["holding"]["perdida_T_C"] = 0.5
    ruta = tmp_path / "fija.yaml"
    ruta.write_text(yaml.safe_dump(c, allow_unicode=True), encoding="utf-8")
    r = correr(config.cargar(ruta))
    assert all(d == 0.5 for d in r.equipo.perdida_tubo.values())
    assert r.equipo.aislacion == {}


def test_clara_por_la_ruta_de_la_soda(resultado):
    # pH 9,3 y el par de FSIS Tabla III.A: 56,7 °C × 4,96 min
    cp = resultado.cfg.productos["CP"]
    assert cp.pH == pytest.approx(9.3)
    assert cp.T_holding == pytest.approx(56.7)
    assert cp.t_holding == pytest.approx(4.96 * 60)
    # la partícula más rápida logra al menos los 5,7 log de la tabla
    assert resultado.equipo.tubo.elegido.productos["CP"].logs_min >= 5.7 - 0.1


def test_agua_de_la_clara_bajo_59(resultado):
    for v in resultado.verificaciones["CP"]:
        assert v.T_agua_necesaria <= 59.0


def test_aviso_de_espesamiento_solo_donde_corresponde(resultado):
    espesa = [a for a in resultado.avisos if "empieza a espesar" in a]
    claves = {a.split(":")[0] for a in espesa}
    assert claves == {"HP", "CP"}   # la yema recién espesa a 65 °C
