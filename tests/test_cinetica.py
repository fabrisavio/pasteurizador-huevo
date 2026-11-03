"""Tests del modelo D–z: las definiciones de D y z se tienen que cumplir."""

import pytest

from pasteurizador.cinetica import Cinetica

c = Cinetica(D_ref_min=1.25, T_ref_C=60.0, z_C=4.84)


def test_D_en_la_referencia():
    assert c.D(60.0) == pytest.approx(1.25)


def test_subir_z_grados_divide_D_por_10():
    assert c.D(60.0 + 4.84) == pytest.approx(0.125)


def test_logs_proporcionales_al_tiempo():
    assert c.logs(2.5, 60.0) == pytest.approx(2.0)
    assert c.logs(5.0, 60.0) == pytest.approx(2 * c.logs(2.5, 60.0))
