# -*- coding: utf-8 -*-
"""Test sintético de las métricas de copa (lli_core.metricas_copa y perfil_anchura) calculadas A MANO.

Franjas de 0,1 m desde 0,3 m (franja j = [0,3 + 0,1 j, 0,4 + 0,1 j)); puntos a la altura central de cada franja y repartidos
en t de -w/2 a +w/2 (anchura estricta = w).
  Rebanada 0 (s = 0,05): franjas 2..7 con anchuras 0,6 0,7 0,8 0,9 1,0 0,5 (10 puntos por franja, 60 en total).
      área = 0,1 x 4,5 = 0,45 m² ; volumen = 0,045 m³ ; base de copa = 0,5 m ; MH = 1,05 ; profundidad = 0,55
      anchura máxima = 1,0 en la franja 6 (centro 0,95 m) ; MW = 4,5 / 6 = 0,75 ; TRV (4 m) = 1,05 x 0,75 x 10 000 / 4 = 1968,75
  Rebanada 1 (s = 0,15): franjas 0, 1, 3 y 4 con anchura 1,0 (15 puntos por franja); la franja 2 está vacía (hueco).
      área = 0,4 ; base = 0,3 ; MH = 0,75 ; profundidad = 0,45 ; anchura máxima 1,0 en la franja 0 (centro 0,35 m)
  Perfil de la hilera (2 rebanadas): franja 0 -> ocupada en 1 de 2, anchura 1,0, contando huecos 0,5;
      franja 2 -> solo la rebanada 0 (0,6), contando huecos 0,3; franja 3 -> (0,7 + 1,0) / 2 = 0,85.
Uso:  python tests\\test_metricas_copa.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils  # noqa: F401,E402
import numpy as np  # noqa: E402
import lli_core as lc  # noqa: E402


def construir():
    s, t, h = [], [], []
    for j, w in zip(range(2, 8), (0.6, 0.7, 0.8, 0.9, 1.0, 0.5)):
        for tt in np.linspace(-w / 2, w / 2, 10):
            s.append(0.05); t.append(tt); h.append(0.3 + 0.1 * j + 0.05)
    for j in (0, 1, 3, 4):
        for tt in np.linspace(-0.5, 0.5, 15):
            s.append(0.15); t.append(tt); h.append(0.3 + 0.1 * j + 0.05)
    return np.array(s), np.array(t), np.array(h)


def test_metricas_copa():
    s, t, h = construir()
    P = lc.perfiles_rebanada(s, t, h, s_ini=0.0, n_reb=2, min_pts_franja=5, min_pts_reb=50, n_franjas=10)
    m = lc.metricas_copa(P, 4.0)
    assert np.allclose(m["area_seccion_m2"], [0.45, 0.4]), m["area_seccion_m2"]
    assert np.allclose(m["vol_copa_m3"], [0.045, 0.04]), m["vol_copa_m3"]
    assert np.allclose(m["h_base_copa_m"], [0.5, 0.3]), m["h_base_copa_m"]
    assert np.allclose(P.MH_max, [1.05, 0.75]), P.MH_max
    assert np.allclose(m["prof_copa_m"], [0.55, 0.45]), m["prof_copa_m"]
    assert np.allclose(m["anchura_max_m"], [1.0, 1.0]), m["anchura_max_m"]
    assert np.allclose(m["h_anchura_max_m"], [0.95, 0.35]), m["h_anchura_max_m"]
    assert np.isclose(m["TRV_m3_ha"][0], 1968.75), m["TRV_m3_ha"]
    assert np.isnan(lc.metricas_copa(P)["TRV_m3_ha"]).all()          # sin distancia entre hileras, sin TRV
    perf = {f["franja"]: f for f in lc.perfil_anchura(P)}
    assert sorted(perf) == list(range(8)), sorted(perf)
    assert perf[0]["n_ocupadas"] == 1 and np.isclose(perf[0]["anchura_media_m"], 1.0)
    assert np.isclose(perf[0]["anchura_media_con_huecos_m"], 0.5)
    assert np.isclose(perf[2]["anchura_media_m"], 0.6) and np.isclose(perf[2]["anchura_media_con_huecos_m"], 0.3)
    assert np.isclose(perf[3]["anchura_media_m"], 0.85) and perf[3]["frac_ocupada"] == 1.0
    assert np.isclose(perf[0]["h_ini_m"], 0.3) and np.isclose(perf[7]["h_fin_m"], 1.1)
    por_tramo = lc.perfil_anchura(P, np.array([0, 1]))
    assert {f["grupo"] for f in por_tramo} == {0, 1}
    assert max(f["franja"] for f in por_tramo if f["grupo"] == 1) == 4      # el tramo 1 llega hasta la franja 4
    print("OK: área %s, base %s, TRV %.2f, perfil franja 3 = %.3f m" % (m["area_seccion_m2"], m["h_base_copa_m"],
                                                                      m["TRV_m3_ha"][0], perf[3]["anchura_media_m"]))


if __name__ == "__main__":
    test_metricas_copa()
