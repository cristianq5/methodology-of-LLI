# -*- coding: utf-8 -*-
"""Test sintético de MH, MW, P y LLI calculados A MANO.

Caso: una rebanada (s en [0, 0,1)) con una copa de 10 franjas (h = 0,3 ... 1,3 m). En cada franja
6 puntos con t = {-1, -0,5, 0, 0,5, 0,75, 1} a la altura central de la franja (0,35 ... 1,25 m).

A mano:
  MH estricta = 1,25 m (punto más alto) ; MH percentil 99 = 1,25 (los 6 puntos más altos valen 1,25)
  MW estricta = 2,0 m (t_max - t_min en todas las franjas)
  MW percentil (2,5-97,5): p97,5 = 0,75 + 0,875*0,25 = 0,96875 ; p2,5 = -1 + 0,125*0,5 = -0,9375
                           anchura = 1,90625
  Rayos izquierda (origen t=+2 m, hacia -t): por franja 10 rayos, 6 sin eco (atraviesan) y 4 con eco en t=-0,5
      -> Bt=10, Bi=6, P_izq = 0,6
  Rayos derecha (origen t=-2 m, hacia +t): por franja 10 rayos, 8 sin eco y 2 con eco en t=+0,5 -> P_der = 0,8
  P = media de lados = 0,7 ; L = 0,3
  LLI estricto   = 1,25 * 2,0     * 0,3 = 0,75
  LLI percentil  = 1,25 * 1,90625 * 0,3 = 0,7148437...
Filtro +-2,5 DE: 20 valores 1,0 y uno 10,0 -> el 10,0 se descarta (4,47 DE) y la media queda 1,0.
Uso:  python tests\\test_lli_sintetico.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils  # noqa: F401,E402
import numpy as np  # noqa: E402
import lli_core as lc  # noqa: E402


def construir():
    s, t, h = [], [], []
    for j in range(10):
        for tt in (-1.0, -0.5, 0.0, 0.5, 0.75, 1.0):
            s.append(0.05)
            t.append(tt)
            h.append(0.3 + 0.1 * j + 0.05)
    return np.array(s), np.array(t), np.array(h)


def rayos():
    O, D, Rr = [], [], []
    for j in range(10):
        z = 0.3 + 0.1 * j + 0.05
        for k in range(10):                      # izquierda: 4 con eco, 6 sin eco
            O.append((0.05, 2.0, z)); D.append((0.0, -1.0, 0.0)); Rr.append(2.5 if k < 4 else np.inf)
        for k in range(10):                      # derecha: 2 con eco, 8 sin eco
            O.append((0.05, -2.0, z)); D.append((0.0, 1.0, 0.0)); Rr.append(2.5 if k < 2 else np.inf)
    return np.array(O), np.array(D), np.array(Rr)


def test_todo():
    s, t, h = construir()
    P = lc.perfiles_rebanada(s, t, h, s_ini=0.0, n_reb=1, min_pts_franja=5, min_pts_reb=50,
                             mh_percentil=99.0, mw_percentil_t=(2.5, 97.5), n_franjas=10)
    assert np.isclose(P.MH_max[0], 1.25), P.MH_max
    assert np.isclose(P.MH_p[0], 1.25), P.MH_p
    assert np.isclose(P.MW_max[0], 2.0), P.MW_max
    assert np.isclose(P.MW_p[0], 1.90625), P.MW_p
    env = lc.envolvente(P, "max")
    O, D, Rr = rayos()
    cl = lc.clasificar_rayos(O, D, Rr, env, P, dtm_fn=lambda s_, t_: np.zeros_like(s_), paso=0.025,
                             az_max_deg=25, el_rango_deg=(-60, 15), origen_t_abs=(0.8, 3.0), alcance_max=30.0)
    cont = lc.acumular_conteos(cl, 1, 10)
    assert cont["Bt"][0].sum() == 100 and cont["Bi"][0].sum() == 60      # izquierda
    assert cont["Bt"][1].sum() == 100 and cont["Bi"][1].sum() == 80      # derecha
    Pl, nfi, Pm = lc.p_por_rebanada(cont["Bt"], cont["Bi"], min_rayos_celda=10, min_franjas=3)
    assert np.allclose(Pl[:, 0], [0.6, 0.8]), Pl
    assert np.isclose(Pm[0], 0.7), Pm
    assert np.isclose(lc.lli(P.MH_max[0], P.MW_max[0], Pm[0]), 0.75)
    assert np.isclose(lc.lli(P.MH_p[0], P.MW_p[0], Pm[0]), 1.25 * 1.90625 * 0.3)
    x = np.r_[np.ones(20), 10.0]
    ok = lc.filtro_atipicos(x, 2.5)
    assert ok.sum() == 20 and not ok[-1] and np.isclose(x[ok].mean(), 1.0)
    print("OK: MH=%.3f MW=%.3f P=%.2f LLI=%.4f (esperado 0,75); LLI percentil=%.4f"
          % (P.MH_max[0], P.MW_max[0], Pm[0], lc.lli(P.MH_max[0], P.MW_max[0], Pm[0]),
             lc.lli(P.MH_p[0], P.MW_p[0], Pm[0])))


if __name__ == "__main__":
    test_todo()
