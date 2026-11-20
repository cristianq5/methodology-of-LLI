# -*- coding: utf-8 -*-
"""Test sintético de la envolvente rellena (2026-11-16).

Caso: una rebanada con 10 franjas (h = 0,3 ... 1,3 m). Las franjas 4 y 5 (h = 0,7-0,9 m) NO tienen puntos: hueco
interior de la copa. El resto, 6 puntos por franja entre t = -1 y t = 1.
Rayos desde la izquierda (t = +2 m hacia -t), 10 por franja, en las 10 franjas:
  - franjas con copa: 4 con eco en t = -0,5 (interceptados) y 6 sin eco (atraviesan)
  - franjas del hueco: los 10 sin eco (atraviesan)
A mano:
  envolvente actual  (relleno=False): el hueco no es copa -> Bt = 8 x 10 = 80, Bi = 8 x 6 = 48 -> P = 0,60
  envolvente rellena (relleno=True) : Bt = 10 x 10 = 100, Bi = 48 + 20 = 68 -> P = 0,68
  En el hueco, tmin/tmax interpolados entre las franjas 3 y 6 = -1 / 1.
Uso: python tests\\test_envolvente_rellena.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils  # noqa: F401,E402
import numpy as np  # noqa: E402
import lli_core as lc  # noqa: E402

HUECO = (4, 5)


def construir():
    s, t, h = [], [], []
    for j in range(10):
        if j in HUECO:
            continue
        for tt in (-1.0, -0.5, 0.0, 0.5, 0.75, 1.0):
            s.append(0.05); t.append(tt); h.append(0.3 + 0.1 * j + 0.05)
    return np.array(s), np.array(t), np.array(h)


def rayos():
    O, D, Rr = [], [], []
    for j in range(10):
        z = 0.3 + 0.1 * j + 0.05
        for k in range(10):
            O.append((0.05, 2.0, z)); D.append((0.0, -1.0, 0.0))
            Rr.append(np.inf if (j in HUECO or k >= 4) else 2.5)
    return np.array(O), np.array(D), np.array(Rr)


def contar(P, relleno):
    env = lc.envolvente(P, "max", relleno=relleno)
    O, D, Rr = rayos()
    cl = lc.clasificar_rayos(O, D, Rr, env, P, dtm_fn=lambda s_, t_: np.zeros_like(s_), paso=0.025,
                             az_max_deg=25, el_rango_deg=(-60, 15), origen_t_abs=(0.8, 3.0), alcance_max=30.0)
    c = lc.acumular_conteos(cl, 1, 10)
    return env, int(c["Bt"][0].sum()), int(c["Bi"][0].sum())


def test_todo():
    s, t, h = construir()
    P = lc.perfiles_rebanada(s, t, h, s_ini=0.0, n_reb=1, min_pts_franja=5, min_pts_reb=40, n_franjas=10)
    assert not P.ocupada[0, 4] and not P.ocupada[0, 5]
    env0, bt0, bi0 = contar(P, False)
    env1, bt1, bi1 = contar(P, True)
    assert (bt0, bi0) == (80, 48), (bt0, bi0)
    assert (bt1, bi1) == (100, 68), (bt1, bi1)
    assert env1[2][0].all() and np.allclose(env1[0][0, 4:6], -1.0) and np.allclose(env1[1][0, 4:6], 1.0)
    # el relleno no toca la envolvente original ni lo que queda fuera (por encima/debajo)
    assert not P.ocupada[0, 4] and np.isnan(P.tmin[0, 4])
    print("OK: actual Bt=%d Bi=%d P=%.2f | rellena Bt=%d Bi=%d P=%.2f" % (bt0, bi0, bi0 / bt0, bt1, bi1, bi1 / bt1))


if __name__ == "__main__":
    test_todo()
