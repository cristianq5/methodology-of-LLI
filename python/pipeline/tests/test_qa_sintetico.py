# -*- coding: utf-8 -*-
"""Caso sintético reconstruido del informe QA (§2) + pruebas de QA-01, QA-05 y QA-04.

2 rebanadas (la segunda vacía), DTM plano a 0,2 m (h = z - 0,2). Franjas de 0,1 m desde 0,3 m.
  Franja 0 (h 0,35): 6 pts t = -0,05 (x3) y +0,05 (x3)      -> anchura 0,10 (percentil 0,10)
  Franjas 1-3 (h 0,45 0,55 0,65): 6 pts t = -0,4 ... 0,6     -> 1,00 (percentil 0,95)
  Franja 4 (h 0,75): 5 pts t = -0,1 ... 0,3                  -> 0,40 (percentil 0,38)
  Franja 5 (h 0,85): 3 pts (< 5) -> no cuenta en MW, sí en MH. 10 pts a h = 0,25 (< umbral): se ignoran.
  MH = 0,85 (p99 = 0,85) ; MW = (0,1 + 3*1,0 + 0,4)/5 = 0,70 ; MW percentil = (0,1 + 3*0,95 + 0,38)/5 = 0,666
Rayos (envolvente de las franjas 1-3: t en [-0,4; 0,6]):
  Izq. franja 2: 5 con eco dentro, 3 sin eco, 2 con eco más allá, 1 con eco ANTES de la copa -> Bt=10, Bi=5
  Izq. franja 1: 4 con eco dentro, 16 sin eco -> Bt=20, Bi=16
  Der. franja 3: 9 con eco dentro, 1 sin eco -> Bt=10, Bi=1
  3 rayos que se deben descartar (azimut 40°, elevación +30°, origen a |t| = 0,5).
  P izq = 21/30 = 0,70 ; P der = 0,10 ; P = 0,40 ; LLI = 0,85*0,70*0,60 = 0,357 ; LLI percentil = 0,85*0,666*0,6 = 0,33966
Sensibilidad 0,6 m (franjas j>=3): izq solo franja 3 (sin rayos) -> NaN con min_franjas=1 si no hay rayos en j>=3.
Cobertura (QA-05): izq. franjas fiables 1 y 2 -> techo 0,3+3*0,1 = 0,60 m ; der. franja 3 -> techo 0,70 m.
Uso:  python tests\\test_qa_sintetico.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils  # noqa: F401,E402
import numpy as np  # noqa: E402
import lli_core as lc  # noqa: E402

DTM = 0.2


def puntos():
    s, t, h = [], [], []

    def add(hag, ts):
        for tt in ts:
            s.append(0.05)
            t.append(tt)
            h.append(hag)
    add(0.35, [-0.05] * 3 + [0.05] * 3)
    for hag in (0.45, 0.55, 0.65):
        add(hag, [-0.4, -0.2, 0.0, 0.2, 0.4, 0.6])
    add(0.75, [-0.1, 0.0, 0.1, 0.2, 0.3])
    add(0.85, [0.0, 0.1, 0.2])
    add(0.25, [0.0] * 10)
    return np.array(s), np.array(t), np.array(h)


def rayos():
    O, D, R = [], [], []

    def add(lado, hag, Rr, n, az=0.0, el=0.0, to=2.0):
        for _ in range(n):
            O.append((0.05, lado * to, hag + DTM))
            ce = np.cos(np.radians(el))
            D.append((ce * np.sin(np.radians(az)), -lado * ce * np.cos(np.radians(az)), np.sin(np.radians(el))))
            R.append(Rr)
    add(+1, 0.55, 2.2, 5)            # eco dentro (t = -0,2)
    add(+1, 0.55, np.inf, 3)         # sin eco
    add(+1, 0.55, 3.0, 2)            # eco más allá (t = -1,0)
    add(+1, 0.55, 1.0, 1)            # eco ANTES de la copa (t = +1,0)
    add(+1, 0.45, 2.2, 4)
    add(+1, 0.45, np.inf, 16)
    add(-1, 0.65, 2.5, 9)            # der.: eco dentro (t = +0,5)
    add(-1, 0.65, np.inf, 1)
    add(+1, 0.55, 2.2, 1, az=40.0)   # descartado por azimut
    add(+1, 0.55, 2.2, 1, el=30.0)   # descartado por elevación
    add(+1, 0.55, 2.2, 1, to=0.5)    # descartado: origen a |t| = 0,5
    return np.array(O), np.array(D), np.array(R)


def test_nucleo():
    s, t, h = puntos()
    P = lc.perfiles_rebanada(s, t, h, 0.0, 2, min_pts_franja=5, min_pts_reb=10, n_franjas=8)
    assert np.isclose(P.MH_max[0], 0.85) and np.isclose(P.MH_p[0], 0.85)
    assert np.isclose(P.MW_max[0], 0.70), P.MW_max
    assert np.isclose(P.MW_p[0], 0.666), P.MW_p
    assert np.isnan(P.MH_max[1]) and np.isnan(P.MW_max[1])                 # rebanada vacía
    # envolvente: se fuerza la de las franjas 1-3 para los rayos
    env = lc.envolvente(P, "max")
    O, D, Rr = rayos()
    cl = lc.clasificar_rayos(O, D, Rr, env, P, dtm_fn=lambda a, b: np.full_like(a, DTM), paso=0.025,
                             az_max_deg=25, el_rango_deg=(-60, 15), origen_t_abs=(0.8, 3.0), alcance_max=30.0, tol_extremo=0.025)
    cont = lc.acumular_conteos(cl, 2, 8)
    Bt, Bi = cont["Bt"], cont["Bi"]
    assert Bt[0, 0, 2] == 10 and Bi[0, 0, 2] == 5, (Bt[0, 0], Bi[0, 0])
    assert Bt[0, 0, 1] == 20 and Bi[0, 0, 1] == 16
    assert Bt[1, 0, 3] == 10 and Bi[1, 0, 3] == 1
    assert Bt.sum() == 40, Bt.sum()                                         # los 3 descartados no cuentan
    Pl, nfi, Pm = lc.p_por_rebanada(Bt, Bi, min_rayos_celda=10, min_franjas=1)
    assert np.isclose(Pl[0, 0], 0.7) and np.isclose(Pl[1, 0], 0.1) and np.isclose(Pm[0], 0.4), (Pl, Pm)
    assert np.isclose(lc.lli(P.MH_max[0], P.MW_max[0], Pm[0]), 0.357)
    assert np.isclose(lc.lli(P.MH_p[0], P.MW_p[0], Pm[0]), 0.85 * 0.666 * 0.6)
    assert np.isnan(Pm[1])                                                  # sin rayos -> NaN
    # solo un lado -> NaN
    Bt1 = Bt.copy(); Bt1[1] = 0
    assert np.isnan(lc.p_por_rebanada(Bt1, Bi, 10, 1)[2][0])
    # QA-05: cobertura vertical
    hm, fr, frm = lc.cobertura_P(Bt, P.MH_max, 0.3, 0.1, min_rayos_celda=10)
    assert np.isclose(hm[0, 0], 0.60) and np.isclose(hm[1, 0], 0.70), hm
    assert np.isclose(fr[0, 0], 2 * 0.1 / (0.85 - 0.3)), fr
    # QA-04: variante desde 0,6 m (j >= 3): izquierda no tiene franja fiable -> NaN
    Pl6, _, Pm6 = lc.p_por_rebanada(Bt, Bi, 10, 1, j_min=3)
    assert np.isnan(Pl6[0, 0]) and np.isclose(Pl6[1, 0], 0.1) and np.isnan(Pm6[0])
    # filtro de atípicos
    x = np.r_[np.ones(30), 1.1, 1.05, 0.95, 9.0]
    ok = lc.filtro_atipicos(x, 2.5)
    assert ok.sum() == 33 and not ok[-1]
    print("OK nucleo: MH=%.2f MW=%.3f P=%.2f LLI=%.4f" % (P.MH_max[0], P.MW_max[0], Pm[0], lc.lli(P.MH_max[0], P.MW_max[0], Pm[0])))


def test_qa01_tramo_completo():
    """Modo prueba: tramo_prueba_s_m. Modo completo con hilera.extension = "gpkg": la línea entera menos las cabeceras.
    Modo completo con hilera.extension = "nube" (por defecto): zona de búsqueda = la línea ampliada margen_busqueda_m por
    cada extremo (la longitud real la fija 08 con la copa)."""
    cfg = utils.cargar_parametros()
    cfg["recorte"]["tramo_prueba_s_m"] = [18.0, 30.0]
    m = utils.MarcoHilera([0.0, 0.0], [0.0, 1.0], largo=46.0)
    c1 = dict(cfg, recorte=dict(cfg["recorte"], modo="prueba"))
    c_gpkg = dict(cfg, recorte=dict(cfg["recorte"], modo="completo"), hilera=dict(cfg["hilera"], extension="gpkg"))
    c_nube = dict(cfg, recorte=dict(cfg["recorte"], modo="completo"), hilera=dict(cfg["hilera"], extension="nube"))
    ext = cfg["hilera"]["recorte_extremos_m"]
    mg = cfg["hilera"]["margen_busqueda_m"]
    assert utils.tramo_s(c1, m) == tuple(cfg["recorte"]["tramo_prueba_s_m"])
    assert utils.tramo_s(c_gpkg, m) == (ext, 46.0 - ext)
    assert utils.tramo_s(c_nube, m) == (-mg, 46.0 + mg)
    print("OK QA-01: prueba %s, completo gpkg %s, completo nube %s" % (utils.tramo_s(c1, m), utils.tramo_s(c_gpkg, m),
                                                                      utils.tramo_s(c_nube, m)))


if __name__ == "__main__":
    test_nucleo()
    test_qa01_tramo_completo()
