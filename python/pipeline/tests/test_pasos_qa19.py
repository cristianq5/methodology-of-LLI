# -*- coding: utf-8 -*-
"""QA-19: los pasos del sensor se cortan por cambio de lado y de sentido y nunca se fusionan; las ventanas no se solapan.

1) Trayectoria sintética: calle t=+2 (sentido +s), media vuelta por la cabecera, calle t=-2 (sentido -s), con la media
   vuelta DENTRO de s1+pad y |t|<4,5 y sin hueco temporal: deben salir exactamente 2 pasos de un solo lado y sentido.
2) Trayectoria real (b1) en modo completo: cada hilera tiene pasos de los dos lados, cada paso un solo lado y sentido, las ventanas
   no se solapan y cada ventana es un único paso.
Uso:  python tests\\test_pasos_qa19.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils  # noqa: F401,E402
import numpy as np  # noqa: E402
import geometria as geo  # noqa: E402


def test_sintetico():
    dt = 0.005
    tr = []
    for s in np.arange(0, 40, 1.0 * dt):          # 1 m/s por la calle t=+2
        tr.append((s, 2.0))
    for k in np.arange(0, 1.0, dt):               # media vuelta de 4 m en 1 s por la cabecera (s=40..42)
        tr.append((40 + 2 * np.sin(np.pi * k), 2.0 - 4.0 * k))
    for s in np.arange(40, 0, -1.0 * dt):         # vuelta por t=-2
        tr.append((s, -2.0))
    tr = np.array(tr)
    t = 1.0e9 + np.arange(len(tr)) * dt
    pasos = geo.segmentar_pasos(t, tr[:, 0], tr[:, 1], 1.0, 39.0, 2.0, (0.8, 4.5), 1.0)
    assert len(pasos) == 2, pasos
    assert [p["lado"] for p in pasos] == [1, -1] and [p["sentido"] for p in pasos] == [1, -1], pasos
    print("OK sintetico: %d pasos, lados %s, sentidos %s" % (len(pasos), [p["lado"] for p in pasos], [p["sentido"] for p in pasos]))


def test_real():
    if not os.environ.get(utils.VAR_ENTORNO_ENSAYO):
        print("(sin ficha de ensayo en %s: se omite el test real)" % utils.VAR_ENTORNO_ENSAYO)
        return
    cfg = utils.cargar_config()
    cache = os.path.join(cfg["rutas"]["salida"], "01_qc", "trayectoria_utm.npz")
    if not os.path.exists(cache):
        print("(sin trayectoria en caché: se omite el test real)")
        return
    from trayectoria import cargar_trayectoria
    cfg["recorte"]["modo"] = "completo"
    tray, _ = cargar_trayectoria(cfg, cache_dir=os.path.dirname(cache))
    filas = utils.cargar_filas(cfg, todas=True)
    for _, f in filas.iterrows():
        m = utils.MarcoHilera.de_linea(f.geometry)
        s0, s1 = utils.tramo_s(cfg, m)
        ps = geo.pasos_sensor(cfg, tray, m, s0, s1)
        assert {p["lado"] for p in ps} == {-1, 1}, "hilera %s: pasos %s" % (f["id_fila"], [(p["lado"], p["sentido"], round(p["t1"] - p["t0"])) for p in ps])
        assert all(p["lado"] in (-1, 1) and p["sentido"] in (-1, 1) for p in ps)   # un lado y un sentido por paso
        for p in ps:
            assert p["t1"] - p["t0"] < 90, "hilera %s: paso demasiado largo (¿fusión de calles?)" % f["id_fila"]
    glob = sorted(geo.ventanas_globales(cfg, tray), key=lambda g: g["t0"])
    for g in glob:
        lados = {l for _, l in g["sirve"]}
        assert len(g["sirve"]) <= 2 and len({f for f, _ in g["sirve"]}) == len(g["sirve"]), g
    ven, pasos = geo.ventanas_tiempo(cfg, tray, filas)
    for a, b in zip(ven[:-1], ven[1:]):
        assert a[1] <= b[0] + 1e-9, "ventanas solapadas %s %s" % (a, b)
    print("OK real: %d hileras con pasos de ambos lados (un lado y sentido por paso); %d ventanas globales sin solape" % (len(filas), len(ven)))


if __name__ == "__main__":
    test_sintetico()
    test_real()
