# -*- coding: utf-8 -*-
"""Tests de casos límite (sin datos reales):

1) Hold-out: los grupos temporales de la validación nunca coinciden con los del ajuste (k_ajuste).
2) Pasos del sensor: dos hileras vecinas digitalizadas en sentidos OPUESTOS comparten la misma ventana de la calle
   (la pasada no se parte por la mitad).
3) pcap: la resincronización encuentra el inicio de registro en un pcap sintético.
Uso:  python tests/test_casos_limite.py
"""
import os
import struct
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils  # noqa: F401,E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from shapely.geometry import LineString  # noqa: E402

import extrinseca as ex  # noqa: E402
import geometria as geo  # noqa: E402
import pcap_hesai as ph  # noqa: E402


def test_holdout_disjunto():
    rng = np.random.default_rng(0)
    n_k = 400                                         # 400 sellos de 1/1024 s con 80 puntos LAS y 2 bloques de pcap cada uno
    k = np.arange(10_000, 10_000 + n_k)
    las_t = np.repeat(k * ex.CUANTO_T, 80) + rng.uniform(-1e-4, 1e-4, n_k * 80)
    las_xyz = rng.normal(size=(len(las_t), 3))
    tb = np.repeat(k * ex.CUANTO_T, 2)
    az = rng.uniform(0, 360, len(tb))
    dist = rng.uniform(1, 10, (len(tb), 32))
    G1 = ex.grupos_temporales(las_t, las_xyz, tb, az, dist, 150, np.random.default_rng(42))
    G2 = ex.grupos_temporales(las_t, las_xyz, tb, az, dist, 300, np.random.default_rng(777), excluir=G1["k"])
    assert len(G1["k"]) == 150 and len(G2["k"]) == 250, (len(G1["k"]), len(G2["k"]))
    assert not set(G1["k"].tolist()) & set(G2["k"].tolist())
    print("OK hold-out: %d grupos de ajuste, %d de validación, 0 en común" % (len(G1["k"]), len(G2["k"])))


class _Tray:
    pass


def test_sentido_opuesto():
    dt = 0.005
    y = np.arange(-2, 52, dt)                         # el sensor camina por la calle x = 2 hacia +y, a 1 m/s
    tr = _Tray()
    tr.t = 1000 + np.arange(len(y)) * dt
    tr.p = np.c_[np.full(len(y), 2.0), y, np.zeros(len(y))]
    cfg = utils.cargar_parametros()
    cfg["recorte"].update(modo="completo", filas_prueba_id=[], tramo_prueba_s_m=[0, 1e9])
    cfg["rutas"] = dict(filas_gpkg="sintetico")
    res = {}
    for nom, l2 in (("mismo", LineString([(4, 0), (4, 50)])), ("opuesto", LineString([(4, 50), (4, 0)]))):
        filas = pd.DataFrame(dict(id_fila=["1", "2"], geometry=[LineString([(0, 0), (0, 50)]), l2]))
        geo._CACHE_FILAS.clear()
        geo._CACHE_FILAS["sintetico|completo"] = filas
        vent, pasos = geo.ventanas_tiempo(cfg, tr, filas)
        res[nom] = vent
        assert len(vent) == 1, (nom, vent)
        assert sorted(p["fila"] for p in pasos) == ["1", "2"], pasos
        # el sentido de cada paso es el de la propia hilera
        sent = {p["fila"]: p["sentido"] for p in pasos}
        assert sent["1"] == 1 and sent["2"] == (1 if nom == "mismo" else -1), sent
    geo._CACHE_FILAS.clear()
    assert np.allclose(res["mismo"], res["opuesto"])
    print("OK sentido: una sola ventana %s con la segunda hilera en cualquier sentido" % np.round(res["opuesto"], 1).tolist())


def _registro(t_host, payload_extra=b""):
    carga = b"\xee\xff" + bytes(818)                    # carga UDP de 820 B que empieza por la firma 0xEE 0xFF
    trama = bytes(ph.OFF_ETH) + carga                    # 42 B de Ethernet/IP/UDP + 820 B = 862 B
    ts = int(t_host)
    return struct.pack("<IIII", ts, int(round((t_host - ts) * 1e6)), 862, 862) + trama


def test_resincronizar_pcap():
    cab = bytes(24)
    regs = b"".join(_registro(100.0 + 0.001 * i) for i in range(50))
    with tempfile.TemporaryDirectory() as d:
        ruta = os.path.join(d, "s.pcap")
        with open(ruta, "wb") as f:
            f.write(cab + regs)
        with open(ruta, "rb") as f:
            o = ph._resincronizar(f, 24 + 5 * ph.TAM_REG + 100)       # desde la mitad del registro 5
            assert o == 24 + 6 * ph.TAM_REG, o
            o2 = ph.buscar_offset(f, 100.030, os.path.getsize(ruta))
            assert ph._tiempo_host(f, o2) <= 100.030 + 1e-9 and o2 > 24, o2
    print("OK pcap: resincroniza en el registro 6 (offset %d) y la bisección avanza (offset %d)" % (o, o2))


if __name__ == "__main__":
    test_holdout_disjunto()
    test_sentido_opuesto()
    test_resincronizar_pcap()
