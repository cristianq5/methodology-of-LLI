# -*- coding: utf-8 -*-
"""02_recorte: extrae del LAS original (solo lectura, por trozos) los bloques de trabajo.

Modo "prueba": un bloque con las hileras filas_prueba_id (tramo tramo_prueba_s_m, buffer buffer_recorte_m).
Modo "completo": un bloque POR HILERA (zona de búsqueda de utils.tramo_s: con hilera.extension = "nube", la línea del
archivo vectorial ampliada margen_busqueda_m por cada extremo; buffer buffer_recorte_m) para que
03-11 trabajen hilera a hilera con memoria acotada (QA-03). Una sola pasada por el LAS; cada bloque se
acumula en RAM (~0,5 GB por hilera) y se escribe en salidas\\02_recorte[\\hilera_<id>]\\recorte.las.
Después se quita el ruido aéreo/bajo con un suelo aproximado (mínimo por celda) y 'margen_z_m'.
NO se submuestrea.
"""
import utils  # noqa: F401
import os
import numpy as np
import laspy
import shapely

from utils import (cargar_config, dir_salida, Paso, cargar_filas, guardar_json, cfg_bloque, lista_bloques, valores_atributos)
from trayectoria import cargar_trayectoria
from geometria import poligono_recorte, ventanas_tiempo

NOMBRE = "02_recorte"


def suelo_aproximado(x, y, z, celda=1.0):
    """Mínimo de z por celda -> función de consulta (vecino más próximo de celdas con datos)."""
    from scipy.spatial import cKDTree
    ix, iy = np.floor(x / celda).astype(np.int64), np.floor(y / celda).astype(np.int64)
    k = ix * 1_000_003 + iy
    o = np.argsort(k)
    ks, zs = k[o], z[o]
    ini = np.r_[0, np.nonzero(np.diff(ks))[0] + 1]
    zmin = np.minimum.reduceat(zs, ini)
    cx, cy = (ix[o][ini] + .5) * celda, (iy[o][ini] + .5) * celda
    tr = cKDTree(np.c_[cx, cy])
    return lambda xq, yq: zmin[tr.query(np.c_[xq, yq])[1]]


def main(cfg=None):
    cfg = cfg or cargar_config()
    bloques = lista_bloques(cfg)
    cfgs = [cfg_bloque(cfg, b) for b in bloques]
    tray, _ = cargar_trayectoria(cfg, cache_dir=dir_salida(cfg, "01_qc"))
    with Paso(cfg, NOMBRE, entrada=cfg["rutas"]["las_entrada"],
              parametros=dict(cfg["recorte"], bloques=[str(b) for b in bloques])) as p:
        polis, filas_b = [], []
        for c in cfgs:
            f = cargar_filas(c)
            filas_b.append(f)
            polis.append(poligono_recorte(c, f))
            shapely.prepare(polis[-1])
        las_in = laspy.open(cfg["rutas"]["las_entrada"])
        cab = las_in.header
        n_tot = cab.point_count
        p.n_antes = n_tot
        acum = [[] for _ in bloques]
        n_leidos = 0
        for tr in las_in.chunk_iterator(cfg["lectura"]["tam_trozo_puntos"]):
            n_leidos += len(tr)
            x, y = np.asarray(tr.x), np.asarray(tr.y)
            for i, poli in enumerate(polis):
                bx0, by0, bx1, by1 = poli.bounds
                m = (x >= bx0) & (x <= bx1) & (y >= by0) & (y <= by1)
                if m.any():
                    idx = np.nonzero(m)[0]
                    dentro = shapely.contains_xy(poli, x[idx], y[idx])
                    if dentro.any():
                        acum[i].append(tr.array[idx[dentro]].copy())
            if (n_leidos // cfg["lectura"]["tam_trozo_puntos"]) % 5 == 0:
                p.info("leidos %d / %d M", n_leidos // 10 ** 6, n_tot // 10 ** 6)
        las_in.close()
        total = 0
        metas = []
        for c, f, poli, trozos in zip(cfgs, filas_b, polis, acum):
            out = dir_salida(c, NOMBRE)
            arr = np.concatenate(trozos)
            trozos.clear()
            pts = laspy.ScaleAwarePointRecord(arr, cab.point_format, cab.scales, cab.offsets)
            x, y, z = np.asarray(pts.x), np.asarray(pts.y), np.asarray(pts.z)
            n_poly = len(z)
            hz = z - suelo_aproximado(x, y, z)(x, y)
            lo, hi = cfg["recorte"]["margen_z_m"]
            ok = (hz >= lo) & (hz <= hi)
            las_out = laspy.LasData(cab)
            las_out.points = laspy.PackedPointRecord(arr[ok], cab.point_format)
            las_out.header.scales, las_out.header.offsets = cab.scales, cab.offsets
            ruta_las = os.path.join(out, "recorte.las")
            las_out.write(ruta_las)
            gt = np.asarray(las_out.gps_time)
            vent, pasos = ventanas_tiempo(c, tray, f)
            meta = dict(bloque=c.get("_bloque"), n_en_poligono=n_poly, n_final=int(ok.sum()),
                        bbox=[float(x[ok].min()), float(y[ok].min()), float(x[ok].max()), float(y[ok].max())],
                        gps_time=[float(gt.min()), float(gt.max())], ventanas_tiempo=vent, pasos=pasos,
                        filas=[dict(id_fila=r["id_fila"], **valores_atributos(c, r), largo=float(r.geometry.length)) for _, r in f.iterrows()],
                        poligono_wkt=poli.wkt)
            guardar_json(meta, os.path.join(out, "recorte_meta.json"))
            metas.append(meta)
            total += int(ok.sum())
            p.info("bloque %s: %d puntos en polígono, %d tras margen z, %d ventanas de tiempo", c.get("_bloque", "prueba"), n_poly, int(ok.sum()), len(vent))
            del arr, pts, las_out
        p.n_despues = total
        p.salida = dir_salida(cfg, NOMBRE)
    return metas


if __name__ == "__main__":
    main()
