# -*- coding: utf-8 -*-
"""05_dtm: ráster de terreno y QA de la altura del sensor sobre el DTM.

QA-11 (confirmado por el usuario: bajo la hilera hay HIERBA, no un caballón): el CSF clasifica como suelo la hierba
de debajo del seto (~0,1 m por encima del suelo de las calles). El DTM se interpola SOLO con el suelo de las calles
(se excluye una banda de +-suelo.excluir_bajo_hilera_m alrededor del eje GPKG de cada hilera) y el hueco bajo el seto se
rellena por interpolación. Los puntos de suelo que quedan más de suelo.rechazo_sobre_calle_m (0,12 m) por encima de ese DTM
se reclasifican como no-suelo (clase 1, hierba) en recorte_suelo_calles.las, que es la entrada de 06.
Salida: dtm.tif (IDW 0,25 m + relleno), recorte_suelo_calles.las, qa_suelo.json, qa_altura_sensor_sobre_dtm.json.
"""
import utils  # noqa: F401
import os
import numpy as np
import laspy
import rasterio
from rasterio.fill import fillnodata
from scipy.ndimage import distance_transform_edt

from utils import (cargar_config, dir_salida, Paso, run_pdal, guardar_json, png_matriz, cargar_filas, MarcoHilera,
                   escribir_las_desde)
from trayectoria import cargar_trayectoria
from geometria import ventanas_tiempo

NOMBRE = "05_dtm"


def muestrear_dtm(ruta_dtm, x, y):
    """Valor del DTM (celda) en puntos x,y; NaN fuera o sin dato."""
    with rasterio.open(ruta_dtm) as r:
        a = r.read(1)
        nd = r.nodata
        rows, cols = rasterio.transform.rowcol(r.transform, x, y)
    rows, cols = np.asarray(rows), np.asarray(cols)
    ok = (rows >= 0) & (rows < a.shape[0]) & (cols >= 0) & (cols < a.shape[1])
    v = np.full(len(x), np.nan)
    v[ok] = a[rows[ok], cols[ok]]
    v[v == nd] = np.nan
    return v


def _idw_fill(ent_las, bruto, out_tif, res, d, epsg, log):
    etapas = [{"type": "readers.las", "filename": ent_las},
              {"type": "writers.gdal", "filename": bruto, "resolution": res, "output_type": "idw",
               "window_size": d["idw_window"], "nodata": -9999, "override_srs": "EPSG:%d" % epsg, "gdaldriver": "GTiff"}]
    n, _ = run_pdal(etapas, log)
    with rasterio.open(bruto) as r:
        a = r.read(1).astype("float32")
        perfil = r.profile
    vacio0 = a == -9999
    base = np.where(vacio0, 0, a).astype("float32")
    rel = fillnodata(base, mask=(~vacio0).astype("uint8"),
                     max_search_distance=max(1, int(round(d["rellenar_huecos_dist_max_m"] / res))))
    sin_dato = vacio0 & (rel == 0)
    if sin_dato.any():       # relleno final por vecino más próximo (suelo suave)
        _, ii = distance_transform_edt(sin_dato, return_indices=True)
        rel = rel[tuple(ii)]
    perfil.update(nodata=-9999, dtype="float32")
    with rasterio.open(out_tif, "w", **perfil) as w:
        w.write(rel.astype("float32"), 1)
    return n, rel, float(100 * vacio0.mean())


def main(cfg=None):
    cfg = cfg or cargar_config()
    ent = os.path.join(dir_salida(cfg, "04_suelo"), "recorte_suelo.las")
    out = dir_salida(cfg, NOMBRE)
    d = cfg["dtm"]
    sc = cfg["suelo"]
    with Paso(cfg, NOMBRE, entrada=ent, parametros=dict(d, banda_excluida_m=sc["excluir_bajo_hilera_m"],
                                                        rechazo_m=sc["rechazo_sobre_calle_m"])) as p:
        res = d["resolucion_m"]
        las = laspy.read(ent)
        x, y, z = np.asarray(las.x), np.asarray(las.y), np.asarray(las.z)
        cls = np.asarray(las.classification)
        suelo = cls == 2
        p.n_antes = int(suelo.sum())
        # --- banda bajo cada hilera (eje del GPKG de todas las hileras)
        banda = np.zeros(len(x), bool)
        for _, f in cargar_filas(cfg, todas=True).iterrows():
            m = MarcoHilera.de_linea(f.geometry)
            s, t = m.st(np.c_[x, y])
            banda |= (np.abs(t) < sc["excluir_bajo_hilera_m"]) & (s > -1.0) & (s < m.largo + 1.0)
        calle = suelo & ~banda
        p.info("suelo CSF: %d puntos; en la banda bajo la hilera: %d; de calle: %d", int(suelo.sum()), int((suelo & banda).sum()), int(calle.sum()))
        # --- DTM con el suelo de las calles
        tmp = os.path.join(out, "_suelo_calles.las")
        escribir_las_desde(las, calle, tmp)
        ruta_dtm = os.path.join(out, "dtm.tif")
        n, rel, pct_vacio = _idw_fill(tmp, os.path.join(out, "dtm_bruto.tif"), ruta_dtm, res, d, cfg["crs"]["epsg"], p.log)
        os.remove(tmp)
        p.extra["celdas_sin_dato_inicial_pct"] = round(pct_vacio, 1)
        # --- escalón bajo la hilera antes/después y reclasificación de la hierba
        zd = muestrear_dtm(ruta_dtm, x, y)
        dz = z - zd
        en_banda = suelo & banda & np.isfinite(dz)
        esc_antes = float(np.median(dz[en_banda])) if en_banda.any() else float("nan")
        recl = en_banda & (dz > sc["rechazo_sobre_calle_m"])
        cls2 = cls.copy()
        cls2[recl] = 1
        suelo2 = cls2 == 2
        en_banda2 = suelo2 & banda & np.isfinite(dz)
        esc_despues = float(np.median(dz[en_banda2])) if en_banda2.any() else float("nan")
        p95_antes = float(np.percentile(dz[en_banda], 95)) if en_banda.any() else float("nan")
        las.classification = cls2
        ruta_las = os.path.join(out, "recorte_suelo_calles.las")
        las.write(ruta_las)
        qa_suelo = dict(puntos_suelo_csf=int(suelo.sum()), puntos_suelo_bajo_hilera=int((suelo & banda).sum()),
                        reclasificados_hierba=int(recl.sum()), pct_reclasificados_de_banda=round(100 * recl.sum() / max(1, (suelo & banda).sum()), 1),
                        escalon_mediano_bajo_hilera_antes_m=round(esc_antes, 3), escalon_p95_antes_m=round(p95_antes, 3),
                        escalon_mediano_bajo_hilera_despues_m=round(esc_despues, 3), banda_m=sc["excluir_bajo_hilera_m"],
                        rechazo_m=sc["rechazo_sobre_calle_m"])
        guardar_json(qa_suelo, os.path.join(out, "qa_suelo.json"))
        p.info("QA suelo: %s", qa_suelo)
        p.n_despues = int(suelo2.sum())
        with rasterio.open(ruta_dtm) as r_:
            ext_ = (r_.bounds.left, r_.bounds.right, r_.bounds.bottom, r_.bounds.top)
        png_matriz(np.where(rel == 0, np.nan, rel), os.path.join(out, "dtm.png"), "DTM (m) desde el suelo de las calles - UTM 30N",
                   extent=ext_, cmap="terrain", cbar="z (m)")
        # --- QA altura del sensor sobre el DTM (pasos del bloque)
        tray, _ = cargar_trayectoria(cfg, cache_dir=dir_salida(cfg, "01_qc"))
        _, pasos = ventanas_tiempo(cfg, tray, cargar_filas(cfg))
        qa = []
        for pa in pasos:
            mm = (tray.t >= pa["t0"]) & (tray.t <= pa["t1"])
            h = tray.p[mm, 2] - muestrear_dtm(ruta_dtm, tray.p[mm, 0], tray.p[mm, 1])
            h = h[np.isfinite(h)]
            qa.append(dict(fila=pa["fila"], lado=pa["lado"], t_med_m=round(pa["t_med"], 2), n=int(len(h)),
                           h_p5=round(float(np.percentile(h, 5)), 2) if len(h) else None,
                           h_med=round(float(np.median(h)), 2) if len(h) else None,
                           h_p95=round(float(np.percentile(h, 95)), 2) if len(h) else None))
        lo, hi = cfg["qc"]["altura_sensor_sobre_dtm_m"]
        meds = [q["h_med"] for q in qa if q["h_med"] is not None]
        res_qa = dict(pasos=qa, rango_esperado=[lo, hi], mediana_global=float(np.median(meds)) if meds else None,
                      dentro_rango=bool(meds and lo <= np.median(meds) <= hi))
        guardar_json(res_qa, os.path.join(out, "qa_altura_sensor_sobre_dtm.json"))
        p.salida = ruta_dtm
    return ruta_dtm


if __name__ == "__main__":
    main()
