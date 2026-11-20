# -*- coding: utf-8 -*-
"""06_normalizacion: altura sobre el suelo (HAG) con PDAL filters.hag_dem y el DTM de 05.
Se conservan TODOS los puntos (también los < 0,3 m, que sirven de fondo de los rayos); el umbral
de 0,3 m se aplica después, en 08/09, a MH, MW y la envolvente."""
import utils  # noqa: F401
import os
import numpy as np
from utils import cargar_config, dir_salida, Paso, run_pdal

NOMBRE = "06_normalizacion"


def main(cfg=None):
    cfg = cfg or cargar_config()
    ent = os.path.join(dir_salida(cfg, "05_dtm"), "recorte_suelo_calles.las")   # suelo reclasificado por 05 (hierba fuera)
    dtm = os.path.join(dir_salida(cfg, "05_dtm"), "dtm.tif")
    out = dir_salida(cfg, NOMBRE)
    with Paso(cfg, NOMBRE, entrada=ent + " + " + dtm, parametros=cfg["normalizacion"]) as p:
        sal = os.path.join(out, "recorte_hag.las")
        etapas = [{"type": "readers.las", "filename": ent},
                  {"type": "filters.hag_dem", "raster": dtm, "zero_ground": True},
                  {"type": "writers.las", "filename": sal, "forward": "all", "minor_version": 4,
                   "dataformat_id": 7, "extra_dims": "all", "a_srs": "EPSG:%d" % cfg["crs"]["epsg"]}]
        n, pipe = run_pdal(etapas, p.log)
        a = pipe.arrays[0]
        h = a["HeightAboveGround"]
        p.n_antes = n
        p.n_despues = n
        u = cfg["normalizacion"]["umbral_altura_m"]
        p.extra["HAG_min_max"] = (float(h.min()), float(h.max()))
        p.extra["puntos_HAG>=umbral"] = "%d (%.1f %%)" % ((h >= u).sum(), 100 * (h >= u).mean())
        p.extra["HAG_p50_p99"] = (float(np.percentile(h, 50)), float(np.percentile(h, 99)))
        p.salida = sal
    return sal


if __name__ == "__main__":
    main()
