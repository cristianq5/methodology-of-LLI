# -*- coding: utf-8 -*-
"""07_chm: modelo de altura de copa (máximo de HAG por celda de 0,1 m, sin ruido, HAG >= umbral)."""
import utils  # noqa: F401
import os
import numpy as np
import rasterio
from utils import cargar_config, dir_salida, Paso, run_pdal, png_matriz

NOMBRE = "07_chm"


def main(cfg=None):
    cfg = cfg or cargar_config()
    ent = os.path.join(dir_salida(cfg, "06_normalizacion"), "recorte_hag.las")
    out = dir_salida(cfg, NOMBRE)
    c = cfg["chm"]
    u = cfg["normalizacion"]["umbral_altura_m"]
    with Paso(cfg, NOMBRE, entrada=ent, parametros=dict(c, umbral=u)) as p:
        sal = os.path.join(out, "chm.tif")
        etapas = [{"type": "readers.las", "filename": ent},
                  {"type": "filters.range", "limits": "Classification![7:7],HeightAboveGround[%s:]" % u},
                  {"type": "writers.gdal", "filename": sal, "dimension": "HeightAboveGround",
                   "resolution": c["resolucion_m"], "output_type": c["estadistico"], "nodata": -9999,
                   "override_srs": "EPSG:%d" % cfg["crs"]["epsg"], "gdaldriver": "GTiff"}]
        n, _ = run_pdal(etapas, p.log)
        p.n_antes = p.n_despues = n
        with rasterio.open(sal) as r:
            a = r.read(1).astype(float)
            ext = (r.bounds.left, r.bounds.right, r.bounds.bottom, r.bounds.top)
        a[a == -9999] = np.nan
        p.extra["CHM_max_m"] = float(np.nanmax(a))
        png_matriz(a, os.path.join(out, "chm.png"), "CHM (m, 0,1 m)", extent=ext, cmap="viridis", cbar="HAG máx (m)")
        p.salida = sal
    return sal


if __name__ == "__main__":
    main()
