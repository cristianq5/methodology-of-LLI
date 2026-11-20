# -*- coding: utf-8 -*-
"""04_suelo: clasificación de suelo (CSF o SMRF de PDAL, según config) ignorando el ruido (clase 7).
El método no lo fija la ficha [SUP]. Suelo = clase 2."""
import utils  # noqa: F401
import os
from utils import cargar_config, dir_salida, Paso, run_pdal

NOMBRE = "04_suelo"


def main(cfg=None):
    cfg = cfg or cargar_config()
    ent = os.path.join(dir_salida(cfg, "03_filtrado_sor"), "recorte_sor.las")
    out = dir_salida(cfg, NOMBRE)
    s = cfg["suelo"]
    with Paso(cfg, NOMBRE, entrada=ent, parametros=s) as p:
        if s["metodo"] == "csf":
            c = s["csf"]
            filt = {"type": "filters.csf", "resolution": c["resolucion"], "rigidness": c["rigidez"],
                    "threshold": c["umbral"], "step": c["paso_tiempo"], "iterations": c["iteraciones"],
                    "smooth": bool(c["slope_smooth"]), "ignore": "Classification[7:7]"}
        else:
            c = s["smrf"]
            filt = {"type": "filters.smrf", "cell": c["celda"], "slope": c["pendiente"], "window": c["ventana"],
                    "threshold": c["umbral"], "scalar": c["escalar"], "ignore": "Classification[7:7]"}
        sal = os.path.join(out, "recorte_suelo.las")
        etapas = [{"type": "readers.las", "filename": ent}, filt,
                  {"type": "writers.las", "filename": sal, "forward": "all", "minor_version": 4,
                   "dataformat_id": 7, "extra_dims": "all", "a_srs": "EPSG:%d" % cfg["crs"]["epsg"]}]
        n, pipe = run_pdal(etapas, p.log)
        cls = pipe.arrays[0]["Classification"]
        p.n_antes = n
        p.n_despues = int((cls == 2).sum())
        p.extra["puntos_suelo"] = "%d (%.1f %%)" % (p.n_despues, 100 * p.n_despues / n)
        p.salida = sal
    return sal


if __name__ == "__main__":
    main()
