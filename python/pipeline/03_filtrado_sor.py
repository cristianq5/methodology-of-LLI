# -*- coding: utf-8 -*-
"""03_filtrado_sor: Statistical Outlier Removal (PDAL filters.outlier). Los puntos de ruido se
MARCAN (clase 7), no se borran, y NO se submuestrea (rompería el recuento de rayos)."""
import utils  # noqa: F401
import os
from utils import cargar_config, dir_salida, Paso, run_pdal

NOMBRE = "03_filtrado_sor"


def main(cfg=None):
    cfg = cfg or cargar_config()
    ent = os.path.join(dir_salida(cfg, "02_recorte"), "recorte.las")
    out = dir_salida(cfg, NOMBRE)
    s = cfg["sor"]
    assert s["submuestreo_voxel"] is None, "no se permite submuestreo por vóxel antes del recuento de rayos"
    with Paso(cfg, NOMBRE, entrada=ent, parametros=s) as p:
        sal = os.path.join(out, "recorte_sor.las")
        etapas = [{"type": "readers.las", "filename": ent},
                  {"type": "filters.outlier", "method": s["metodo"], "mean_k": s["mean_k"],
                   "multiplier": s["multiplicador"]},
                  {"type": "writers.las", "filename": sal, "forward": "all", "minor_version": 4,
                   "dataformat_id": 7, "extra_dims": "all", "a_srs": "EPSG:%d" % cfg["crs"]["epsg"]}]
        n, pipe = run_pdal(etapas, p.log)
        cls = pipe.arrays[0]["Classification"]
        p.n_antes = n
        p.n_despues = int((cls != 7).sum())
        p.extra["marcados_ruido"] = int((cls == 7).sum())
        p.salida = sal
    return sal


if __name__ == "__main__":
    main()
