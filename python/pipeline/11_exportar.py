# -*- coding: utf-8 -*-
"""11_exportar: GPKG (hileras, tramos, rebanadas), CSV (incluido el perfil de anchura por altura), GeoTIFF (DTM, CHM, LLI por rebanada),
LAS/LAZ con dimensiones extra (HAG, id_hilera, s, t, lado) y copia de la configuración usada."""
import utils  # noqa: F401
import os
import shutil
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.features import rasterize
from shapely.geometry import LineString, Polygon

from utils import cargar_config, dir_salida, Paso, MarcoHilera, cargar_json, copiar_config

NOMBRE = "11_exportar"


def main(cfg=None):
    cfg = cfg or cargar_config()
    out = dir_salida(cfg, NOMBRE)
    e = cfg["exportar"]
    epsg = cfg["crs"]["epsg"]
    d10, d09, d08 = (dir_salida(cfg, n) for n in ("10_lli", "09_metricas", "08_segmentacion"))
    with Paso(cfg, NOMBRE, entrada=d10, parametros=e) as p:
        H = pd.read_csv(os.path.join(d10, "lli_hileras.csv"), dtype={"id_fila": str})
        T = pd.read_csv(os.path.join(d10, "lli_tramos.csv"), dtype={"id_fila": str})
        Rb = pd.read_csv(os.path.join(d10, "lli_rebanadas.csv"), dtype={"id_fila": str})
        marcos = cargar_json(os.path.join(d08, "marcos.json"))["marcos"]
        p.n_antes = len(Rb)
        # --- geometrías
        def marco(k):
            return MarcoHilera.de_dict(marcos[str(k)]["marco"])
        gh = []
        for _, r in H.iterrows():
            m = marco(r.id_fila)
            s0, s1 = marcos[str(r.id_fila)]["s_rango"]
            gh.append(LineString([m.xy(s0, 0.0), m.xy(s1, 0.0)]))
        GH = gpd.GeoDataFrame(H.assign(id_fila=H.id_fila.astype(str)), geometry=gh, crs=epsg)
        gt = [LineString([marco(r.id_fila).xy(r.s_ini, 0.0), marco(r.id_fila).xy(r.s_fin, 0.0)]) for _, r in T.iterrows()]
        GT = gpd.GeoDataFrame(T.assign(id_fila=T.id_fila.astype(str)), geometry=gt, crs=epsg)
        polis = []
        for _, r in Rb.iterrows():
            m = marco(r.id_fila)
            w = float(r.MW_max) / 2 if np.isfinite(r.MW_max) else 0.05
            polis.append(Polygon([m.xy(r.s_ini, -w), m.xy(r.s_ini + 0.1, -w), m.xy(r.s_ini + 0.1, w), m.xy(r.s_ini, w)]))
        GR = gpd.GeoDataFrame(Rb.assign(id_fila=Rb.id_fila.astype(str)), geometry=polis, crs=epsg)
        ruta_gpkg = os.path.join(out, "lli_resultados.gpkg")
        if os.path.exists(ruta_gpkg):
            os.remove(ruta_gpkg)
        if "gpkg" in e["formatos"]:
            GH.to_file(ruta_gpkg, layer="hileras", driver="GPKG")
            GT.to_file(ruta_gpkg, layer="tramos_1_5m", driver="GPKG")
            GR.to_file(ruta_gpkg, layer="rebanadas", driver="GPKG")
        # --- CSV
        if "csv" in e["formatos"]:
            for nom in ("lli_hileras.csv", "lli_tramos.csv", "lli_rebanadas.csv", "lli_resumen.json"):
                shutil.copy2(os.path.join(d10, nom), out)
            for nom in ("extrinseca.json", "comparacion_P.json", "metricas_rebanadas.csv", "perfil_anchura.csv"):
                shutil.copy2(os.path.join(d09, nom), out)
        # --- GeoTIFF
        if "geotiff" in e["formatos"]:
            shutil.copy2(os.path.join(dir_salida(cfg, "05_dtm"), "dtm.tif"), os.path.join(out, "dtm.tif"))
            shutil.copy2(os.path.join(dir_salida(cfg, "07_chm"), "chm.tif"), os.path.join(out, "chm.tif"))
            with rasterio.open(os.path.join(out, "chm.tif")) as r:
                tr = rasterio.transform.from_origin(r.bounds.left, r.bounds.top, 0.1, 0.1)
                shape = (int(np.ceil((r.bounds.top - r.bounds.bottom) / 0.1)), int(np.ceil((r.bounds.right - r.bounds.left) / 0.1)))
            ok = GR["LLI_estricta"].notna()
            ras = rasterize(zip(GR.geometry[ok], GR["LLI_estricta"][ok]), out_shape=shape, transform=tr, fill=-9999, dtype="float32")
            with rasterio.open(os.path.join(out, "lli_estricta_rebanadas.tif"), "w", driver="GTiff", height=shape[0], width=shape[1],
                               count=1, dtype="float32", crs=f"EPSG:{epsg}", transform=tr, nodata=-9999) as w:
                w.write(ras, 1)
        # --- nube con dimensiones extra
        if e["las_con_dimensiones_extra"]:
            src = os.path.join(d08, "segmentacion.las")
            if e["comprimir_laz"]:
                import laspy
                las = laspy.read(src)
                las.write(os.path.join(out, "nube_hileras.laz"), do_compress=True)
            else:
                shutil.copy2(src, os.path.join(out, "nube_hileras.las"))
        if e["guardar_copia_config"]:
            copiar_config(cfg, out)
        p.n_despues = len(GR)
        p.extra["capas_gpkg"] = "hileras(%d), tramos_1_5m(%d), rebanadas(%d)" % (len(GH), len(GT), len(GR))
        p.salida = out
    return out


if __name__ == "__main__":
    main()
