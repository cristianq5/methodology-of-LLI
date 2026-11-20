# -*- coding: utf-8 -*-
"""Test sintético del archivo vectorial genérico (utils.cargar_filas), sin datos reales.

Caso: un Shapefile en EPSG:4326 (no en el CRS de trabajo) con 4 líneas y campos propios ('linea', 'cultivar', 'trat').
- los ids se repiten entre tratamientos -> sin id_compuesto debe fallar con un mensaje claro;
- con filtro {trat: "Secano"} (mayúsculas y espacios en el dato) quedan 2 hileras;
- con id_compuesto quedan 4 ids únicos y válidos para nombres de carpeta;
- la errata 'arbosona' se corrige en memoria y el archivo vectorial se reproyecta a EPSG:25830.
"""
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import utils  # noqa: E402  (primero: entorno)
import geopandas as gpd  # noqa: E402
from shapely.geometry import LineString  # noqa: E402


def _cfg(ruta, **croquis):
    cfg = utils.cargar_parametros()
    cq = dict(capa=None, campo_id="linea", atributos=["cultivar", "trat"], filtro={}, normalizar_texto=True,
              reemplazos={"cultivar": {"arbosona": "arbosana"}}, id_compuesto=False)
    cq.update(croquis)
    cfg["croquis"] = cq
    cfg["rutas"] = dict(filas_gpkg=ruta, filas_capa=None)
    cfg["crs"] = dict(cfg["crs"], epsg=25830)
    cfg["recorte"] = dict(cfg["recorte"], modo="completo")
    return cfg


def test_croquis():
    d = tempfile.mkdtemp()
    ruta = os.path.join(d, "croquis.shp")
    lon0, lat0 = -4.43, 37.47
    lineas = [LineString([(lon0 + i * 5e-5, lat0), (lon0 + i * 5e-5, lat0 + 4e-4)]) for i in range(4)]
    gpd.GeoDataFrame(dict(linea=["A1", "A2", "A1", "A2"], cultivar=[" Arbosona", "picual", "picual", "arbequina"],
                          trat=["Secano ", "secano", "riego", "riego"]), geometry=lineas, crs=4326).to_file(ruta)
    try:
        utils.cargar_filas(_cfg(ruta), todas=True)
        raise AssertionError("debía fallar por id repetido")
    except ValueError as e:
        assert "se repite" in str(e), e
    f = utils.cargar_filas(_cfg(ruta, filtro={"trat": "SECANO"}), todas=True)
    assert list(f.id_fila) == ["A1", "A2"], list(f.id_fila)
    assert list(f.cultivar) == ["arbosana", "picual"], list(f.cultivar)
    assert f.crs.to_epsg() == 25830 and 30 < f.geometry.iloc[0].length < 60, (f.crs, f.geometry.iloc[0].length)
    g = utils.cargar_filas(_cfg(ruta, id_compuesto=True), todas=True)
    assert g.id_fila.is_unique and g.id_fila.str.match(r"^[0-9A-Za-z_.\-]+$").all(), list(g.id_fila)
    assert utils.valores_atributos(_cfg(ruta), f.iloc[0]) == {"cultivar": "arbosana", "trat": "secano"}
    try:
        utils.cargar_filas(_cfg(ruta, atributos=["no_existe"], filtro={"trat": "riego"}), todas=True)
        raise AssertionError("debía fallar por campo inexistente")
    except KeyError as e:
        assert "no_existe" in str(e)
    print("OK croquis genérico: filtro %s, ids compuestos %s, largo %.1f m" % (list(f.id_fila), list(g.id_fila), f.geometry.iloc[0].length))


if __name__ == "__main__":
    test_croquis()
