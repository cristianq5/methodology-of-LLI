# -*- coding: utf-8 -*-
"""Utilidades comunes del flujo LLI (configuración, rutas, archivo vectorial de hileras, logs, LAS).

IMPORTAR ESTE MÓDULO EL PRIMERO: corrige el PATH de las DLL del entorno conda y define
GDAL_DATA / PROJ_DATA. Sin ello, en esta máquina numpy.linalg.svd y GDAL pueden cerrar el
proceso sin mensaje (DLL de LAPACK/PROJ no encontradas al lanzar python.exe sin activar el entorno).
"""
import os
import sys

# --- Entorno (antes de importar numpy/gdal/pdal) ---------------------------------
_PREFIJO = sys.prefix
for _sub in ("", r"Library\mingw-w64\bin", r"Library\usr\bin", r"Library\bin", "Scripts"):
    _d = os.path.join(_PREFIJO, _sub)
    if os.path.isdir(_d) and _d not in os.environ.get("PATH", ""):
        os.environ["PATH"] = _d + os.pathsep + os.environ.get("PATH", "")
os.environ.setdefault("GDAL_DATA", os.path.join(_PREFIJO, "Library", "share", "gdal"))
os.environ.setdefault("PROJ_DATA", os.path.join(_PREFIJO, "Library", "share", "proj"))
os.environ.setdefault("PROJ_LIB", os.environ["PROJ_DATA"])

import json
import logging
import time
import shutil
import numpy as np
import yaml

DIR_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
# Raíz del repositorio: <repo>/python/pipeline -> <repo>
RAIZ = os.path.dirname(os.path.dirname(DIR_SCRIPTS))
DIR_CONFIG = os.path.join(RAIZ, "config")
RUTA_PARAMETROS = os.path.join(DIR_CONFIG, "parametros.yaml")   # parámetros del método (comunes a todos los ensayos)
VAR_ENTORNO_ENSAYO = "LLI_ENSAYO"     # ficha del ensayo por defecto (para lanzar un paso suelto: python 05_dtm.py)

# Ficheros de la captura del CHCNAV RS30 que usa el pipeline y dónde están dentro de la carpeta de brutos.
# (patrón relativo a la carpeta de brutos, obligatorio)
ARCHIVOS_BRUTOS = {
    "trayectoria_posi": ("SLAM_DATA/*.PosI", True),                                  # 200 Hz, lat/lon/h + RPY
    "trayectoria_hf": ("POST_SLAM_DATA/result_trajectory3d_high_freq.txt", True),    # post-SLAM, t x y z qx qy qz qw
    "rtk_gga": ("GPS/Rtk/rtk.gga", False),                                           # solo control de calidad
    "pcap": ("SCAN/Scanner1/*.pcap", True),                                          # disparos del Hesai XT32 (Bt, Bi)
    "calib_ep": ("SCAN/Scanner1/*.EP", False),                                       # calibración del escáner (tabla de elevaciones; modo global_bloque)
}


# --- Configuración y rutas -----------------------------------------------------
def _fusionar(base, extra):
    """Copia de 'base' con los valores de 'extra' encima (diccionarios anidados, recursivo)."""
    out = dict(base)
    for k, v in (extra or {}).items():
        out[k] = _fusionar(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def _ruta(base_dir, r):
    """Ruta absoluta: las relativas se resuelven desde la carpeta de la ficha del ensayo."""
    if r in (None, ""):
        return None
    r = os.path.expanduser(str(r))
    return os.path.normpath(r if os.path.isabs(r) else os.path.join(base_dir, r))


def _buscar_unico(carpeta, patron, nombre, obligatorio):
    import glob
    c = sorted(glob.glob(os.path.join(carpeta, patron)))
    if len(c) == 1:
        return os.path.normpath(c[0])
    if not c:
        if obligatorio:
            raise FileNotFoundError("No se encuentra '%s' (%s) en %s. Indícalo en entradas.brutos_archivos.%s"
                                    % (nombre, patron, carpeta, nombre))
        return None
    raise ValueError("Hay %d candidatos para '%s' en %s: %s. Indica uno en entradas.brutos_archivos.%s"
                     % (len(c), nombre, carpeta, [os.path.basename(x) for x in c], nombre))


def _epsg_las(ruta):
    import laspy
    with laspy.open(ruta) as f:
        crs = f.header.parse_crs()
    return crs.to_epsg() if crs is not None else None


def cargar_parametros(ruta=None):
    """Solo los parámetros del método (parametros.yaml), sin ensayo. Lo usan los tests sintéticos."""
    with open(ruta or RUTA_PARAMETROS, encoding="utf-8") as f:
        return yaml.safe_load(f)


def cargar_config(ruta_ensayo=None, ruta_parametros=None):
    """Configuración efectiva = parametros.yaml + ficha del ensayo (entradas, archivo vectorial, CRS, modo, salidas)
    + 'ajustes' de la ficha. Devuelve un único diccionario con la estructura que usan los pasos 01-11.
    Comprueba que existen todas las entradas antes de empezar."""
    ruta_ensayo = ruta_ensayo or os.environ.get(VAR_ENTORNO_ENSAYO)
    if not ruta_ensayo:
        raise ValueError("Falta la ficha del ensayo: run_pipeline.py --ensayo <ficha.yaml> "
                         "(o la variable de entorno %s). Plantilla: config/ensayo_plantilla.yaml" % VAR_ENTORNO_ENSAYO)
    ruta_ensayo = os.path.abspath(ruta_ensayo)
    with open(ruta_ensayo, encoding="utf-8") as f:
        e = yaml.safe_load(f) or {}
    for sec in ("ensayo", "entradas", "croquis"):
        if sec not in e:
            raise ValueError("La ficha %s no tiene la sección '%s' (ver config/ensayo_plantilla.yaml)" % (ruta_ensayo, sec))
    base = os.path.dirname(ruta_ensayo)
    cfg = _fusionar(cargar_parametros(_ruta(base, e.get("parametros")) or ruta_parametros), e.get("ajustes") or {})

    # --- entradas: LAS + ficheros de la carpeta de brutos del RS30
    ent = e["entradas"]
    rutas = dict(las_entrada=_ruta(base, ent.get("las")), brutos=_ruta(base, ent.get("brutos")))
    elegidos = ent.get("brutos_archivos") or {}
    for nombre, (patron, oblig) in ARCHIVOS_BRUTOS.items():
        r = _ruta(base, elegidos.get(nombre))
        if r is None:
            if not rutas["brutos"]:
                raise ValueError("Falta entradas.brutos (o entradas.brutos_archivos.%s)" % nombre)
            r = _buscar_unico(rutas["brutos"], patron, nombre, oblig)
        rutas[nombre] = r
    cq = _fusionar(dict(capa=None, campo_id="id_fila", atributos=[], filtro={}, normalizar_texto=True,
                        reemplazos={}, id_compuesto=False), e["croquis"])
    cq["ruta"] = _ruta(base, cq.get("ruta"))
    rutas["filas_gpkg"], rutas["filas_capa"] = cq["ruta"], cq.get("capa")
    for k in ("las_entrada", "filas_gpkg", "trayectoria_posi", "trayectoria_hf", "pcap"):
        if not rutas.get(k) or not os.path.exists(rutas[k]):
            raise FileNotFoundError("No existe la entrada '%s': %s" % (k, rutas.get(k)))
    sal = e.get("salidas") or {}
    rutas["salida"] = _ruta(base, sal.get("carpeta") or "salidas")
    rutas["logs"] = _ruta(base, sal.get("logs")) or os.path.join(rutas["salida"], "logs")
    cfg["pcap"] = dict(cfg["pcap"], ruta=rutas.pop("pcap"))
    cfg["rutas"], cfg["croquis"] = rutas, cq

    # --- CRS de trabajo (proyectado, en metros)
    epsg = (e.get("crs") or {}).get("epsg") or _epsg_las(rutas["las_entrada"])
    if not epsg:
        raise ValueError("El LAS no trae CRS en la cabecera: indica crs.epsg en la ficha del ensayo")
    cfg["crs"] = dict(cfg.get("crs") or {}, epsg=int(epsg))

    # --- modo de ejecución y hileras de prueba
    ej = e.get("ejecucion_ensayo") or {}
    cfg["recorte"] = dict(cfg["recorte"], modo=ej.get("modo", "prueba"),
                          filas_prueba_id=[str(x) for x in (ej.get("hileras_prueba") or [])],
                          tramo_prueba_s_m=list(ej.get("tramo_prueba_s_m") or [0.0, 1e9]))
    if cfg["recorte"]["modo"] not in ("prueba", "completo"):
        raise ValueError("ejecucion_ensayo.modo debe ser 'prueba' o 'completo'")
    cfg["proyecto"] = dict(cfg["proyecto"], nombre=str(e["ensayo"].get("nombre", "ensayo")))
    cfg["ensayo"] = e["ensayo"]
    cfg["_ruta_ensayo"] = ruta_ensayo
    return cfg


PASOS_GLOBALES = ("01_qc",)          # salidas comunes a todos los bloques (trayectoria, QC)


def dir_salida(cfg, nombre):
    """salidas\\NN_nombre (se crea si no existe). Si el cfg lleva '_bloque' (modo completo, un bloque =
    una hilera), los pasos 02-11 escriben en salidas\\NN_nombre\\hilera_<id> (QA-03)."""
    d = os.path.join(cfg["rutas"]["salida"], nombre)
    if cfg.get("_bloque") and nombre not in PASOS_GLOBALES:
        d = os.path.join(d, "hilera_%s" % cfg["_bloque"])
    os.makedirs(d, exist_ok=True)
    return d


def cfg_bloque(cfg, bloque):
    """Copia del cfg para un bloque (None = recorte de prueba completo con sus hileras de prueba)."""
    c = dict(cfg)
    if bloque is not None:
        c["_bloque"] = str(bloque)
        c["_filas_bloque"] = [str(bloque)]
    return c


def lista_bloques(cfg):
    """[None] en modo prueba (un único bloque = las hileras de prueba); un bloque por hilera en modo completo."""
    if cfg["recorte"]["modo"] == "prueba":
        return [None]
    ids = list(cargar_filas(cfg, todas=True)["id_fila"])
    if cfg.get("_bloques_sel"):
        ids = [i for i in ids if str(i) in {str(b) for b in cfg["_bloques_sel"]}]
    return ids


def extension_desde_nube(cfg):
    """True si la longitud de la hilera se toma de la nube (hilera.extension = "nube") y no de la línea del archivo vectorial."""
    return cfg["hilera"].get("extension", "gpkg") == "nube" and cfg["recorte"]["modo"] != "prueba"


def tramo_s(cfg, marco):
    """Tramo [s0, s1] de la hilera, en m a lo largo de la línea del archivo vectorial (QA-01).
    Modo prueba: tramo_prueba_s_m (recortado a los extremos de la línea).
    Modo completo con hilera.extension = "gpkg": la línea entera menos recorte_extremos_m en cada cabecera.
    Modo completo con hilera.extension = "nube": ZONA DE BÚSQUEDA = la línea ampliada margen_busqueda_m por cada extremo;
    la longitud real de la hilera la fija 08_segmentacion con la copa de la nube (extension_copa)."""
    ext = cfg["hilera"]["recorte_extremos_m"]
    lo, hi = ext, marco.largo - ext
    if cfg["recorte"]["modo"] == "prueba":
        a, b = cfg["recorte"]["tramo_prueba_s_m"]
        return max(a, lo), min(b, hi)
    if extension_desde_nube(cfg):
        mg = cfg["hilera"]["margen_busqueda_m"]
        return -mg, marco.largo + mg
    return lo, hi


def extension_copa(cfg, s_copa, lim, ref):
    """Inicio y final reales de la hilera a partir de la copa (s de los puntos con HAG >= extension_hag_min_m).

    s_copa: s de esos puntos (eje reajustado); lim = (s0, s1) zona de búsqueda (la línea +- margen_busqueda_m);
    ref = (g0, g1) extremos de la línea del archivo vectorial en el mismo eje (solo referencia).
    Copa: un tramo de extension_celda_m es copa si tiene >= max(extension_min_puntos, extension_frac_min x mediana de
    los tramos con puntos a lo largo de la línea) puntos (las ramas sueltas y la vegetación no llegan).
    1. Dentro de la línea: primer y último tramo con copa (si la copa termina antes que la línea, la hilera termina ahí;
       los huecos internos, p. ej. un olivo que falte, no cortan la hilera).
    2. Fuera de la línea, si la copa llega a su extremo: barrido en pasos de extension_paso_m (0,1, 0,2 ... hasta
       margen_busqueda_m = 1 m); el límite es el último paso con copa pegada a la hilera (el primer paso sin copa corta).
    Después se quita recorte_extremos_m en cada extremo (0 = ninguno). Devuelve (s0, s1, info) con el barrido de cada extremo."""
    h = cfg["hilera"]
    c, paso = h["extension_celda_m"], h["extension_paso_m"]
    M = min(h["margen_busqueda_m"], ref[0] - lim[0], lim[1] - ref[1])
    s_ = np.sort(np.asarray(s_copa, dtype=float))
    contar = lambda a, b: int(np.searchsorted(s_, b, "left") - np.searchsorted(s_, a, "left"))   # puntos en [a, b)
    g0, g1 = float(ref[0]), float(ref[1])
    bordes = np.arange(g0, g1, c)
    cnt = np.array([contar(a, min(a + c, g1)) for a in bordes])
    tipico = float(np.median(cnt[cnt > 0])) if (cnt > 0).any() else 0.0
    umbral = max(h["extension_min_puntos"], h.get("extension_frac_min", 0.0) * tipico)    # por tramo de longitud c
    es_copa = lambda a, b: contar(a, b) >= umbral * (b - a) / c
    occ = np.nonzero([es_copa(a, min(a + c, g1)) for a in bordes])[0]
    if len(occ) == 0:
        return g0, g1, dict(metodo="gpkg (sin copa en la nube)")
    s0, s1 = float(bordes[occ[0]]), float(min(bordes[occ[-1]] + c, g1))
    eps = 1e-6
    barrido = {"inicio": [], "final": []}
    pasos = [round(paso * i, 3) for i in range(1, int(round(M / paso)) + 1)]
    if s1 >= g1 - c + eps:                     # la copa llega al final de la línea: se busca el límite real fuera
        for m in pasos:
            a, b = g1 + m - paso, g1 + m
            ok = es_copa(a, b)
            barrido["final"].append([m, contar(a, b), ok])
            if not ok:
                break
            s1 = b
    if s0 <= g0 + c - eps:                     # idem al inicio
        for m in pasos:
            a, b = g0 - m, g0 - m + paso
            ok = es_copa(a, b)
            barrido["inicio"].append([m, contar(a, b), ok])
            if not ok:
                break
            s0 = a
    ext = h["recorte_extremos_m"]
    info = dict(metodo="nube", s_copa=[round(s0, 2), round(s1, 2)], s_gpkg=[round(g0, 2), round(g1, 2)],
                dif_inicio_m=round(g0 - s0, 2), dif_final_m=round(s1 - g1, 2),
                largo_nube_m=round(s1 - s0, 2), largo_gpkg_m=round(g1 - g0, 2),
                barrido_inicio=barrido["inicio"], barrido_final=barrido["final"],
                toca_limite_busqueda=bool(g0 - s0 >= M - eps or s1 - g1 >= M - eps), recorte_extremos_m=ext,
                puntos_tramo_tipico=round(tipico, 1), umbral_puntos_tramo=round(umbral, 1))
    return s0 + ext, s1 - ext, info


def huella(obj):
    """Hash corto y estable de un objeto serializable (QA-02: validación de cachés)."""
    import hashlib
    t = json.dumps(obj, sort_keys=True, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    return hashlib.sha1(t.encode("utf-8")).hexdigest()[:16]


def huella_fichero(ruta):
    st = os.stat(ruta)
    return "%s|%d|%d" % (os.path.basename(ruta), st.st_size, int(st.st_mtime))


def huella_codigo(*nombres):
    """Hash del texto de módulos .py de python/pipeline (versión del código que generó una caché)."""
    import hashlib
    h = hashlib.sha1()
    for n in nombres:
        with open(os.path.join(DIR_SCRIPTS, n), "rb") as f:
            h.update(f.read())
    return h.hexdigest()[:16]


def dir_logs(cfg):
    d = cfg["rutas"]["logs"]
    os.makedirs(d, exist_ok=True)
    return d


class Paso:
    """Contexto de un paso: escribe en logs\\NN_nombre.log entrada, parámetros,
    nº de puntos antes/después, tiempo y salida.

    Uso:  with Paso(cfg, "02_recorte", entrada=..., parametros={...}) as p:
              ...; p.n_antes = ...; p.n_despues = ...; p.salida = ...
    """

    def __init__(self, cfg, nombre, entrada="", parametros=None):
        self.cfg, self.nombre = cfg, nombre
        self.entrada, self.parametros = entrada, parametros or {}
        self.n_antes = self.n_despues = None
        self.salida = ""
        self.extra = {}
        self.log = logging.getLogger(nombre)

    def __enter__(self):
        self.t0 = time.time()
        suf = ("_hilera_%s" % self.cfg["_bloque"]) if self.cfg.get("_bloque") else ""
        ruta = os.path.join(dir_logs(self.cfg), self.nombre + suf + ".log")
        self.log.setLevel(logging.INFO)
        self.log.propagate = False
        self.log.handlers.clear()
        fh = logging.FileHandler(ruta, mode="w", encoding="utf-8")
        fh.setFormatter(logging.Formatter("%(asctime)s %(message)s"))
        sh = logging.StreamHandler(sys.stdout)
        sh.setFormatter(logging.Formatter("[%(name)s] %(message)s"))
        self.log.addHandler(fh)
        self.log.addHandler(sh)
        self.log.info("INICIO")
        self.log.info("entrada: %s", self.entrada)
        self.log.info("parametros: %s", json.dumps(self.parametros, ensure_ascii=False, default=str))
        return self

    def info(self, *a):
        self.log.info(*a)

    def __exit__(self, tipo, val, tb):
        dt = time.time() - self.t0
        if tipo is not None:
            self.log.error("FALLO: %r", val, exc_info=(tipo, val, tb))
        self.log.info("puntos antes: %s | despues: %s", self.n_antes, self.n_despues)
        for k, v in self.extra.items():
            self.log.info("%s: %s", k, v)
        self.log.info("salida: %s", self.salida)
        self.log.info("tiempo: %.1f s", dt)
        self.log.info("FIN" if tipo is None else "FIN CON ERROR")
        for h in list(self.log.handlers):
            h.close()
            self.log.removeHandler(h)
        return False


def guardar_json(obj, ruta):
    def _c(o):
        if isinstance(o, (np.floating, np.integer)):
            return o.item()
        if isinstance(o, np.ndarray):
            return o.tolist()
        return str(o)
    with open(ruta, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False, default=_c)


def cargar_json(ruta):
    with open(ruta, encoding="utf-8") as f:
        return json.load(f)


# --- Filas del ensayo ------------------------------------------------------------
def _texto(v):
    """Valor de un campo como texto estable (1.0 -> '1')."""
    if isinstance(v, (float, np.floating)) and np.isfinite(v) and float(v).is_integer():
        return str(int(v))
    return str(v).strip()


def atributos_croquis(cfg):
    """Nombres de los atributos del archivo vectorial que se copian a los resultados (croquis.atributos)."""
    return [a for a in (cfg.get("croquis", {}).get("atributos") or []) if a != "id_fila"]


def valores_atributos(cfg, fila):
    """{atributo: valor} de una hilera (fila de cargar_filas o dict de marcos.json)."""
    if isinstance(fila, dict) and "atributos" in fila:
        return dict(fila["atributos"])
    return {a: (fila[a].item() if hasattr(fila[a], "item") else fila[a]) for a in atributos_croquis(cfg)}


def cargar_filas(cfg, solo_prueba=None, todas=False):
    """Hileras del archivo vectorial del ensayo (capa de LÍNEAS, una por hilera), según la sección 'croquis' de la ficha:
    - columna 'id_fila' = croquis.campo_id como texto (o, si id_compuesto, campo_id + atributos);
    - atributos de texto normalizados en memoria (espacios, minúsculas) y 'reemplazos' de erratas;
    - 'filtro' {campo: valor | [valores]} para quedarse con las hileras del ensayo;
    - el archivo vectorial se reproyecta en memoria al EPSG de trabajo si viene en otro CRS.
    Selección: todas=True -> todas las del filtro; si el cfg lleva '_filas_bloque' (modo completo) -> esas;
    si no, modo 'prueba' -> filas_prueba_id; modo 'completo' -> todas las del filtro.
    El fichero original no se modifica."""
    import geopandas as gpd
    cq = cfg["croquis"]
    g = gpd.read_file(cfg["rutas"]["filas_gpkg"], layer=cfg["rutas"]["filas_capa"]) if cfg["rutas"]["filas_capa"] \
        else gpd.read_file(cfg["rutas"]["filas_gpkg"])
    campos = [cq["campo_id"]] + atributos_croquis(cfg) + list((cq.get("filtro") or {}).keys())
    faltan = [c for c in campos if c not in g.columns]
    if faltan:
        raise KeyError("El archivo vectorial no tiene los campos %s. Campos disponibles: %s"
                       % (faltan, [c for c in g.columns if c != "geometry"]))
    tipos = set(g.geom_type.unique())
    if not tipos <= {"LineString", "MultiLineString"}:
        raise ValueError("El archivo vectorial debe ser de líneas (una por hilera); tiene %s" % sorted(tipos))
    if g.geom_type.eq("MultiLineString").any():
        from shapely.ops import linemerge
        g["geometry"] = [linemerge(x) if x.geom_type == "MultiLineString" else x for x in g.geometry]
        if not g.geom_type.eq("LineString").all():
            raise ValueError("Hay hileras multiparte que no se pueden unir en una sola línea")
    if g.crs is not None and g.crs.to_epsg() != cfg["crs"]["epsg"]:
        g = g.to_crs(cfg["crs"]["epsg"])
    for a in set(atributos_croquis(cfg)) | set((cq.get("filtro") or {}).keys()):
        es_texto = g[a].dtype == object or str(g[a].dtype).startswith("str")
        if es_texto:
            g[a] = g[a].astype(str).str.strip()
            if cq.get("normalizar_texto", True):
                g[a] = g[a].str.lower()
        if (cq.get("reemplazos") or {}).get(a):
            g[a] = g[a].replace(cq["reemplazos"][a])
    for a, v in (cq.get("filtro") or {}).items():
        vals = v if isinstance(v, (list, tuple)) else [v]
        norm = (lambda s: s.lower()) if cq.get("normalizar_texto", True) else (lambda s: s)
        vals = {norm(_texto(x)) for x in vals}
        g = g[g[a].map(lambda x: norm(_texto(x))).isin(vals)]
    g = g.copy()
    if g.empty:
        raise ValueError("Ninguna hilera del archivo vectorial cumple el filtro %s" % cq.get("filtro"))
    ids = g[cq["campo_id"]].map(_texto)
    if cq.get("id_compuesto"):
        ids = ids.str.cat([g[a].map(_texto) for a in atributos_croquis(cfg)], sep="_")
    g["id_fila"] = ids.str.replace(r"[^0-9A-Za-z_.\-]", "_", regex=True)   # se usa en nombres de carpeta
    if not g["id_fila"].is_unique:
        rep = sorted(g["id_fila"][g["id_fila"].duplicated()].unique())
        raise ValueError("El id de hilera '%s' se repite (%s): añade un filtro en croquis.filtro o pon "
                         "croquis.id_compuesto: true" % (cq["campo_id"], rep))
    if not todas:
        if cfg.get("_filas_bloque"):
            g = g[g["id_fila"].isin([str(i) for i in cfg["_filas_bloque"]])]
        elif cfg["recorte"]["modo"] == "prueba":
            pid = [str(i) for i in cfg["recorte"]["filas_prueba_id"]]
            g = g[g["id_fila"].isin(pid)]
            if g.empty:
                raise ValueError("Modo prueba: ninguna de las hileras_prueba %s está en el archivo vectorial filtrado" % pid)
    return g.reset_index(drop=True)


class MarcoHilera:
    """Marco (s, t) de una hilera: s a lo largo del eje (m desde 'origen'), t perpendicular
    (positivo a la izquierda de la dirección de avance)."""

    def __init__(self, origen_xy, direccion_xy, largo=None):
        self.o = np.asarray(origen_xy, float)
        d = np.asarray(direccion_xy, float)
        self.e = d / np.linalg.norm(d)
        self.n = np.array([-self.e[1], self.e[0]])
        self.largo = largo

    @classmethod
    def de_linea(cls, geom):
        c = np.array(geom.coords)[:, :2]
        d = c[-1] - c[0]
        return cls(c[0], d, float(np.linalg.norm(d)))

    def st(self, xy):
        q = np.asarray(xy, float)[..., :2] - self.o
        return q @ self.e, q @ self.n

    def xy(self, s, t):
        return self.o + np.asarray(s)[..., None] * self.e + np.asarray(t)[..., None] * self.n

    def a_dict(self):
        return dict(origen=self.o.tolist(), e=self.e.tolist(), largo=self.largo)

    @classmethod
    def de_dict(cls, d):
        return cls(d["origen"], d["e"], d.get("largo"))


def reajustar_eje_pca(marco, xy):
    """Reajusta el eje con PCA sobre los puntos de copa (xy). Devuelve (nuevo marco,
    giro respecto al GPKG en grados, desplazamiento lateral del eje en el origen en m)."""
    c = xy.mean(0)
    _, _, vt = np.linalg.svd(xy - c, full_matrices=False)
    d = vt[0]
    if d @ marco.e < 0:
        d = -d
    o2 = c + d * ((marco.o - c) @ d)
    m2 = MarcoHilera(o2, d, marco.largo)
    giro = float(np.degrees(np.arccos(np.clip(d @ marco.e, -1, 1))))
    return m2, giro, float(np.linalg.norm(o2 - marco.o))


# --- LAS -------------------------------------------------------------------------
def escribir_las_desde(las_ref, idx, ruta, dims_extra=None):
    """Escribe el subconjunto 'idx' (máscara o índices) de un laspy.LasData con dimensiones
    extra opcionales (dict nombre -> array del tamaño del subconjunto)."""
    import laspy
    h = laspy.LasHeader(point_format=las_ref.header.point_format, version=las_ref.header.version)
    h.offsets, h.scales = las_ref.header.offsets, las_ref.header.scales
    try:
        crs = las_ref.header.parse_crs()
        if crs is not None:
            h.add_crs(crs)
    except Exception:
        pass
    for nom, arr in (dims_extra or {}).items():
        h.add_extra_dim(laspy.ExtraBytesParams(name=nom, type=np.asarray(arr).dtype))
    sub = las_ref.points[idx]
    n = len(sub)
    out = laspy.LasData(h)
    out.points = laspy.ScaleAwarePointRecord.zeros(n, header=h)
    for dim in las_ref.point_format.dimension_names:
        if dim in out.point_format.dimension_names and dim not in (dims_extra or {}):
            out[dim] = sub[dim]
    for nom, arr in (dims_extra or {}).items():
        out[nom] = arr
    out.write(ruta)
    return n


def run_pdal(etapas, log=None):
    """Ejecuta un pipeline PDAL (lista de etapas). Devuelve (nº de puntos, pipeline)."""
    import pdal
    p = pdal.Pipeline(json.dumps(etapas))
    n = p.execute()
    if log:
        log.info("PDAL: %d puntos", n)
    return n, p


# --- Figuras ---------------------------------------------------------------------
def png_matriz(a, ruta, titulo="", extent=None, cmap="viridis", cbar="", vmin=None, vmax=None):
    """PNG de una matriz 2D (backend Agg; si matplotlib falla, Pillow)."""
    a = np.asarray(a, float)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 6), dpi=110)
        im = ax.imshow(a, extent=extent, cmap=cmap, vmin=vmin, vmax=vmax, origin="upper")
        ax.set_title(titulo)
        fig.colorbar(im, ax=ax, label=cbar)
        fig.tight_layout()
        fig.savefig(ruta)
        plt.close(fig)
    except Exception:
        from PIL import Image
        m = np.isfinite(a)
        lo = float(np.nanmin(a[m])) if vmin is None and m.any() else (vmin or 0)
        hi = float(np.nanmax(a[m])) if vmax is None and m.any() else (vmax or 1)
        g = np.clip((np.nan_to_num(a, nan=lo) - lo) / max(hi - lo, 1e-9), 0, 1)
        Image.fromarray((g * 255).astype(np.uint8)).resize((a.shape[1] * 4, a.shape[0] * 4), Image.NEAREST).save(ruta)


def png_lineas(series, ruta, titulo="", xlabel="", ylabel=""):
    """series: lista de (x, y, etiqueta)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 4.5), dpi=110)
    for x, y, lab in series:
        ax.plot(x, y, label=lab, lw=1.2)
    ax.set_title(titulo)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(alpha=.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(ruta)
    plt.close(fig)


def copiar_config(cfg, destino):
    """Guarda la ficha del ensayo tal cual y la configuración efectiva (parámetros + ficha + ajustes)."""
    if cfg.get("_ruta_ensayo"):
        shutil.copy2(cfg["_ruta_ensayo"], os.path.join(destino, "ficha_ensayo_usada.yaml"))
    efectiva = {k: v for k, v in cfg.items() if not k.startswith("_")}
    with open(os.path.join(destino, "config_efectiva.yaml"), "w", encoding="utf-8") as f:
        yaml.safe_dump(efectiva, f, allow_unicode=True, sort_keys=False)
