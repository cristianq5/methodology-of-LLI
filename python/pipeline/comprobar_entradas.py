# -*- coding: utf-8 -*-
"""Comprueba, ANTES de lanzar el pipeline, que están todos los datos de un ensayo y que se pueden leer.

    python comprobar_entradas.py --ensayo ficha.yaml

Solo lee (nunca modifica las entradas) y tarda segundos: lee cabeceras y unas pocas líneas o paquetes, no las nubes
enteras. Revisa la ficha, la nube LAS, cada fichero de la carpeta de brutos del RS30 (PosI, high_freq, rtk.gga, pcap,
LAS-1.EP), que los tiempos de todos se solapen, el archivo vectorial de hileras y la carpeta de salida.
Cada comprobación sale como OK, AVISO (se puede lanzar, pero conviene mirarlo) o ERROR (el pipeline fallaría).
Escribe el informe en <salida>/logs/comprobar_entradas.json (si la carpeta de salida se puede crear) y termina con
código 1 si hay algún ERROR. Qué debe contener la carpeta de entradas: docs/ENTRADAS.md.
"""
import utils  # noqa: F401 (primero: entorno)
import argparse
import glob
import json
import os
import struct
import sys

import numpy as np
import yaml

MAGIA_PCAP = (0xA1B2C3D4, 0xA1B23C4D)        # pcap clásico (µs / ns), little endian


class Informe:
    def __init__(self):
        self.filas = []

    def add(self, nivel, que, detalle=""):
        self.filas.append(dict(nivel=nivel, comprobacion=que, detalle=str(detalle)))
        print("%-5s  %-42s %s" % (nivel, que, detalle), flush=True)

    def ok(self, q, d=""):
        self.add("OK", q, d)

    def aviso(self, q, d=""):
        self.add("AVISO", q, d)

    def error(self, q, d=""):
        self.add("ERROR", q, d)

    @property
    def n_errores(self):
        return sum(f["nivel"] == "ERROR" for f in self.filas)


def _primera_ultima_linea(ruta, saltar=0):
    """Primera línea de datos (tras 'saltar' líneas) y última línea no vacía, sin leer el fichero entero."""
    with open(ruta, "rb") as f:
        for _ in range(saltar):
            f.readline()
        prim = f.readline().decode("latin-1").strip()
        f.seek(0, os.SEEK_END)
        tam = f.tell()
        f.seek(max(0, tam - 4096))
        ult = [l for l in f.read().decode("latin-1").splitlines() if l.strip()][-1].strip()
    return prim, ult


def _numeros(linea):
    """Valores numéricos del principio de la línea (se para en el primer campo de texto, p. ej. 'Fixed' del PosI)."""
    out = []
    for x in linea.replace(",", " ").split():
        try:
            out.append(float(x))
        except ValueError:
            break
    return out


def comprobar_texto_tiempo(inf, nombre, ruta, saltar, min_cols):
    """Fichero de trayectoria en texto: nº de columnas numéricas y rango de tiempo (primera columna)."""
    try:
        a, b = _primera_ultima_linea(ruta, saltar)
        na, nb = _numeros(a), _numeros(b)
    except Exception as e:
        inf.error(nombre + ": lectura", "%s (%r)" % (ruta, e))
        return None
    if len(na) < min_cols:
        inf.error(nombre + ": columnas", "%d columnas numéricas; se esperan al menos %d" % (len(na), min_cols))
        return None
    t0, t1 = na[0], nb[0]
    if not t1 > t0:
        inf.error(nombre + ": tiempo", "el tiempo no avanza (%.3f -> %.3f)" % (t0, t1))
        return None
    inf.ok(nombre, "%d columnas; t %.1f -> %.1f (%.0f s)" % (len(na), t0, t1, t1 - t0))
    return t0, t1


def _registros_xt32(f, off, n_bytes):
    """Inicios de los registros pcap de una trama del XT32 (862 B, payload que empieza por EE FF) en [off, off + n_bytes)."""
    import pcap_hesai as ph
    f.seek(off)
    buf = f.read(n_bytes)
    pat = struct.pack("<II", 862, 862)            # incl_len y orig_len (bytes 8-15 de la cabecera del registro)
    out, i = [], buf.find(pat)
    while i >= 0:
        p = i + 8 + ph.OFF_ETH                      # payload = registro + 16 B de cabecera + 42 B Ethernet/IP/UDP
        if i >= 8 and p + 2 <= len(buf) and buf[p:p + 2] == b"\xee\xff" and i - 8 + ph.TAM_REG <= len(buf):
            out.append(off + i - 8)
        i = buf.find(pat, i + 1)
    return out


def comprobar_pcap(inf, ruta):
    """Cabecera global del pcap, registros de trama del XT32 (862 B con payload EE FF) al principio y al final y rango de
    tiempo (cabecera de cada registro)."""
    import pcap_hesai as ph
    try:
        tam = os.path.getsize(ruta)
        with open(ruta, "rb") as f:
            magia = struct.unpack("<I", f.read(4))[0]
            if magia not in MAGIA_PCAP:
                inf.error("pcap: formato", "no es un pcap clásico little endian (magia 0x%08X)" % magia)
                return None
            ini = _registros_xt32(f, 24, 2 ** 20)
            fin = _registros_xt32(f, max(24, tam - 2 ** 20), 2 ** 20)
            if not ini or not fin:
                inf.error("pcap: paquetes del XT32", "no hay registros de 862 B con payload EE FF en el primer o el último MB")
                return None
            t0, t1 = ph._tiempo_host(f, ini[0]), ph._tiempo_host(f, fin[-1])
    except Exception as e:
        inf.error("pcap: lectura", repr(e))
        return None
    inf.ok("pcap", "%.2f GB; paquetes del Hesai XT32; t (host) %.1f -> %.1f (%.0f s)" % (tam / 2 ** 30, t0, t1, t1 - t0))
    return t0, t1


def comprobar_las(inf, ruta, epsg_ficha):
    import laspy
    try:
        with laspy.open(ruta) as fl:
            h = fl.header
            dims = list(h.point_format.dimension_names)
            crs = h.parse_crs()
            n = int(h.point_count)
            if "gps_time" not in dims:
                inf.error("LAS: gps_time", "el formato de punto %d no tiene gps_time" % h.point_format.id)
                return None, None
            gt = None
            for tr in fl.chunk_iterator(1_000_000):
                gt = np.asarray(tr.gps_time)
                break
    except Exception as e:
        inf.error("LAS: lectura", "%s (%r)" % (ruta, e))
        return None, None
    epsg = crs.to_epsg() if crs is not None else None
    inf.ok("LAS", "%d puntos, formato %d, versión %s, gps_time presente" % (n, h.point_format.id, h.version))
    if epsg_ficha:
        if epsg and int(epsg) != int(epsg_ficha):
            inf.aviso("LAS: CRS", "cabecera EPSG:%s, ficha EPSG:%s (manda la ficha; comprueba que es el correcto)" % (epsg, epsg_ficha))
        else:
            inf.ok("LAS: CRS", "EPSG:%s (ficha)" % epsg_ficha)
    elif epsg:
        inf.ok("LAS: CRS", "EPSG:%s (cabecera)" % epsg)
    else:
        inf.error("LAS: CRS", "la cabecera no trae CRS: indica crs.epsg en la ficha")
    if gt is None or len(gt) == 0:
        inf.error("LAS: puntos", "la nube está vacía")
        return None, epsg
    rango = (float(gt.min()), float(gt.max()))
    if np.all(gt == 0):
        inf.error("LAS: gps_time", "gps_time vale 0 en todos los puntos leídos")
        return None, epsg
    return rango, (epsg_ficha or epsg)


def comprobar_ep(inf, ruta):
    import calibracion_global as cg
    try:
        el = cg.leer_elevaciones_ep(ruta)
    except Exception as e:
        inf.error("LAS-1.EP", "no se puede leer la tabla de elevaciones (%r)" % (e,))
        return
    if not (np.all(np.isfinite(el)) and -30 < el.min() < el.max() < 30):
        inf.error("LAS-1.EP", "elevaciones fuera de rango: %.2f .. %.2f°" % (el.min(), el.max()))
        return
    inf.ok("LAS-1.EP", "32 elevaciones de %.2f a %.2f°" % (el.min(), el.max()))


def comprobar_rtk(inf, ruta):
    n, n_fijo = 0, 0
    with open(ruta, encoding="latin-1") as f:
        for i, l in enumerate(f):
            p = l.split(",")
            if len(p) > 6 and p[0].endswith("GGA") and p[6].isdigit():
                n += 1
                n_fijo += p[6] == "4"
    if n == 0:
        inf.aviso("rtk.gga", "no hay sentencias GGA (el 01_qc no podrá dar el % de RTK fijo)")
    else:
        inf.ok("rtk.gga", "%d épocas GGA, %.1f %% con RTK fijo" % (n, 100 * n_fijo / n))


def solape(inf, nombre, a, b, minimo=0.5):
    """a, b = (t0, t1). Fracción de 'a' que cae dentro de 'b'."""
    if a is None or b is None:
        return
    s = max(0.0, min(a[1], b[1]) - max(a[0], b[0]))
    frac = s / max(a[1] - a[0], 1e-9)
    texto = "%.0f %% dentro (%.1f-%.1f frente a %.1f-%.1f)" % (100 * frac, a[0], a[1], b[0], b[1])
    if frac >= minimo:
        inf.ok(nombre, texto)
    elif frac > 0:
        inf.aviso(nombre, texto)
    else:
        inf.error(nombre, "no se solapan en el tiempo: " + texto + ". ¿Son de la misma captura?")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ensayo", default=os.environ.get(utils.VAR_ENTORNO_ENSAYO), help="ficha del ensayo")
    ap.add_argument("--parametros", default=None, help="parametros.yaml alternativo")
    a = ap.parse_args(argv)
    inf = Informe()
    print("Comprobación de entradas del pipeline LLI\n" + "-" * 78)
    if not a.ensayo or not os.path.exists(a.ensayo):
        inf.error("ficha del ensayo", "no existe: %s (usa --ensayo)" % a.ensayo)
        return 1
    ruta_f = os.path.abspath(a.ensayo)
    base = os.path.dirname(ruta_f)
    try:
        with open(ruta_f, encoding="utf-8") as f:
            e = yaml.safe_load(f) or {}
    except Exception as ex_:
        inf.error("ficha del ensayo", "YAML no válido (%r)" % ex_)
        return 1
    faltan = [s for s in ("ensayo", "entradas", "croquis") if s not in e]
    if faltan:
        inf.error("ficha: secciones", "faltan %s (ver config/ensayo_plantilla.yaml)" % faltan)
        return 1
    inf.ok("ficha", ruta_f)
    try:
        par = utils.cargar_parametros(utils._ruta(base, e.get("parametros")) or a.parametros)
        par = utils._fusionar(par, e.get("ajustes") or {})
        inf.ok("parametros.yaml + ajustes de la ficha", "modo de calibración: %s" % par["pcap"].get("modo_calibracion"))
    except Exception as ex_:
        inf.error("parametros.yaml", repr(ex_))
        return 1
    ent = e["entradas"]

    # ---- nube LAS
    ruta_las = utils._ruta(base, ent.get("las"))
    rango_las, epsg = (None, None)
    if not ruta_las or not os.path.exists(ruta_las):
        inf.error("LAS", "no existe: %s" % ruta_las)
    else:
        rango_las, epsg = comprobar_las(inf, ruta_las, (e.get("crs") or {}).get("epsg"))

    # ---- ficheros de la carpeta de brutos
    brutos = utils._ruta(base, ent.get("brutos"))
    elegidos = ent.get("brutos_archivos") or {}
    if brutos and not os.path.isdir(brutos):
        inf.error("carpeta de brutos", "no existe: %s" % brutos)
        brutos = None
    rutas = {}
    cg = par["pcap"].get("calibracion_global") or {}
    necesita_ep = par["pcap"].get("modo_calibracion") == "global_bloque" and cg.get("elevaciones") == "EP"
    for nombre, (patron, oblig) in utils.ARCHIVOS_BRUTOS.items():
        oblig = oblig or (nombre == "calib_ep" and necesita_ep)
        r = utils._ruta(base, elegidos.get(nombre))
        if r:
            if os.path.exists(r):
                rutas[nombre] = r
                inf.ok("brutos: " + nombre, r + " (indicado en la ficha)")
            else:
                inf.error("brutos: " + nombre, "indicado en la ficha pero no existe: %s" % r)
            continue
        if not brutos:
            (inf.error if oblig else inf.aviso)("brutos: " + nombre, "sin carpeta de brutos ni ruta en brutos_archivos")
            continue
        c = sorted(glob.glob(os.path.join(brutos, patron)))
        if len(c) == 1:
            rutas[nombre] = c[0]
            inf.ok("brutos: " + nombre, os.path.relpath(c[0], brutos))
        elif not c:
            (inf.error if oblig else inf.aviso)("brutos: " + nombre, "no se encuentra %s%s" % (
                patron, "" if oblig else " (opcional)"))
        else:
            inf.error("brutos: " + nombre, "%d candidatos para %s: %s. Indica uno en entradas.brutos_archivos.%s"
                      % (len(c), patron, [os.path.basename(x) for x in c], nombre))

    t_posi = comprobar_texto_tiempo(inf, "PosI", rutas["trayectoria_posi"], 2, 7) if "trayectoria_posi" in rutas else None
    t_hf = comprobar_texto_tiempo(inf, "high_freq", rutas["trayectoria_hf"], 0, 8) if "trayectoria_hf" in rutas else None
    t_pcap = comprobar_pcap(inf, rutas["pcap"]) if "pcap" in rutas else None
    if "rtk_gga" in rutas:
        comprobar_rtk(inf, rutas["rtk_gga"])
    if "calib_ep" in rutas:
        comprobar_ep(inf, rutas["calib_ep"])
    elif not necesita_ep:
        inf.ok("LAS-1.EP", "no hace falta con esta configuración")

    # ---- tiempos: todo debe ser de la misma captura
    solape(inf, "tiempo: LAS (primer trozo) en high_freq", rango_las, t_hf, minimo=0.99)
    solape(inf, "tiempo: high_freq en PosI", t_hf, t_posi)
    solape(inf, "tiempo: high_freq en pcap", t_hf, t_pcap)

    # ---- archivo vectorial de hileras
    cq = utils._fusionar(dict(capa=None, campo_id="id_fila", atributos=[], filtro={}, normalizar_texto=True,
                              reemplazos={}, id_compuesto=False), e["croquis"])
    ruta_v = utils._ruta(base, cq.get("ruta"))
    if not ruta_v or not os.path.exists(ruta_v):
        inf.error("archivo vectorial", "no existe: %s" % ruta_v)
    elif not epsg:
        inf.error("archivo vectorial", "sin CRS de trabajo (LAS sin CRS y ficha sin crs.epsg): no se puede comprobar")
    else:
        cq["ruta"] = ruta_v
        ej = e.get("ejecucion_ensayo") or {}
        cfg = dict(par, croquis=cq, rutas=dict(filas_gpkg=ruta_v, filas_capa=cq.get("capa")), crs=dict(epsg=int(epsg)),
                   recorte=dict(par["recorte"], modo=ej.get("modo", "prueba"),
                                filas_prueba_id=[str(x) for x in (ej.get("hileras_prueba") or [])]))
        try:
            g = utils.cargar_filas(cfg, todas=True)
            inf.ok("archivo vectorial", "%d hileras tras el filtro %s; ids %s" % (
                len(g), cq.get("filtro") or "{}", ", ".join(map(str, list(g["id_fila"])[:15])) + (" ..." if len(g) > 15 else "")))
            largo = g.geometry.length
            if (largo < 2).any():
                inf.aviso("archivo vectorial: longitud", "hileras de menos de 2 m: %s" % list(g.loc[largo < 2, "id_fila"]))
            if cfg["recorte"]["modo"] == "prueba":
                pid = cfg["recorte"]["filas_prueba_id"]
                no = [i for i in pid if i not in set(g["id_fila"])]
                if not pid:
                    inf.error("hileras de prueba", "modo prueba sin ejecucion_ensayo.hileras_prueba")
                elif no:
                    inf.error("hileras de prueba", "no están en el archivo vectorial filtrado: %s" % no)
                else:
                    inf.ok("hileras de prueba", pid)
            if ruta_las and os.path.exists(ruta_las):
                import laspy
                with laspy.open(ruta_las) as fl:
                    x0, y0 = fl.header.mins[:2]
                    x1, y1 = fl.header.maxs[:2]
                gx0, gy0, gx1, gy1 = g.total_bounds
                if gx0 >= x0 - 1 and gx1 <= x1 + 1 and gy0 >= y0 - 1 and gy1 <= y1 + 1:
                    inf.ok("hileras dentro de la nube", "la extensión de las hileras cae dentro de la del LAS")
                elif gx1 < x0 or gx0 > x1 or gy1 < y0 or gy0 > y1:
                    inf.error("hileras dentro de la nube", "las hileras no tocan la nube: ¿CRS distinto o captura equivocada?")
                else:
                    inf.aviso("hileras dentro de la nube", "parte de las hileras queda fuera de la extensión del LAS")
        except Exception as ex_:
            inf.error("archivo vectorial", str(ex_))

    # ---- salida
    sal = e.get("salidas") or {}
    ruta_s = utils._ruta(base, sal.get("carpeta") or "salidas")
    ruta_logs = utils._ruta(base, sal.get("logs")) or os.path.join(ruta_s, "logs")
    if brutos and os.path.commonpath([os.path.abspath(ruta_s), os.path.abspath(brutos)]) == os.path.abspath(brutos):
        inf.error("carpeta de salida", "está dentro de la carpeta de brutos: %s" % ruta_s)
    else:
        try:
            os.makedirs(ruta_logs, exist_ok=True)
            import shutil
            libre = shutil.disk_usage(ruta_s).free / 2 ** 30
            (inf.ok if libre > 20 else inf.aviso)("carpeta de salida", "%s (%.0f GB libres%s)" % (
                ruta_s, libre, "" if libre > 20 else "; un bloque completo ocupa decenas de GB"))
        except Exception as ex_:
            inf.error("carpeta de salida", "no se puede crear %s (%r)" % (ruta_s, ex_))
            ruta_logs = None
    print("-" * 78)
    n_av = sum(f["nivel"] == "AVISO" for f in inf.filas)
    print("RESULTADO: %d errores, %d avisos. %s" % (inf.n_errores, n_av,
          "Se puede lanzar el pipeline." if inf.n_errores == 0 else "Corrige los ERROR antes de lanzar el pipeline."))
    if ruta_logs and os.path.isdir(ruta_logs):
        with open(os.path.join(ruta_logs, "comprobar_entradas.json"), "w", encoding="utf-8") as f:
            json.dump(dict(ficha=ruta_f, errores=inf.n_errores, avisos=n_av, comprobaciones=inf.filas), f, indent=2, ensure_ascii=False)
    return 1 if inf.n_errores else 0


if __name__ == "__main__":
    sys.exit(main())
