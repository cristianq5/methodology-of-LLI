# -*- coding: utf-8 -*-
"""Geometría del recorte: polígono por tramo de hilera y pasos del sensor por cada calle.

QA-19/QA-20: un PASO es un recorrido continuo del sensor por UNA calle en UN sentido. Los pasos se cortan por
hueco temporal, por cambio de lado respecto al eje (signo de t) y por inversión del rumbo a lo largo de la hilera
(media vuelta en la cabecera). Nunca se fusionan pasos de lado o sentido distinto. Los pasos de hileras vecinas que
son el mismo recorrido (misma calle, mismo sentido, intervalos solapados) se unifican en una VENTANA GLOBAL, de modo
que cada calle se decodifica y calibra una sola vez y se comparte entre las dos hileras que la bordean.
"""
import numpy as np
from shapely.geometry import Polygon
from shapely.ops import unary_union

from utils import MarcoHilera, tramo_s, cargar_filas

_CACHE_FILAS = {}


def tramo_rectangulo(marco, s0, s1, ancho):
    """Rectángulo de 's0..s1' a lo largo del eje y ±ancho de t."""
    pts = [marco.xy(s0, -ancho), marco.xy(s1, -ancho), marco.xy(s1, ancho), marco.xy(s0, ancho)]
    return Polygon(pts)


def poligono_recorte(cfg, filas):
    ancho = cfg["recorte"]["buffer_recorte_m"]
    polis = []
    for _, f in filas.iterrows():
        m = MarcoHilera.de_linea(f.geometry)
        s0, s1 = tramo_s(cfg, m)
        polis.append(tramo_rectangulo(m, s0, s1, ancho))
    return unary_union(polis)


def segmentar_pasos(t_tray, s, t, s0, s1, pad_s, t_abs, hueco_s, min_muestras=20, min_tramo_s=3.0, suav_s=1.0):
    """Pasos del sensor (s, t = coordenadas respecto a la hilera) por el tramo [s0, s1].

    Se consideran las muestras con s en [s0-pad_s, s1+pad_s] y t_abs[0] <= |t| <= t_abs[1] (calle adyacente).
    Se corta por: hueco temporal > hueco_s, cambio de signo de t (lado) y cambio de signo de la velocidad a lo
    largo de la hilera (sentido). Devuelve lista de dict(t0, t1, lado, sentido, t_med, s_ini, s_fin)."""
    dt_med = float(np.median(np.diff(t_tray)))
    k = max(1, int(round(suav_s / dt_med / 2)))
    v = np.zeros_like(s)
    v[k:-k] = s[2 * k:] - s[:-2 * k]                      # variación de s en ~suav_s (sentido de avance)
    m = (s > s0 - pad_s) & (s < s1 + pad_s) & (np.abs(t) >= t_abs[0]) & (np.abs(t) <= t_abs[1])
    idx = np.nonzero(m)[0]
    if len(idx) == 0:
        return []
    lado = np.sign(t[idx]).astype(int)
    sent = np.where(v[idx] >= 0, 1, -1)
    corte = (np.diff(t_tray[idx]) > hueco_s) | (np.diff(lado) != 0) | (np.diff(sent) != 0)
    out = []
    for sg in np.split(idx, np.nonzero(corte)[0] + 1):
        if len(sg) < min_muestras or abs(s[sg[-1]] - s[sg[0]]) < min_tramo_s:
            continue
        out.append(dict(t0=float(t_tray[sg[0]]), t1=float(t_tray[sg[-1]]), lado=int(np.sign(np.median(t[sg]))),
                        sentido=int(np.sign(s[sg[-1]] - s[sg[0]])), t_med=float(np.median(t[sg])),
                        s_ini=float(s[sg[0]]), s_fin=float(s[sg[-1]])))
    return out


def pasos_sensor(cfg, tray, marco, s0, s1):
    """Pasos del sensor por el tramo [s0,s1] de una hilera (ver segmentar_pasos)."""
    r = cfg["recorte"]
    s, t = marco.st(tray.p[:, :2])
    return segmentar_pasos(tray.t, s, t, s0, s1, r["ventana_tiempo_s_pad_m"],
                           (r["ventana_t_abs_min_m"], r["ventana_tiempo_t_max_m"]), r["ventana_hueco_s"])


def _todas(cfg):
    clave = cfg["rutas"]["filas_gpkg"] + "|" + str(cfg["recorte"]["modo"])
    if clave not in _CACHE_FILAS:
        c = {k: v for k, v in cfg.items() if k not in ("_bloque", "_filas_bloque", "_bloques_sel")}
        _CACHE_FILAS[clave] = cargar_filas(c, todas=True)
    return _CACHE_FILAS[clave]


def ventanas_globales(cfg, tray):
    """Pasos de TODAS las hileras unificados por recorrido (misma calle y sentido, intervalos solapados > 50 %).
    Cada ventana lleva la lista de (fila, lado) a las que sirve."""
    pasos = []
    for _, f in _todas(cfg).iterrows():
        m = MarcoHilera.de_linea(f.geometry)
        s0, s1 = tramo_s(cfg, m)
        for p in pasos_sensor(cfg, tray, m, s0, s1):
            p["fila"] = str(f["id_fila"])
            # rumbo del paso en el MUNDO: el 'sentido' es relativo a la línea de cada hilera, y dos hileras vecinas pueden
            # estar digitalizadas en sentidos opuestos; para unificar pasos se compara el rumbo real del sensor
            p["rumbo"] = p["sentido"] * m.e
            pasos.append(p)
    pasos.sort(key=lambda p: p["t0"])
    glob = []
    for p in pasos:
        for g in glob:
            ov = min(g["t1"], p["t1"]) - max(g["t0"], p["t0"])
            if float(g["rumbo"] @ p["rumbo"]) > 0 and ov > 0.5 * min(g["t1"] - g["t0"], p["t1"] - p["t0"]) and \
                    np.sign(g["t_med"]) * np.sign(p["t_med"]) != 0 and abs(g["t0"] - p["t0"]) < 30:
                # misma calle: el mismo recorrido visto desde dos hileras (lados opuestos respecto a cada eje)
                g["t0"], g["t1"] = min(g["t0"], p["t0"]), max(g["t1"], p["t1"])
                g["sirve"].append((p["fila"], p["lado"]))
                g["sentido_fila"][p["fila"]] = p["sentido"]
                break
        else:
            glob.append(dict(t0=p["t0"], t1=p["t1"], sentido=p["sentido"], rumbo=p["rumbo"], t_med=p["t_med"],
                             sirve=[(p["fila"], p["lado"])], sentido_fila={p["fila"]: p["sentido"]}))
    return glob


def ventanas_tiempo(cfg, tray, filas):
    """Ventanas de tiempo (sin solapes) de los pasos que sirven a las hileras 'filas' y lista de pasos por fila.
    Los límites se calculan sobre la lista GLOBAL de pasos (iguales para todos los bloques, así la calibración de
    una calle se comparte). Cada ventana es un único paso (un lado, un sentido). Devuelve (ventanas, pasos)."""
    pad_t = cfg["recorte"]["ventana_tiempo_pad_s"]
    ids = {str(i) for i in filas["id_fila"]}
    glob = sorted(ventanas_globales(cfg, tray), key=lambda g: g["t0"])
    ven = []
    for i, g in enumerate(glob):
        a, b = g["t0"] - pad_t, g["t1"] + pad_t
        if i > 0 and a < ven[-1][1]:                      # nunca solapar ni fusionar ventanas contiguas
            mid = 0.5 * (glob[i - 1]["t1"] + g["t0"])
            ven[-1][1] = mid
            a = mid
        ven.append([a, b])
    sel, pasos = [], []
    for g, w in zip(glob, ven):
        sirve = [(f, lado) for f, lado in g["sirve"] if f in ids]
        if sirve:
            sel.append(w)
            pasos.extend(dict(t0=g["t0"], t1=g["t1"], lado=lado, sentido=g["sentido_fila"][f], t_med=g["t_med"], fila=f) for f, lado in sirve)
    return sel, pasos
