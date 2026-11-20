# -*- coding: utf-8 -*-
"""01_qc: control de calidad previo (no recorta ni modifica nada).

- Cabecera del LAS (n puntos, extensión, CRS) y muestra de gps_time frente a la trayectoria.
- Trayectoria post-SLAM (hf) llevada al CRS del ensayo: ajuste frente al PosI, huecos, % RTK fijo (rtk.gga, si existe).
- Cobertura BILATERAL por hilera: la trayectoria debe recorrer las dos calles de cada fila.
- Mapa de control (trayectoria + filas + recorte).
La comprobación de altura del sensor sobre el DTM se hace en 05_dtm (cuando existe el DTM).
"""
import utils  # noqa: F401 (primero)
import os
import numpy as np
import pandas as pd
import laspy

from utils import cargar_config, dir_salida, Paso, MarcoHilera, cargar_filas, guardar_json, valores_atributos
from trayectoria import cargar_trayectoria
from geometria import poligono_recorte, ventanas_tiempo

NOMBRE = "01_qc"


def pct_rtk_fijo(ruta):
    q = []
    with open(ruta, encoding="latin-1") as f:
        for l in f:
            p = l.split(",")
            if len(p) > 6 and p[0].endswith("GGA") and p[6].isdigit():
                q.append(int(p[6]))
    q = np.array(q)
    return (float(100 * np.mean(q == 4)) if len(q) else float("nan")), len(q)


def main(cfg=None):
    cfg = cfg or cargar_config()
    out = dir_salida(cfg, NOMBRE)
    with Paso(cfg, NOMBRE, entrada=cfg["rutas"]["las_entrada"], parametros=cfg["qc"]) as p:
        res = {}
        # --- cabecera LAS y muestra de tiempos (solo lectura, un trozo)
        with laspy.open(cfg["rutas"]["las_entrada"]) as fl:
            h = fl.header
            res["las"] = dict(n_puntos=int(h.point_count), mins=h.mins.tolist(), maxs=h.maxs.tolist(),
                              version=str(h.version), formato=int(h.point_format.id),
                              crs=str(h.parse_crs().to_epsg()) if h.parse_crs() else None)
            tr0 = None
            for trozo in fl.chunk_iterator(2_000_000):
                tr0 = (float(trozo.gps_time.min()), float(trozo.gps_time.max()))
                break
            res["las"]["gps_time_primer_trozo"] = tr0
        p.n_antes = res["las"]["n_puntos"]
        p.info("LAS: %s", res["las"])

        # --- trayectoria
        tray, info = cargar_trayectoria(cfg, cache_dir=out)
        dt = np.diff(tray.t)
        res["trayectoria"] = dict(fuente=cfg["trayectoria_modelo"]["fuente"], n=len(tray.t),
                                  t0=float(tray.t[0]), t1=float(tray.t[-1]), dt_mediana=float(np.median(dt)),
                                  huecos_mayores_tol=int((dt > cfg["tiempo"]["tolerancia_hueco_trayectoria_s"]).sum()),
                                  ajuste_vs_posi=info)
        res["trayectoria"]["gps_time_las_dentro"] = bool(tr0 and tr0[0] >= tray.t[0] and tr0[1] <= tray.t[-1])
        if cfg["rutas"].get("rtk_gga"):
            pct, nq = pct_rtk_fijo(cfg["rutas"]["rtk_gga"])
            res["rtk"] = dict(pct_fijo=pct, n_epocas=nq, umbral=cfg["qc"]["min_pct_rtk_fijo"],
                              aviso=bool(pct < cfg["qc"]["min_pct_rtk_fijo"]))
        else:
            res["rtk"] = dict(pct_fijo=None, aviso=True, motivo="no hay rtk.gga en los brutos")
        p.info("trayectoria: %s", res["trayectoria"])
        p.info("RTK: %s", res["rtk"])

        # --- cobertura bilateral por hilera (todas las del archivo vectorial filtrado)
        filas = cargar_filas(cfg, todas=True)
        lo, hi = cfg["qc"]["distancia_tray_eje_esperada_m"]
        filas_csv = []
        for _, f in filas.iterrows():
            m = MarcoHilera.de_linea(f.geometry)
            s, t = m.st(tray.p[:, :2])
            fila = dict(id_fila=f["id_fila"], **valores_atributos(cfg, f), largo_m=round(m.largo, 2))
            for nom, sg in (("izq", 1), ("der", -1)):
                ok = (np.abs(t) >= lo) & (np.abs(t) <= hi + 1.0) & (np.sign(t) == sg) & (s > 0) & (s < m.largo)
                celdas = np.unique(np.floor(s[ok]).astype(int))
                fila["cobertura_" + nom + "_pct"] = round(100 * len(celdas) / int(np.ceil(m.largo)), 1)
                fila["t_med_" + nom] = round(float(np.median(t[ok])), 2) if ok.any() else None
            fila["bilateral"] = fila["cobertura_izq_pct"] > 50 and fila["cobertura_der_pct"] > 50
            filas_csv.append(fila)
        df = pd.DataFrame(filas_csv)
        df.to_csv(os.path.join(out, "qc_cobertura_bilateral.csv"), index=False)
        p.info("cobertura bilateral:\n%s", df.to_string(index=False))
        res["n_filas"] = len(df)
        res["filas_bilaterales"] = int(df["bilateral"].sum())
        res["filas_prueba_bilaterales"] = df[df.id_fila.isin(cfg["recorte"]["filas_prueba_id"])].to_dict("records")

        # --- pasos del recorte
        fp = cargar_filas(cfg)
        vent, pasos = ventanas_tiempo(cfg, tray, fp)
        res["recorte_ventanas_tiempo"] = vent
        res["recorte_pasos"] = pasos
        p.info("ventanas de tiempo del recorte: %s", vent)
        guardar_json(res, os.path.join(out, "qc_resumen.json"))

        # --- mapa
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(8, 8), dpi=110)
            x0, y0, x1, y1 = filas.total_bounds          # ventana del mapa = extensión del archivo vectorial + margen
            mg = cfg["qc"].get("margen_mapa_m", 10.0)
            m_ = (tray.p[:, 0] > x0 - mg) & (tray.p[:, 0] < x1 + mg) & (tray.p[:, 1] > y0 - mg) & (tray.p[:, 1] < y1 + mg)
            ax.plot(tray.p[m_, 0], tray.p[m_, 1], ".", ms=0.5, color="tab:blue", label="trayectoria")
            for _, f in filas.iterrows():
                x, y = f.geometry.xy
                ax.plot(x, y, "g-", lw=1)
                ax.text(x[0], y[0] - 2, f["id_fila"], fontsize=7)
            pol = poligono_recorte(cfg, fp)
            xx, yy = pol.exterior.xy if pol.geom_type == "Polygon" else pol.geoms[0].exterior.xy
            ax.plot(xx, yy, "r-", lw=1.2, label="recorte")
            ax.set_aspect("equal")
            ax.legend()
            ax.set_title("QC %s: trayectoria (hf), hileras y recorte" % cfg["proyecto"]["nombre"])
            fig.savefig(os.path.join(out, "mapa_trayectoria_filas.png"))
            plt.close(fig)
        except Exception as e:  # pragma: no cover
            p.info("PNG omitido: %r", e)
        p.n_despues = p.n_antes
        p.salida = out
        p.extra["filas_bilaterales"] = f'{res["filas_bilaterales"]}/{res["n_filas"]}'
    return res


if __name__ == "__main__":
    main()
