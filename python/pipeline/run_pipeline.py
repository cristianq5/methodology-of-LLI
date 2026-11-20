# -*- coding: utf-8 -*-
"""Lanza el flujo completo o desde un paso.

  python run_pipeline.py --ensayo ficha.yaml                     # 01 .. 11 (modo de la ficha)
  python run_pipeline.py --ensayo ficha.yaml --desde 08          # 08 .. 11
  python run_pipeline.py --ensayo ficha.yaml --solo 09           # solo el paso 09
  python run_pipeline.py --ensayo ficha.yaml --limpiar-cache     # borra pcap_ventanas.npz, extrinseca.json y la caché por calle
  python run_pipeline.py --ensayo ficha.yaml --modo completo --bloque 5 --desde 02   # solo la hilera 5 entera
  python run_pipeline.py --ensayo ficha.yaml --modo completo     # todas las hileras, ejecucion.n_procesos en paralelo
  python run_pipeline.py --ensayo ficha.yaml --procesos 2        # sobrescribe ejecucion.n_procesos

"bloque" = unidad de trabajo del pipeline (NO el bloque del ensayo): en modo completo, una hilera.
Modo "prueba" (ejecucion_ensayo.modo): un único bloque con las hileras de prueba.
Modo "completo": 01 y 02 una vez (02 reparte el LAS en un bloque por hilera) y 03..11 hilera a hilera, con
memoria acotada; los bloques se reparten entre `ejecucion.n_procesos` procesos (por defecto 3). La calibración de
cada calle se hace una sola vez y se comparte entre bloques (caché por calle con huella y .lock). Al final se
consolidan los CSV en salidas\\11_exportar\\consolidado. Se registra la RAM pico (suma de procesos) en
logs\\ram_pico.json. NO lanzar "completo" sin la QA del recorte.
"""
import utils  # noqa: F401 (primero: entorno)
import argparse
import glob
import importlib
import json
import os
import threading
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp

import pandas as pd

PASOS = ["01_qc", "02_recorte", "03_filtrado_sor", "04_suelo", "05_dtm", "06_normalizacion",
         "07_chm", "08_segmentacion", "09_metricas", "10_lli", "11_exportar"]
GLOBALES = ("01_qc", "02_recorte")       # se ejecutan una vez para todos los bloques


def limpiar_cache(cfg):
    n = 0
    base = os.path.join(cfg["rutas"]["salida"], "09_metricas")
    for pat in ("pcap_ventanas.npz", "extrinseca.json", os.path.join("_calles", "*.json"), os.path.join("_calles", "*.lock"),
                os.path.join("_calles", "*_grupos.pkl"), os.path.join("_global", "montaje_bloque.json")):
        for f in glob.glob(os.path.join(base, "**", pat), recursive=True) + glob.glob(os.path.join(base, pat)):
            try:
                os.remove(f)
                n += 1
            except OSError:
                pass
    print("cachés borradas:", n)


def consolidar(cfg, bloques):
    out = os.path.join(utils.dir_salida(cfg, "11_exportar"), "consolidado")
    os.makedirs(out, exist_ok=True)
    for nom in ("lli_hileras.csv", "lli_tramos.csv", "lli_rebanadas.csv", "perfil_anchura.csv"):
        dfs = []
        for b in bloques:
            r = os.path.join(utils.dir_salida(utils.cfg_bloque(cfg, b), "11_exportar"), nom)
            if os.path.exists(r):
                dfs.append(pd.read_csv(r, dtype={"id_fila": str}))
        if dfs:
            pd.concat(dfs, ignore_index=True).to_csv(os.path.join(out, nom), index=False)
    print("consolidado en", out)


def resumen_final(cfg, bloques):
    """Resumen por hilera (completa / incompleta y motivo). Nunca silencioso."""
    filas = []
    for b in bloques:
        r = os.path.join(utils.dir_salida(utils.cfg_bloque(cfg, b), "11_exportar"), "lli_hileras.csv")
        if not os.path.exists(r):
            filas.append(dict(id_fila=b, estado_P="NO PROCESADA (falta lli_hileras.csv)"))
            continue
        filas.extend(pd.read_csv(r, dtype={"id_fila": str}).to_dict("records"))
    df = pd.DataFrame(filas)
    cols = [c for c in ("id_fila", *utils.atributos_croquis(cfg), "n_ventanas_validas", "n_rebanadas_con_P", "n_rebanadas", "P_media", "LLI_estricta_media",
                        "LLI_percentil_media", "LLI_estricta_mejor_media", "origen_P_mejor", "estado_P") if c in df]
    print("=" * 70)
    print("RESUMEN FINAL POR HILERA")
    print(df[cols].to_string(index=False))
    inc = df[~df["estado_P"].astype(str).str.startswith("completa")]
    print("Hileras completas: %d/%d" % (len(df) - len(inc), len(df)))
    for _, r in inc.iterrows():
        print("  INCOMPLETA hilera %s: %s" % (r["id_fila"], r["estado_P"]))
    out = os.path.join(utils.dir_salida(cfg, "11_exportar"), "consolidado")
    os.makedirs(out, exist_ok=True)
    df.to_csv(os.path.join(out, "resumen_hileras.csv"), index=False)
    return df


def _trabajo_bloque(cfg, b, pasos):
    """Ejecuta los pasos por bloque (03..11) de una hilera. Se llama en el proceso principal o en un worker."""
    cb = utils.cfg_bloque(cfg, b)
    t0 = time.time()
    tiempos = {}
    for nom in pasos:
        t1 = time.time()
        print("=" * 70, "\nPASO", nom, "| bloque", b if b is not None else "prueba", flush=True)
        mod, _, fase = nom.partition(":")             # "09_metricas:calibrar" = paso 09 solo en su fase de calibración
        if fase:
            importlib.import_module(mod).main(cb, fase=fase)
        else:
            importlib.import_module(mod).main(cb)
        tiempos[nom] = round(time.time() - t1, 1)
        print("  -> %s [%s]: %.1f s" % (nom, b, tiempos[nom]), flush=True)
    return b, round(time.time() - t0, 1), tiempos


class MonitorRAM(threading.Thread):
    """Muestrea la memoria residente del proceso y de sus hijos; guarda el pico."""

    def __init__(self, cada=2.0):
        super().__init__(daemon=True)
        import psutil
        self.ps = psutil.Process()
        self.cada, self.pico, self.parar = cada, 0, False

    def run(self):
        while not self.parar:
            try:
                tot = self.ps.memory_info().rss + sum(c.memory_info().rss for c in self.ps.children(recursive=True))
                self.pico = max(self.pico, tot)
            except Exception:
                pass
            time.sleep(self.cada)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--desde", default="01", help="número del primer paso (p. ej. 05)")
    ap.add_argument("--hasta", default="11")
    ap.add_argument("--solo", default=None)
    ap.add_argument("--ensayo", default=None, help="ficha del ensayo (ver config/ensayo_plantilla.yaml); "
                    "por defecto, la variable de entorno LLI_ENSAYO")
    ap.add_argument("--parametros", default=None, help="parametros.yaml alternativo (por defecto config/parametros.yaml)")
    ap.add_argument("--limpiar-cache", action="store_true")
    ap.add_argument("--modo", default=None, choices=["prueba", "completo"], help="sobrescribe recorte.modo")
    ap.add_argument("--bloque", default=None, help="(modo completo) ids de hilera separados por comas")
    ap.add_argument("--procesos", type=int, default=None, help="sobrescribe ejecucion.n_procesos")
    a = ap.parse_args()
    cfg = utils.cargar_config(a.ensayo, a.parametros)
    # los pasos sueltos que lancen los procesos hijos usan la misma ficha
    os.environ[utils.VAR_ENTORNO_ENSAYO] = cfg["_ruta_ensayo"]
    if a.modo:
        cfg["recorte"]["modo"] = a.modo
    if a.bloque is not None:
        cfg["_bloques_sel"] = [x.strip() for x in str(a.bloque).split(",")]
    if a.limpiar_cache:
        limpiar_cache(cfg)
    sel = [p for p in PASOS if p.startswith(a.solo)] if a.solo else [p for p in PASOS if a.desde <= p[:2] <= a.hasta]
    bloques = utils.lista_bloques(cfg)
    n_proc = max(1, a.procesos or cfg.get("ejecucion", {}).get("n_procesos", 1))
    mon = MonitorRAM()
    mon.start()
    t0 = time.time()
    for nom in [p for p in sel if p in GLOBALES]:
        t1 = time.time()
        print("=" * 70, "\nPASO", nom, flush=True)
        importlib.import_module(nom).main(cfg)
        print("  -> %s: %.1f s" % (nom, time.time() - t1), flush=True)
    por_bloque = [p for p in sel if p not in GLOBALES]
    resumen = {}

    def ejecutar(pasos, etiqueta=""):
        if not pasos:
            return
        if n_proc > 1 and len(bloques) > 1:
            print("Ejecución en paralelo%s: %d procesos, %d bloques" % (etiqueta, n_proc, len(bloques)), flush=True)
            with ProcessPoolExecutor(max_workers=n_proc, mp_context=mp.get_context("spawn")) as ex:
                futs = [ex.submit(_trabajo_bloque, cfg, b, pasos) for b in bloques]
                for f in as_completed(futs):
                    b, dt, tt = f.result()
                    r = resumen.setdefault(str(b), dict(tiempo_s=0, pasos={}))
                    r["tiempo_s"] += dt; r["pasos"].update(tt)
        else:
            for b in bloques:
                b_, dt, tt = _trabajo_bloque(cfg, b, pasos)
                r = resumen.setdefault(str(b), dict(tiempo_s=0, pasos={}))
                r["tiempo_s"] += dt; r["pasos"].update(tt)

    if "09_metricas" in por_bloque and cfg["pcap"].get("modo_calibracion", "por_calle") == "global_bloque":
        # Calibración global del bloque: (1) pasos previos + 09 en fase "calibrar" (calibración por calle y grupos) en
        # todas las hileras; (2) ajuste conjunto de un único montaje con todas las calles; (3) 09 (P con ese montaje) .. 11.
        ejecutar([p for p in por_bloque if p < "09"] + ["09_metricas:calibrar"], " (fase 1: hasta la calibración por calle)")
        t1 = time.time()
        print("=" * 70, "\nPASO 09b_calibracion_global (montaje único del bloque)", flush=True)
        importlib.import_module("calibracion_global").main(cfg)
        print("  -> 09b_calibracion_global: %.1f s" % (time.time() - t1), flush=True)
        resumen["09b_calibracion_global"] = dict(tiempo_s=round(time.time() - t1, 1))
        ejecutar([p for p in por_bloque if p >= "09"], " (fase 2: P con el montaje global, LLI y exportación)")
    else:
        ejecutar(por_bloque)
    if bloques != [None] and "11_exportar" in sel:
        consolidar(cfg, bloques)
        resumen_final(cfg, bloques)
    mon.parar = True
    total = time.time() - t0
    info = dict(tiempo_total_s=round(total, 1), ram_pico_gb=round(mon.pico / 2 ** 30, 2), n_procesos=n_proc, bloques=resumen,
                modo=cfg["recorte"]["modo"], pasos=sel)
    with open(os.path.join(utils.dir_logs(cfg), "ram_pico.json"), "w", encoding="utf-8") as f:
        json.dump(info, f, indent=2, ensure_ascii=False)
    print("TOTAL %.1f s | RAM pico %.2f GB" % (total, mon.pico / 2 ** 30))


if __name__ == "__main__":
    main()
