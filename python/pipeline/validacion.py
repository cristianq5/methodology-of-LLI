# -*- coding: utf-8 -*-
"""QA-22: validación de la extrínseca de cada calle con datos que NO se usaron en el ajuste (hold-out).

El ajuste (calibrar) ya tiene muchos grados de libertad (nudos de fase y de dW, deriva); medir el residuo con los mismos
grupos es optimista. Aquí se sortean otros grupos temporales (otra semilla) EXCLUYENDO todos los que intervinieron en el
ajuste (k_ajuste, que guarda la calibración), solo en tramos estables, y se evalúa el MISMO criterio
(mediana < residuo_mediana_max_m y p90 < residuo_p90_max_m, sin cambiar el umbral).
Una calibración que no registre k_ajuste (versión anterior) se declara NO válida: hay que recalibrarla.
Nota: los grupos de validación están intercalados en el tiempo con los del ajuste (sellos de 1/1024 s distintos).
"""
import numpy as np

import extrinseca as ex


def holdout(cfg, tray, inf, w, las_t, las_xyz, n_grupos=600, semilla=777):
    """Devuelve 'inf' con residuos de ajuste conservados (ajuste_*) y el criterio de validez evaluado en hold-out."""
    ce = cfg["pcap"]["extrinseca"]
    par = ex.par_desde_informe(inf)
    om_max = ce.get("omega_max_deg_s", 0)
    if om_max:
        est = ex.tasa_giro(tray, las_t) <= om_max
        las_t, las_xyz = las_t[est], las_xyz[est]
    inf = dict(inf)
    inf.setdefault("ajuste_mediana_m", inf["residuo_mediana_m"])
    inf.setdefault("ajuste_p90_m", inf["residuo_p90_m"])
    if "k_ajuste" not in inf:
        # calibración antigua que no registró sus grupos: no se puede garantizar que el hold-out no los reutilice
        inf.update(holdout_mediana_m=None, holdout_p90_m=None, valida=False, holdout_n=0,
                   holdout_motivo="la calibración no registra sus grupos (k_ajuste): recalibrar (--limpiar-cache)")
        return inf
    dist_m = w["dist_raw"].astype(np.float64) * cfg["pcap"]["unidad_distancia_m"]
    # grupos sorteados entre los que NO intervinieron en el ajuste (excluir = k_ajuste)
    G = ex.grupos_temporales(las_t, las_xyz, w["t"], w["az"], dist_m, n_grupos, np.random.default_rng(semilla),
                             excluir=inf["k_ajuste"])
    if G is None:
        inf.update(holdout_mediana_m=None, holdout_p90_m=None, valida=False, holdout_n=0,
                   holdout_motivo="sin grupos temporales disponibles fuera del ajuste")
        return inf
    assert not set(int(k) for k in G["k"]) & set(int(k) for k in inf["k_ajuste"]), "hold-out con grupos del ajuste"
    d = ex.residuos_grupos(G, tray, par, 2.0)
    med, p90 = float(np.median(d)), float(np.percentile(d, 90))
    inf["valida_ajuste"] = bool(inf["ajuste_mediana_m"] < ce["residuo_mediana_max_m"] and inf["ajuste_p90_m"] < ce["residuo_p90_max_m"])
    inf.update(holdout_mediana_m=med, holdout_p90_m=p90, holdout_n=int(len(d)),
               residuo_mediana_m=med, residuo_p90_m=p90,
               valida=bool(med < ce["residuo_mediana_max_m"] and p90 < ce["residuo_p90_max_m"]))
    return inf
