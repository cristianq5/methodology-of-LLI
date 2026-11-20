# -*- coding: utf-8 -*-
"""10_lli: LLI = MH * MW * (1 - P) por rebanada; filtro de +-2,5 DE por hilera; media por hilera y por
tramos de 1,5 m. Se calcula con MH/MW estrictos y con percentil, y con P de cada vía (principal, pcap, rayos).
QA de signo: L = 1 - P debe correlacionar positivamente con la densidad de puntos de copa."""
import utils  # noqa: F401
import os
import numpy as np
import pandas as pd

from utils import cargar_config, dir_salida, Paso, guardar_json, png_lineas, atributos_croquis
import lli_core as lc

NOMBRE = "10_lli"

# Métricas de copa de 09 (lli_core.metricas_copa): se agregan sin filtro de atípicos (las rebanadas sin copa ya son NaN).
COPA_MEDIA = ("area_seccion_m2", "area_seccion_p_m2", "h_base_copa_m", "prof_copa_m", "anchura_max_m", "h_anchura_max_m",
              "TRV_m3_ha", "TRV_p_m3_ha")
COPA_SUMA = ("vol_copa_m3", "vol_copa_p_m3")


def agregar_copa(g, cfg):
    """Medias (_media) y volúmenes totales (_total, m³, y _por_m, m³ por metro de hilera) de un grupo de rebanadas."""
    f = {}
    for c_ in COPA_MEDIA:
        if c_ in g:
            f[c_ + "_media"] = g[c_].mean()
    largo = len(g) * cfg["metricas"]["rebanada_m"]
    for c_ in COPA_SUMA:
        if c_ in g:
            f[c_ + "_total"] = g[c_].sum(min_count=1)
            f[c_ + "_por_m"] = f[c_ + "_total"] / largo if largo else np.nan
    if "dist_hileras_m" in g:
        f["dist_hileras_m"] = g["dist_hileras_m"].iloc[0]
    return f


def main(cfg=None):
    cfg = cfg or cargar_config()
    out = dir_salida(cfg, NOMBRE)
    c = cfg["lli"]
    ent = os.path.join(dir_salida(cfg, "09_metricas"), "metricas_rebanadas.csv")
    with Paso(cfg, NOMBRE, entrada=ent, parametros=c) as p:
        df = pd.read_csv(ent, dtype={"id_fila": str})
        atr = [a for a in atributos_croquis(cfg) if a in df.columns]   # atributos del archivo vectorial (variedad, bloque...)
        p.n_antes = len(df)
        variantes = {"estricta": ("MH_max", "MW_max"), "percentil": ("MH_p", "MW_p")}
        fuentes = ["P"] + [c_ for c_ in ("P_h06", "P_mejor") if c_ in df.columns] + [f"P_{v}" for v in ("pcap", "rayos", "pcap_todas") if f"P_{v}" in df.columns]
        cols = []
        for vn, (mh, mw) in variantes.items():
            for pf in fuentes:
                nom = f"LLI_{vn}" if pf == "P" else f"LLI_{vn}_{pf[2:]}"
                df[nom] = lc.lli(df[mh], df[mw], df[pf])
                cols.append(nom)
        df["L"] = 1.0 - df["P"]
        # --- filtro de atípicos (por hilera) y agregación
        filtradas = []
        for k, g in df.groupby("id_fila", sort=False):
            g = g.copy()
            for col in cols:
                ok = lc.filtro_atipicos(g[col].values, c["filtro_atipicos_de"])
                g[col + "_ok"] = ok
            filtradas.append(g)
        df = pd.concat(filtradas)
        df["tramo"] = (np.floor((df["s_ini"] - df.groupby("id_fila")["s_ini"].transform("min") + 1e-9) / cfg["hilera"]["tramo_m"])).astype(int)
        agg_h, agg_t = [], []
        for k, g in df.groupby("id_fila", sort=False):
            n_p = int(g["P"].notna().sum())
            motivos = g.loc[g["P"].isna(), "P_motivo_nan"].fillna("").value_counts() if "P_motivo_nan" in g else []
            if n_p == len(g):
                estado = "completa"
            elif n_p == 0:
                estado = "SIN P: " + (motivos.index[0] if len(motivos) else "motivo no registrado")
            else:
                estado = "parcial: %d/%d rebanadas con P; sin P por: %s" % (n_p, len(g), "; ".join("%s (%d)" % (i, v) for i, v in motivos.items()))
            fila = dict(id_fila=k, **{a: g[a].iloc[0] for a in atr}, n_rebanadas=len(g), n_rebanadas_con_P=n_p, estado_P=estado,
                        n_ventanas_validas=int(g["n_ventanas_validas"].iloc[0]) if "n_ventanas_validas" in g else None,
                        n_rebanadas_P_mejor=int(g["P_mejor"].notna().sum()) if "P_mejor" in g else None,
                        origen_P_mejor="; ".join("%s=%d" % (a_, b_) for a_, b_ in g["via_P_mejor"].value_counts().items() if a_) if "via_P_mejor" in g else None,
                        P_media=g["P"].mean(), P_h06_media=g["P_h06"].mean(),
                        h_max_P_izq_media=g["h_max_P_izq"].mean(), h_max_P_der_media=g["h_max_P_der"].mean(),
                        h_max_P_izq_mediana=g["h_max_P_izq"].median(), h_max_P_der_mediana=g["h_max_P_der"].median(),
                        cob_MH_media=g["cob_MH"].mean(),
                        MH_max_media=g.MH_max.mean(), MW_max_media=g.MW_max.mean(), MH_p_media=g.MH_p.mean(), MW_p_media=g.MW_p.mean(),
                        **agregar_copa(g, cfg))
            for col in cols:
                v = g.loc[g[col + "_ok"], col]
                fila[col + "_media"] = v.mean()
                fila[col + "_n"] = int(len(v))
                fila[col + "_descartadas"] = int(g[col].notna().sum() - len(v))
            agg_h.append(fila)
            for tr, gt in g.groupby("tramo"):
                if len(gt) < c["min_rebanadas_tramo"]:      # tramo incompleto
                    continue
                f2 = dict(id_fila=k, **{a: g[a].iloc[0] for a in atr}, tramo=tr, s_ini=gt.s_ini.min(), s_fin=gt.s_ini.max() + cfg["metricas"]["rebanada_m"], n_rebanadas=len(gt), P_media=gt["P"].mean(), P_h06_media=gt["P_h06"].mean(),
                      h_max_P_izq_media=gt["h_max_P_izq"].mean(), h_max_P_der_media=gt["h_max_P_der"].mean(), cob_MH_media=gt["cob_MH"].mean(),
                      **agregar_copa(gt, cfg))
                for col in cols:
                    f2[col] = gt.loc[gt[col + "_ok"], col].mean()
                agg_t.append(f2)
        H, T = pd.DataFrame(agg_h), pd.DataFrame(agg_t)
        df.to_csv(os.path.join(out, "lli_rebanadas.csv"), index=False)
        H.to_csv(os.path.join(out, "lli_hileras.csv"), index=False)
        T.to_csv(os.path.join(out, "lli_tramos.csv"), index=False)
        # --- resumen y QA
        res = dict(hileras=H.to_dict("records"))
        rng = {}
        for k, g in T.groupby("id_fila"):
            rng[k] = {col: [float(g[col].min()), float(g[col].max())] for col in ("LLI_estricta", "LLI_percentil", "LLI_estricta_h06") if col in g}
        res["rango_tramos_1_5m"] = rng
        lo, hi = c["rango_esperado_m2"]
        res["avisos_rango"] = [f"{r['id_fila']}: LLI_estricta {r['LLI_estricta_media']:.2f} fuera de [{lo},{hi}]" for r in agg_h
                               if not (lo <= r["LLI_estricta_media"] <= hi)]
        if c["qa_signo_L"]:
            # L = 1 - P debe crecer con la frondosidad: se contrasta con proxies independientes de P
            # (puntos de copa por rebanada, MH y MW: los huecos entre árboles bajan MH/MW y suben P)
            dens = df["n_pts"] / (df["MH_max"] * df["MW_max"])
            cor = {}
            for nom, x in (("densidad_puntos", dens), ("MH", df["MH_max"]), ("MW", df["MW_max"])):
                ok = np.isfinite(x) & np.isfinite(df["L"])
                cor[nom] = float(np.corrcoef(x[ok], df["L"][ok])[0, 1]) if ok.sum() > 5 else None
            vals = [v for v in cor.values() if v is not None]
            res["qa_signo_L"] = dict(corr_L=cor, ok=bool(vals and np.mean(vals) > 0))
        guardar_json(res, os.path.join(out, "lli_resumen.json"))
        p.info("LLI por hilera:\n%s", H[["id_fila", *atr, "P_media", "P_h06_media", "LLI_estricta_media", "LLI_estricta_h06_media", "LLI_percentil_media", "h_max_P_izq_mediana", "h_max_P_der_mediana", "cob_MH_media"]].to_string(index=False))
        p.info("rango por tramo 1,5 m: %s", rng)
        p.info("QA signo L: %s | avisos: %s", res.get("qa_signo_L"), res["avisos_rango"])
        for k, g in df.groupby("id_fila", sort=False):
            png_lineas([(g.s_centro, g["LLI_estricta"], "LLI estricta"), (g.s_centro, g["LLI_percentil"], "LLI percentil")],
                       os.path.join(out, f"lli_hilera{k}.png"), f"Hilera {k}: LLI por rebanada", "s (m)", "m2")
        p.n_despues = len(df)
        p.salida = out
    return H, T


if __name__ == "__main__":
    main()
