# -*- coding: utf-8 -*-
"""09_metricas: MH, MW y P por rebanada (0,1 m) y por lado de cada hilera.

  MH, MW : estrictos (máximo) y con percentil (config metricas.mh_percentil / mw_percentil_t).
  P      : vía principal 'pcap' (disparos decodificados, incluidos los de distancia 0, con la pose de la
           trayectoria y la extrínseca estimada por ventana frente al LAS). La vía 'rayos' (trazado hasta los
           puntos del LAS) se calcula como CONTRASTE de forma: está sesgada a la baja por construcción (QA-08).
  Bi = rayos que ATRAVIESAN la copa; P = Bi/Bt (porosidad); L = 1 - P. P usa SOLO las franjas medidas
  (>= min_rayos_celda rayos); no se extrapola. Se informa la cobertura vertical de P (QA-05).
  Oficial: franjas desde umbral_altura_m (0,3 m). Sensibilidad: P desde porosidad.umbral_P_sensibilidad_m (0,6 m).
Cachés (QA-02): pcap_ventanas.npz y extrinseca.json llevan una huella (ventanas, parámetros, entrada, código)
y se recalculan si no coincide. Ventanas con extrínseca fuera de criterio se excluyen de P (QA-06).
Métricas de copa (lli_core.metricas_copa, misma malla rebanada x franja): área de la sección, volumen, altura de inicio y
profundidad de copa, anchura máxima y su altura, TRV; y perfil de anchura por altura por hilera y por tramo (perfil_anchura.csv).
Salidas: metricas_rebanadas.csv, perfil_anchura.csv, conteos_rayos.npz, extrinseca.json, comparacion_P.json, PNG de control.
"""
import utils  # noqa: F401
import glob
import os
import time
import datetime
import numpy as np
import pandas as pd
import laspy

from utils import (cargar_config, dir_salida, Paso, MarcoHilera, guardar_json, cargar_json, png_lineas,
                   huella, huella_fichero, huella_codigo, valores_atributos)
from trayectoria import cargar_trayectoria
import pcap_hesai as ph
import extrinseca as ex
import validacion
import lli_core as lc
import porosidad as po
import calibracion_global as cglob

NOMBRE = "09_metricas"
VERSION_CACHE = 3          # 3: la calibración registra sus grupos (k_ajuste) y el hold-out los excluye


def decodificar_ventanas(cfg, vent, out, p):
    """Decodifica el pcap en las ventanas; caché validada por huella (ventanas + pcap + unidad)."""
    ruta = os.path.join(out, "pcap_ventanas.npz")
    h = huella(dict(v=VERSION_CACHE, vent=[list(map(float, w)) for w in vent], pcap=huella_fichero(cfg["pcap"]["ruta"])))
    if os.path.exists(ruta):
        z = np.load(ruta, allow_pickle=False)
        if "huella" in z.files and str(z["huella"]) == h:
            n = int(z["n"])
            p.info("caché pcap_ventanas.npz válida (huella %s, %d ventanas)", h, n)
            return [dict(t=z["t%d" % i], az=z["az%d" % i], dist_raw=z["d%d" % i]) for i in range(n)], h
        p.info("caché pcap_ventanas.npz NO coincide con las ventanas/entrada actuales: se recalcula")
    ws = []
    for a, b in vent:                     # una ventana cada vez: memoria acotada (~0,3 GB por ventana de 25 s)
        t0 = time.time()
        d = ph.leer_disparos(cfg["pcap"]["ruta"], a, b)
        ws.append(dict(t=d["t"], az=d["az_deg"], dist_raw=d["dist_raw"]))
        p.info("pcap decodificado [%.1f, %.1f] s: %d bloques en %.0f s", a, b, len(d["t"]), time.time() - t0)
    np.savez(ruta, n=len(ws), huella=h, **{**{"t%d" % i: w["t"] for i, w in enumerate(ws)},
                                            **{"az%d" % i: w["az"] for i, w in enumerate(ws)},
                                            **{"d%d" % i: w["dist_raw"] for i, w in enumerate(ws)}})
    return ws, h


def _leer_valida(ruta, h):
    if os.path.exists(ruta):
        try:
            inf = cargar_json(ruta)
            if inf.get("huella") == h:
                return inf
        except Exception:
            return None
    return None


def extrinseca_cache(cfg, tray, ventanas, lims, h_vent, tl, xyz, out, p):
    """Extrínseca por ventana (= un paso por una calle). QA-20: la caché es POR CALLE y global
    (salidas\\09_metricas\\_calles\\calle_<t0>_<t1>.json): la calibración de una calle se hace una sola vez y la
    reutilizan las dos hileras vecinas (aunque se procesen en paralelo: un fichero .lock evita calibrar dos veces).
    Cada caché lleva huella (ventana, parámetros, unidad, signo, pcap, código, semilla) y se valida al cargar (QA-02)."""
    calles = os.path.join(cfg["rutas"]["salida"], NOMBRE, "_calles")
    os.makedirs(calles, exist_ok=True)
    codigo = huella_codigo("extrinseca.py", "pcap_hesai.py")
    pars, infs = [], []

    def semillas_validas(excluir):
        """(rotvec_e [rad], omega) de las calles ya calibradas y válidas (cualquier versión de la caché)."""
        out_ = []
        for f in sorted(glob.glob(os.path.join(calles, "calle_*.json"))):
            if os.path.abspath(f) == os.path.abspath(excluir):
                continue
            try:
                j = cargar_json(f)
                if j.get("valida"):
                    out_.append(dict(rotvec_e=np.radians(j["rotvec_e_deg"]).tolist(), omega=float(np.mean(j["omega"]))))
            except Exception:
                pass
        return out_
    for i, (w, (a, b)) in enumerate(zip(ventanas, lims)):
        h = huella(dict(v=VERSION_CACHE, win=[round(float(a), 2), round(float(b), 2)], ce=cfg["pcap"]["extrinseca"],
                        unidad=cfg["pcap"]["unidad_distancia_m"], signo=cfg["pcap"]["signo_azimut"], codigo=codigo,
                        pcap=huella_fichero(cfg["pcap"]["ruta"]), semilla=cfg["proyecto"]["semilla_aleatoria"]))
        ruta = os.path.join(calles, "calle_%d_%d.json" % (round(a), round(b)))
        usar = cfg["pcap"]["extrinseca"].get("reutilizar", True)
        inf = _leer_valida(ruta, h) if usar else None
        if inf is not None:
            p.info("ventana %d: calibración de la calle reutilizada (%s, huella %s)", i, os.path.basename(ruta), h)
        else:
            lock = ruta + ".lock"
            propietario = True
            try:
                os.close(os.open(lock, os.O_CREAT | os.O_EXCL))
            except FileExistsError:
                propietario = False
            if not propietario:                      # otro proceso la está calibrando: se espera (máx. 2 h)
                p.info("ventana %d: otra hilera calibra esta calle; se espera", i)
                t_w = time.time()
                while time.time() - t_w < 7200 and _leer_valida(ruta, h) is None and os.path.exists(lock):
                    time.sleep(10)
                inf = _leer_valida(ruta, h)
            if inf is None:
                t0 = time.time()
                p.info("ventana %d [%.1f, %.1f]: calibrando la calle", i, a, b)
                m = (tl >= a) & (tl <= b)
                sem = semillas_validas(ruta)
                p.info("ventana %d: %d semillas de calles válidas", i, len(sem))
                par, inf = ex.calibrar_por_ventana(cfg, tray, [w], [(a, b)], tl[m], xyz[m], log=p.info, semillas=sem)
                inf = dict(inf[0], huella=h, fecha=datetime.datetime.now().isoformat(timespec="seconds"),
                           tiempo_s=round(time.time() - t0), codigo=codigo, ventana_global=[float(a), float(b)])
                guardar_json(inf, ruta)
                if propietario and os.path.exists(lock):
                    os.remove(lock)
        if "holdout_mediana_m" not in inf:           # QA-22: validez con datos no usados en el ajuste (se guarda en la caché de la calle)
            mh_ = (tl >= a) & (tl <= b)
            inf = validacion.holdout(cfg, tray, inf, w, tl[mh_], xyz[mh_])
            guardar_json(inf, ruta)
            p.info("ventana %d: validación hold-out: mediana %s m, p90 %s m -> %s (ajuste: %.3f / %.3f)", i,
                   inf["holdout_mediana_m"], inf["holdout_p90_m"], "VALIDA" if inf["valida"] else "NO VALIDA",
                   inf["ajuste_mediana_m"], inf["ajuste_p90_m"])
        if modo_global(cfg) and not cglob.grupos_validos(ruta, inf):   # grupos para el ajuste global del bloque
            mg_ = (tl >= a) & (tl <= b)
            if cglob.guardar_grupos(cfg, ruta, inf, w, tl[mg_], xyz[mg_], cfg.get("_bloque")):
                p.info("ventana %d: grupos de ajuste y validación guardados para la calibración global", i)
        infs.append(inf)
        pars.append(ex.par_desde_informe(inf))
    resumen = dict(ventanas=infs, residuo_mediana_m=max(i["residuo_mediana_m"] for i in infs),
                   residuo_p90_m=max(i["residuo_p90_m"] for i in infs),
                   residuo_rms_inliers_m=max((i["residuo_rms_inliers_m"] or 9) for i in infs),
                   frac_menor_0_15=min(i["frac_menor_0_15"] for i in infs),
                   n_ventanas_validas=int(sum(i["valida"] for i in infs)))
    guardar_json(resumen, os.path.join(out, "extrinseca.json"))
    return pars, resumen


_CACHE_LINEAS = {}


def distancia_hileras(cfg, d, p):
    """Distancia entre hileras (m) para el TRV. metricas.distancia_entre_hileras_m si está en la configuración; si es null,
    se mide en el archivo vectorial: distancia del punto medio de la línea de esta hilera a la línea más próxima de OTRA
    hilera (todas las líneas del fichero, sin el filtro del ensayo, para no saltarse vecinas de otro tratamiento).
    Devuelve (distancia, origen)."""
    v = cfg["metricas"].get("distancia_entre_hileras_m")
    if v:
        return float(v), "config"
    import geopandas as gpd
    from shapely.geometry import Point
    clave = (cfg["rutas"]["filas_gpkg"], cfg["rutas"]["filas_capa"], cfg["crs"]["epsg"])
    if clave not in _CACHE_LINEAS:
        g = gpd.read_file(clave[0], layer=clave[1]) if clave[1] else gpd.read_file(clave[0])
        if g.crs is not None and g.crs.to_epsg() != clave[2]:
            g = g.to_crs(clave[2])
        _CACHE_LINEAS[clave] = g.geometry
    m0 = MarcoHilera.de_dict(d["marco_gpkg"])
    dist = _CACHE_LINEAS[clave].distance(Point(m0.xy(0.5 * m0.largo, 0.0)))
    dist = dist[dist > cfg["metricas"]["distancia_min_otra_hilera_m"]]       # fuera la propia línea
    if dist.empty:
        p.info("AVISO: no hay otra hilera en el archivo vectorial: TRV sin calcular (indica metricas.distancia_entre_hileras_m)")
        return float("nan"), "sin_hilera_vecina"
    return float(dist.min()), "vectorial"


def modo_global(cfg):
    m = cfg["pcap"].get("modo_calibracion", "por_calle")
    if m not in ("por_calle", "global_bloque"):
        raise ValueError("pcap.modo_calibracion debe ser 'por_calle' o 'global_bloque'")
    return m == "global_bloque"


def extrinseca_global(cfg, lims, p):
    """Extrínseca de cada ventana a partir del montaje global del bloque (calibracion_global.py). La validez de cada
    ventana es la de su calle con el montaje global (hold-out, mismo umbral)."""
    mont = cglob.cargar_montaje(cfg)
    idx = {c["calle"]: i for i, c in enumerate(mont["calles"])}
    pars, infs = [], []
    for i, (a, b) in enumerate(lims):
        nom = "calle_%d_%d.json" % (round(a), round(b))
        if nom not in idx:
            raise RuntimeError("La ventana %d (%s) no está en el montaje global del bloque: repite la calibración global" % (i, nom))
        c = mont["calles"][idx[nom]]
        pars.append(cglob.par_ventana(mont, idx[nom], c["valida"]))
        infs.append(dict(c, ventana_t=[float(a), float(b)], residuo_mediana_m=c["holdout_mediana_m"],
                         residuo_p90_m=c["holdout_p90_m"]))
        p.info("ventana %d: %s con el montaje global (hold-out %.3f / %.3f)", i, nom, c["holdout_mediana_m"], c["holdout_p90_m"])
    resumen = dict(modo_calibracion="global_bloque", montaje=dict(huella=mont["huella"], lectura=mont["lectura"],
                   elevaciones=mont["elevaciones"], calle_referencia=mont["calle_referencia"], comun=mont["comun"],
                   n_validas_bloque="%d/%d" % (mont["n_validas"], len(mont["calles"]))),
                   ventanas=infs, residuo_mediana_m=max(i["residuo_mediana_m"] for i in infs),
                   residuo_p90_m=max(i["residuo_p90_m"] for i in infs),
                   residuo_rms_inliers_m=max((i["residuo_rms_inliers_m"] or 9) for i in infs),
                   frac_menor_0_15=min(i["frac_menor_0_15"] for i in infs),
                   n_ventanas_validas=int(sum(i["valida"] for i in infs)))
    return pars, resumen


def main(cfg=None, fase="todo"):
    """fase = "todo" (MH, MW, P) | "calibrar" (solo pcap + calibración por calle + grupos; la usa run_pipeline en el modo
    global_bloque antes del ajuste conjunto del bloque)."""
    cfg = cfg or cargar_config()
    out = dir_salida(cfg, NOMBRE)
    seg_dir = dir_salida(cfg, "08_segmentacion")
    m_ = cfg["metricas"]
    pc = cfg["porosidad"]
    with Paso(cfg, NOMBRE, entrada=seg_dir, parametros=dict(metricas=m_, porosidad=pc, pcap_extrinseca=cfg["pcap"]["extrinseca"])) as p:
        seg = np.load(os.path.join(seg_dir, "hileras.npz"))
        info_seg = cargar_json(os.path.join(seg_dir, "marcos.json"))
        meta = cargar_json(os.path.join(dir_salida(cfg, "02_recorte"), "recorte_meta.json"))
        tray, _ = cargar_trayectoria(cfg, cache_dir=dir_salida(cfg, "01_qc"))
        dtm = po.DTMLookup(os.path.join(dir_salida(cfg, "05_dtm"), "dtm.tif"))
        u = cfg["normalizacion"]["umbral_altura_m"]
        u06 = pc["umbral_P_sensibilidad_m"]
        j06 = int(round((u06 - u) / m_["franja_m"]))
        ids = list(info_seg["marcos"].keys())
        hmax = max(float(seg[f"{k}__hag"].max()) for k in ids)
        nf = int(np.ceil((hmax - u) / m_["franja_m"])) + 1
        # ---------------- MH, MW y envolvente por hilera
        hileras = []
        for k in ids:
            d = info_seg["marcos"][k]
            marco = MarcoHilera.de_dict(d["marco"])
            sl, sh = d["s_rango"]
            n_reb = int(np.floor((sh - sl) / m_["rebanada_m"] + 1e-9))
            s, t, h = seg[f"{k}__s"], seg[f"{k}__t"], seg[f"{k}__hag"]
            perfil = lc.perfiles_rebanada(s, t, h, sl, n_reb, m_["rebanada_m"], m_["franja_m"], u,
                                          m_["min_puntos_franja"], m_["min_puntos_rebanada"],
                                          m_["mh_percentil"], tuple(m_["mw_percentil_t"]), n_franjas=nf,
                                          env_percentil_t=tuple(pc["envolvente_percentil_t"]))
            env = lc.envolvente(perfil, "max" if pc["envolvente_estadistico"] == "max" else "percentil",
                                relleno=bool(pc.get("envolvente_relleno", False)))
            dist_h, origen_dist = distancia_hileras(cfg, d, p)
            hileras.append(dict(id=k, marco=marco, s_ini=sl, s_fin=sh, perfil=perfil, env=env,
                                dtm_fn=po.dtm_fn_hilera(dtm, marco), atributos=valores_atributos(cfg, d),
                                copa=lc.metricas_copa(perfil, dist_h), dist_hileras=dist_h, origen_dist=origen_dist))
            p.info("hilera %s: %d rebanadas, MH_max medio %.2f, MW_max medio %.2f; distancia entre hileras %.2f m (%s)", k, n_reb,
                   np.nanmean(perfil.MH_max), np.nanmean(perfil.MW_max), dist_h, origen_dist)
        # ---------------- datos LAS (recorte normalizado, sin ruido)
        las = laspy.read(os.path.join(dir_salida(cfg, "06_normalizacion"), "recorte_hag.las"))
        ok = np.asarray(las.classification) != 7
        xyz = np.c_[np.asarray(las.x), np.asarray(las.y), np.asarray(las.z)][ok]
        tl = np.asarray(las.gps_time)[ok]
        del las
        p.n_antes = int(len(tl))
        # ---------------- pcap + extrínseca por ventana (con huellas)
        metodo = pc["metodo"]
        lims = [tuple(v) for v in meta["ventanas_tiempo"]]
        ventanas, h_vent = decodificar_ventanas(cfg, lims, out, p)
        assert len(ventanas) == len(lims)
        pars, inf_ext = extrinseca_cache(cfg, tray, ventanas, lims, h_vent, tl, xyz, out, p)
        if fase == "calibrar":
            p.info("fase 'calibrar' terminada: %d ventanas calibradas por calle", len(pars))
            p.salida = os.path.join(cfg["rutas"]["salida"], NOMBRE, "_calles")
            return None
        if modo_global(cfg):                     # P con el montaje global del bloque y la lectura de orientación asociada
            pars, inf_ext = extrinseca_global(cfg, lims, p)
            tray = cglob.trayectoria_global(cfg)
            guardar_json(inf_ext, os.path.join(out, "extrinseca.json"))
        ce = cfg["pcap"]["extrinseca"]
        for i, iw in enumerate(inf_ext["ventanas"]):
            p.info("ventana %d [%.1f, %.1f]: residuo mediana %.3f m, p90 %.3f m, frac<0,15 %.2f -> %s", i, *iw["ventana_t"],
                   iw["residuo_mediana_m"], iw["residuo_p90_m"], iw["frac_menor_0_15"],
                   "VALIDA" if iw["valida"] else "NO VALIDA: excluida de P")
        n_val = int(sum(iw["valida"] for iw in inf_ext["ventanas"]))
        pcap_valida = n_val > 0
        p.info("criterio de validez: mediana < %.2f m Y p90 < %.2f m; ventanas válidas %d/%d", ce["residuo_mediana_max_m"],
               ce["residuo_p90_max_m"], n_val, len(pars))
        p.extra["ventanas_validas"] = "%d/%d" % (n_val, len(pars))
        p.extra["residuo_peor_ventana mediana/p90"] = (round(inf_ext["residuo_mediana_m"], 4), round(inf_ext["residuo_p90_m"], 4))

        def nuevos():
            return ({h["id"]: dict(Bt=np.zeros((2, h["perfil"].n_reb, nf), np.int64),
                                   Bi=np.zeros((2, h["perfil"].n_reb, nf), np.int64)) for h in hileras},
                    {h["id"]: dict(rayos_candidatos=0, rayos_validos=0, rayos_entran=0) for h in hileras})
        vias = {}
        if metodo in ("pcap", "ambos"):
            t0 = time.time()
            cont, est = nuevos()
            po.rayos_pcap(cfg, tray, pars, ventanas, hileras, cont, est, log=p.info)
            vias["pcap"] = (cont, est)
            p.info("vía pcap lista en %.0f s: %s", time.time() - t0, est)
            if n_val == len(pars):
                vias["pcap_todas"] = (cont, est)     # todas las ventanas son válidas -> la alternativa coincide con la oficial
            else:
                # Columnas ALTERNATIVAS (no oficiales): P con TODAS las ventanas, también las que no cumplen el criterio de la
                # extrínseca. Sirven para ver la magnitud cuando el criterio deja una hilera sin P oficial.
                t0 = time.time()
                cont2, est2 = nuevos()
                po.rayos_pcap(cfg, tray, [dict(pp, valida=True) for pp in pars], ventanas, hileras, cont2, est2, log=lambda *a_: None)
                vias["pcap_todas"] = (cont2, est2)
                p.info("vía pcap con todas las ventanas (alternativa, no oficial) lista en %.0f s", time.time() - t0)
        if metodo in ("rayos", "ambos"):
            t0 = time.time()
            cont, est = nuevos()
            po.rayos_las(cfg, tray, pars, lims, tl, xyz, hileras, cont, est)
            vias["rayos"] = (cont, est)
            p.info("vía rayos-LAS (solo contraste) lista en %.0f s: %s", time.time() - t0, est)
        # ---------------- tabla por rebanada
        rc = pc["rayos"]
        filas = []
        resumen_P = {}
        for h in hileras:
            P_, n_reb = h["perfil"], h["perfil"].n_reb
            sc = h["s_ini"] + (np.arange(n_reb) + 0.5) * m_["rebanada_m"]
            xy = h["marco"].xy(sc, np.zeros(n_reb))
            df = pd.DataFrame(dict(id_fila=h["id"], **h["atributos"], rebanada=np.arange(n_reb),
                                   s_ini=h["s_ini"] + np.arange(n_reb) * m_["rebanada_m"], s_centro=sc, x_eje=xy[:, 0], y_eje=xy[:, 1],
                                   n_pts=P_.n_pts, MH_max=P_.MH_max, MH_p=P_.MH_p, MW_max=P_.MW_max, MW_p=P_.MW_p,
                                   **h["copa"], dist_hileras_m=h["dist_hileras"], origen_dist_hileras=h["origen_dist"]))
            for via, (cont, est) in vias.items():
                c = cont[h["id"]]
                for modo, suf in (("suma_rayos", ""), ("media_franjas", "_mf")):
                    Pl, nfi, Pm = lc.p_por_rebanada(c["Bt"], c["Bi"], rc["min_rayos_celda"], pc["min_franjas_fiables"], modo)
                    df[f"P_izq_{via}{suf}"], df[f"P_der_{via}{suf}"], df[f"P_{via}{suf}"] = Pl[0], Pl[1], Pm
                # sensibilidad (QA-04): P solo con las franjas desde 0,6 m
                Pl6, nfi6, Pm6 = lc.p_por_rebanada(c["Bt"], c["Bi"], rc["min_rayos_celda"], pc["min_franjas_fiables"], "suma_rayos", j_min=j06)
                df[f"P_izq_{via}_h06"], df[f"P_der_{via}_h06"], df[f"P_{via}_h06"] = Pl6[0], Pl6[1], Pm6
                # cobertura vertical de P (QA-05)
                hm, fr, frm = lc.cobertura_P(c["Bt"], P_.MH_max, u, m_["franja_m"], rc["min_rayos_celda"])
                df[f"h_max_P_izq_{via}"], df[f"h_max_P_der_{via}"] = hm[0], hm[1]
                df[f"cob_MH_izq_{via}"], df[f"cob_MH_der_{via}"], df[f"cob_MH_{via}"] = fr[0], fr[1], frm
                df[f"Bt_izq_{via}"], df[f"Bt_der_{via}"] = c["Bt"][0].sum(1), c["Bt"][1].sum(1)
                df[f"Bi_izq_{via}"], df[f"Bi_der_{via}"] = c["Bi"][0].sum(1), c["Bi"][1].sum(1)
                df[f"franjas_fiables_izq_{via}"] = (c["Bt"][0] >= rc["min_rayos_celda"]).sum(1)
                df[f"franjas_fiables_der_{via}"] = (c["Bt"][1] >= rc["min_rayos_celda"]).sum(1)
            filas.append(df)
            resumen_P[h["id"]] = {via: dict(P_media=float(np.nanmean(df[f"P_{via}"])) if df[f"P_{via}"].notna().any() else None,
                                            P_media_h06=float(np.nanmean(df[f"P_{via}_h06"])) if df[f"P_{via}_h06"].notna().any() else None,
                                            h_max_P_mediana_izq=float(np.nanmedian(df[f"h_max_P_izq_{via}"])) if df[f"h_max_P_izq_{via}"].notna().any() else None,
                                            h_max_P_mediana_der=float(np.nanmedian(df[f"h_max_P_der_{via}"])) if df[f"h_max_P_der_{via}"].notna().any() else None,
                                            cob_MH_media=float(np.nanmean(df[f"cob_MH_{via}"])) if df[f"cob_MH_{via}"].notna().any() else None,
                                            n_reb_con_P=int(df[f"P_{via}"].notna().sum()), n_reb=int(n_reb), rayos=vias[via][1][h["id"]])
                                  for via in vias}
        tabla = pd.concat(filas, ignore_index=True)
        via_ppal = "pcap" if "pcap" in vias else "rayos"      # P OFICIAL = pcap con ventanas válidas; si falta, NaN con motivo (sin sustituciones silenciosas)
        if n_val < len(pars):
            p.info("AVISO: %d/%d ventanas no cumplen el criterio de la extrínseca: P oficial solo donde hay rayos de ventanas válidas "
                   "(ver P_mejor / via_P_mejor para la alternativa)", len(pars) - n_val, len(pars))
        tabla["P"] = tabla[f"P_{via_ppal}"]
        tabla["P_h06"] = tabla[f"P_{via_ppal}_h06"]
        for c_ in ("h_max_P_izq", "h_max_P_der", "cob_MH_izq", "cob_MH_der", "cob_MH"):
            tabla[c_] = tabla[f"{c_}_{via_ppal}"]
        tabla["via_P"] = via_ppal
        # mejor valor disponible (explícito en via_P_mejor): pcap oficial > pcap con ventana no válida > rayos-LAS (sesgo ~ -0,05)
        pm = tabla["P"].copy()
        vm = np.where(pm.notna(), "pcap", "")
        for alt, nom in (("pcap_todas", "pcap_ventana_no_valida"), ("rayos", "rayos_LAS")):
            if f"P_{alt}" in tabla:
                falta = pm.isna() & tabla[f"P_{alt}"].notna()
                pm[falta] = tabla.loc[falta, f"P_{alt}"]
                vm = np.where(falta, nom, vm)
        tabla["P_mejor"] = pm
        tabla["via_P_mejor"] = np.where(pm.notna(), vm, "sin_P")
        tabla["n_ventanas_validas"] = n_val
        # motivo de cada NaN en P (nunca silencioso)
        mot = np.full(len(tabla), "", dtype=object)
        ppal = tabla["P"].isna().values
        sin_mh = tabla["MH_max"].isna().values
        if via_ppal == "pcap" and not pcap_valida:
            mot[ppal] = "extrinseca no valida: 0/%d ventanas cumplen mediana<%.2f y p90<%.2f" % (len(pars), ce["residuo_mediana_max_m"], ce["residuo_p90_max_m"])
        else:
            izq_na = tabla[f"P_izq_{via_ppal}"].isna().values
            der_na = tabla[f"P_der_{via_ppal}"].isna().values
            mot[ppal & izq_na & ~der_na] = "lado izquierdo sin franjas fiables (<%d franjas con >=%d rayos)" % (pc["min_franjas_fiables"], rc["min_rayos_celda"])
            mot[ppal & der_na & ~izq_na] = "lado derecho sin franjas fiables (<%d franjas con >=%d rayos)" % (pc["min_franjas_fiables"], rc["min_rayos_celda"])
            mot[ppal & izq_na & der_na] = "ambos lados sin franjas fiables (sin rayos suficientes)"
            if n_val < len(pars):
                mot[ppal] = [(x + " | ") if x else "" for x in mot[ppal]]
                mot[ppal] = [x + "%d/%d ventanas excluidas por extrinseca" % (len(pars) - n_val, len(pars)) for x in mot[ppal]]
        mot[sin_mh & ppal] = "rebanada sin puntos de copa (MH/MW NaN)"
        tabla["P_motivo_nan"] = mot
        tabla.to_csv(os.path.join(out, "metricas_rebanadas.csv"), index=False)
        # ---------------- perfil de anchura por altura (hilera entera y tramos de hilera.tramo_m)
        perf = []
        for h in hileras:
            P_ = h["perfil"]
            tramo = np.floor(np.arange(P_.n_reb) * m_["rebanada_m"] / cfg["hilera"]["tramo_m"] + 1e-9).astype(int)
            for nivel, grupos in (("hilera", None), ("tramo", tramo)):
                for f_ in lc.perfil_anchura(P_, grupos):
                    g_ = f_.pop("grupo")
                    perf.append(dict(id_fila=h["id"], **h["atributos"], nivel=nivel,
                                     tramo=g_ if nivel == "tramo" else None, **f_))
            ph_ = [f_ for f_ in perf if f_["id_fila"] == h["id"] and f_["nivel"] == "hilera"]
            if ph_:
                hc = np.array([0.5 * (f_["h_ini_m"] + f_["h_fin_m"]) for f_ in ph_])
                png_lineas([(np.array([f_["anchura_media_m"] for f_ in ph_]), hc, "anchura estricta (franjas ocupadas)"),
                            (np.array([f_["anchura_media_p_m"] for f_ in ph_]), hc, "anchura robusta (p2,5-p97,5)"),
                            (np.array([f_["anchura_media_con_huecos_m"] for f_ in ph_]), hc, "anchura contando huecos")],
                           os.path.join(out, f"perfil_anchura_hilera{h['id']}.png"), f"Hilera {h['id']}: anchura por altura",
                           "anchura (m)", "altura sobre el suelo (m)")
        pd.DataFrame(perf).to_csv(os.path.join(out, "perfil_anchura.csv"), index=False)
        np.savez_compressed(os.path.join(out, "conteos_rayos.npz"),
                            **{f"{via}__{k}__{n}": a for via, (cont, _) in vias.items() for k, c in cont.items() for n, a in c.items()})
        # ---------------- comparación pcap vs rayos
        comp = dict(via_principal=via_ppal, por_hilera=resumen_P, extrinseca=inf_ext,
                    umbral_oficial_m=u, umbral_sensibilidad_m=u06)
        if "pcap" in vias and "rayos" in vias:
            a, b = tabla["P_pcap"].values, tabla["P_rayos"].values
            ok2 = np.isfinite(a) & np.isfinite(b)
            comp["P_pcap_vs_rayos"] = dict(n=int(ok2.sum()), media_pcap=float(a[ok2].mean()) if ok2.any() else None,
                                           media_rayos=float(b[ok2].mean()) if ok2.any() else None,
                                           corr=float(np.corrcoef(a[ok2], b[ok2])[0, 1]) if ok2.sum() > 3 else None,
                                           dif_media=float((a[ok2] - b[ok2]).mean()) if ok2.any() else None)
        guardar_json(comp, os.path.join(out, "comparacion_P.json"))
        p.info("comparación P: %s", comp.get("P_pcap_vs_rayos"))
        for k in ids:
            sub = tabla[tabla.id_fila == k]
            series = [(sub.s_centro, sub[f"P_{v}"], f"P {v}") for v in vias] + [(sub.s_centro, sub.P_h06, f"P {via_ppal} desde {u06} m"),
                                                                                 (sub.s_centro, sub.MH_max / 4, "MH/4"), (sub.s_centro, sub.MW_max / 4, "MW/4")]
            png_lineas(series, os.path.join(out, f"perfil_P_hilera{k}.png"), f"Hilera {k}: P por rebanada", "s (m)", "P (0-1)")
        p.n_despues = int(len(tabla))
        p.salida = os.path.join(out, "metricas_rebanadas.csv")
    return tabla


if __name__ == "__main__":
    main()
